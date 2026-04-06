/* cfx_bge_crypto.c -- crypto operations for BGE stores */

#include "cfx_bge_internal.h"

/*  wipe helpers  */

void bge_v4_store_wipe(bge_v4_store *s) {
    cfx_memzero_s(s->dek, sizeof(s->dek));
    cfx_memzero_s(&s->hdr, sizeof(s->hdr));
    cfx_memzero_s(s->slots, sizeof(s->slots));
    cfx_memzero_s(s->data_nonce, sizeof(s->data_nonce));
    if (s->file_buf) {
        cfx_memzero_s(s->file_buf, s->file_len);
        free(s->file_buf);
        s->file_buf = NULL;
    }
    s->file_len = 0;
    s->matched_slot = -1;
}

void bge_ustore_wipe(bge_ustore *us) {
    bge_v4_store_wipe(&us->v4);
}

/*  v4 slot serialization (field-by-field, no struct packing)  */

void bge_v4_slot_from_buf(bge_v4_slot *slot, const uint8_t buf[BGE_STORE_SLOT_LEN]) {
    const uint8_t *p = buf;
    slot->m_cost = cfx_load32_le(p); p += 4;
    slot->t_cost = cfx_load32_le(p); p += 4;
    slot->p_cost = cfx_load32_le(p); p += 4;
    memcpy(slot->salt,        p, 16); p += 16;
    memcpy(slot->slot_nonce,  p, 24); p += 24;
    memcpy(slot->verifier,    p, 16); p += 16;
    memcpy(slot->wrapped_dek, p, 32); p += 32;
    memcpy(slot->wrap_tag,    p, 16);
}

void bge_v4_slot_to_buf(uint8_t buf[BGE_STORE_SLOT_LEN], const bge_v4_slot *slot) {
    uint8_t *p = buf;
    cfx_store32_le(p, slot->m_cost); p += 4;
    cfx_store32_le(p, slot->t_cost); p += 4;
    cfx_store32_le(p, slot->p_cost); p += 4;
    memcpy(p, slot->salt,        16); p += 16;
    memcpy(p, slot->slot_nonce,  24); p += 24;
    memcpy(p, slot->verifier,    16); p += 16;
    memcpy(p, slot->wrapped_dek, 32); p += 32;
    memcpy(p, slot->wrap_tag,    16);
}

/* derive KEK from password via Argon2, wrap DEK into slot */
int bge_v4_wrap_dek(const char *pwd, size_t pwd_len,
                    uint32_t m, uint32_t t, uint32_t p,
                    const uint8_t dek[BGE_STORE_DEK_LEN], bge_v4_slot *slot) {
    /* Argon2 -> 48 bytes: KEK[32] | verifier[16] */
    uint8_t kdf_out[48];

    cfx_srand_os();
    cfx_rand_bytes(slot->salt, sizeof(slot->salt));
    cfx_rand_bytes(slot->slot_nonce, sizeof(slot->slot_nonce));

    slot->m_cost = m;
    slot->t_cost = t;
    slot->p_cost = p;


    int rc = cfx_argon2id(kdf_out, 48,
        (const uint8_t *)pwd, pwd_len,
        slot->salt, sizeof(slot->salt), m, t, p);
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        return -1;
    }

    memcpy(slot->verifier, kdf_out + 32, 16);

    /* wrap DEK with KEK via XChaCha20-Poly1305 (no AAD) */
    rc = cfx_xchacha20_poly1305_encrypt(
        slot->wrapped_dek, slot->wrap_tag,
        dek, BGE_STORE_DEK_LEN,
        NULL, 0,
        kdf_out, slot->slot_nonce);

    cfx_memzero_s(kdf_out, sizeof(kdf_out));
    if (rc != 0) {
        fprintf(stderr, "error: DEK wrapping failed\n");
        return -1;
    }
    return 0;
}

/* serialize header + slots + data_nonce into AAD buffer.
 * returns total AAD length. caller provides buf of sufficient size. */
static size_t bge_v4_build_aad(uint8_t *aad,
                                const bge_v4_fixed_header *hdr,
                                const bge_v4_slot *slots, int slot_count,
                                const uint8_t data_nonce[24]) {
    uint8_t *w = aad;
    memcpy(w, hdr, BGE_STORE_HDR_LEN);
    w += BGE_STORE_HDR_LEN;

    for (int i = 0; i < slot_count; i++) {
        bge_v4_slot_to_buf(w, &slots[i]);
        w += BGE_STORE_SLOT_LEN;
    }

    memcpy(w, data_nonce, 24);
    w += 24;

    return (size_t)(w - aad);
}

/* atomic write for v4: same .tmp -> fsync -> rename pattern */
int bge_v4_safe_write(const char *path,
                      const bge_v4_fixed_header *hdr,
                      const bge_v4_slot *slots, int slot_count,
                      const uint8_t data_nonce[24],
                      const uint8_t *ct, size_t ct_len,
                      const uint8_t tag[BGE_TAG_LEN]) {
    char tmppath[1088];
    int n = snprintf(tmppath, sizeof(tmppath), "%s.tmp", path);
    if (n < 0 || (size_t)n >= sizeof(tmppath)) {
        fprintf(stderr, "error: path too long\n");
        return -1;
    }

    int lock_fd = -1;
    int ret = -1;

#ifndef _WIN32
    lock_fd = open(path, O_RDONLY);
    if (lock_fd >= 0) {
        if (flock(lock_fd, LOCK_EX) != 0) {
            fprintf(stderr, "error: cannot lock %s: %s\n", path, strerror(errno));
            close(lock_fd);
            return -1;
        }
    }

    int fd = open(tmppath, O_CREAT | O_WRONLY | O_TRUNC, 0600);
    if (fd < 0) {
        fprintf(stderr, "error: cannot create %s: %s\n", tmppath, strerror(errno));
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }
    FILE *f = fdopen(fd, "wb");
    if (!f) {
        close(fd);
        unlink(tmppath);
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }
#else
    FILE *f = fopen(tmppath, "wb");
    if (!f) {
        fprintf(stderr, "error: cannot create %s: %s\n", tmppath, strerror(errno));
        return -1;
    }
#endif

    size_t ok = 1;

    /* fixed header */
    ok = ok && fwrite(hdr, 1, BGE_STORE_HDR_LEN, f) == BGE_STORE_HDR_LEN;

    /* slots (serialized field-by-field) */
    for (int i = 0; i < slot_count && ok; i++) {
        uint8_t sbuf[BGE_STORE_SLOT_LEN];
        bge_v4_slot_to_buf(sbuf, &slots[i]);
        ok = ok && fwrite(sbuf, 1, BGE_STORE_SLOT_LEN, f) == BGE_STORE_SLOT_LEN;
    }

    /* data nonce */
    ok = ok && fwrite(data_nonce, 1, 24, f) == 24;

    /* ciphertext */
    if (ct_len > 0)
        ok = ok && fwrite(ct, 1, ct_len, f) == ct_len;

    /* tag */
    ok = ok && fwrite(tag, 1, BGE_TAG_LEN, f) == BGE_TAG_LEN;

    if (!ok) {
        fprintf(stderr, "error: write failed: %s\n", strerror(errno));
        fclose(f);
        unlink(tmppath);
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }

#ifndef _WIN32
    fsync(fileno(f));
#else
    fflush(f);
#endif
    fclose(f);

#ifdef _WIN32
    if (!MoveFileExA(tmppath, path, MOVEFILE_REPLACE_EXISTING)) {
        fprintf(stderr, "error: rename failed: %lu\n", GetLastError());
#else
    if (rename(tmppath, path) != 0) {
        fprintf(stderr, "error: rename failed: %s\n", strerror(errno));
#endif
        unlink(tmppath);
        if (lock_fd >= 0) close(lock_fd);
        return -1;
    }

    ret = 0;
    if (lock_fd >= 0) close(lock_fd);
    return ret;
}

int bge_v4_encrypt_write(const char *path,
                         const uint8_t *pt, size_t pt_len,
                         const bge_v4_fixed_header *hdr,
                         const bge_v4_slot *slots, int slot_count,
                         const uint8_t dek[BGE_STORE_DEK_LEN]) {
    int ret = -1;
    uint8_t data_nonce[24];

    cfx_srand_os();
    cfx_rand_bytes(data_nonce, sizeof(data_nonce));

    /* build AAD = header + all slots + data_nonce */
    size_t aad_len = BGE_STORE_HDR_LEN + (size_t)slot_count * BGE_STORE_SLOT_LEN + 24;
    uint8_t *aad = malloc(aad_len);
    if (!aad) {
        fprintf(stderr, "error: allocation failed\n");
        return -1;
    }
    bge_v4_build_aad(aad, hdr, slots, slot_count, data_nonce);

    uint8_t *ct = malloc(pt_len);
    uint8_t tag[BGE_TAG_LEN];
    if (pt_len > 0 && !ct) {
        fprintf(stderr, "error: allocation failed\n");
        free(aad);
        return -1;
    }

    int rc = cfx_xchacha20_poly1305_encrypt(
        ct, tag, pt, pt_len,
        aad, aad_len,
        dek, data_nonce);
    if (rc != 0) {
        fprintf(stderr, "error: encryption failed\n");
        free(ct);
        free(aad);
        return -1;
    }

    ret = bge_v4_safe_write(path, hdr, slots, slot_count, data_nonce, ct, pt_len, tag);

    cfx_memzero_s(ct, pt_len);
    free(ct);
    free(aad);
    return ret;
}

/* create a v4 store with 1 slot and write to disk */
int bge_v4_init_write(const char *path, const uint8_t *pt, size_t pt_len,
                      const char *pwd, size_t pwd_len,
                      uint32_t m, uint32_t t, uint32_t p) {
    uint8_t dek[BGE_STORE_DEK_LEN];
    cfx_srand_os();
    cfx_rand_bytes(dek, sizeof(dek));

    bge_v4_fixed_header hdr;
    memcpy(hdr.magic, BGE_MAGIC, 3);
    hdr.version = BGE_STORE_VERSION;
    hdr.slot_count = 1;
    memset(hdr.reserved, 0, sizeof(hdr.reserved));

    bge_v4_slot slot;
    memset(&slot, 0, sizeof(slot));
    int rc = bge_v4_wrap_dek(pwd, pwd_len, m, t, p, dek, &slot);
    if (rc != 0) {
        cfx_memzero_s(dek, sizeof(dek));
        return -1;
    }

    rc = bge_v4_encrypt_write(path, pt, pt_len, &hdr, &slot, 1, dek);
    cfx_memzero_s(dek, sizeof(dek));
    cfx_memzero_s(&slot, sizeof(slot));
    return rc;
}

/* create a v4 store with 1 slot in a caller-provided buffer (for tests) */
int bge_v4_init_to_buf(uint8_t **out, size_t *out_len,
                       const uint8_t *pt, size_t pt_len,
                       const char *pwd, size_t pwd_len,
                       uint32_t m, uint32_t t, uint32_t p) {
    uint8_t dek[BGE_STORE_DEK_LEN];
    cfx_srand_os();
    cfx_rand_bytes(dek, sizeof(dek));

    bge_v4_fixed_header hdr;
    memcpy(hdr.magic, BGE_MAGIC, 3);
    hdr.version = BGE_STORE_VERSION;
    hdr.slot_count = 1;
    memset(hdr.reserved, 0, sizeof(hdr.reserved));

    bge_v4_slot slot;
    memset(&slot, 0, sizeof(slot));
    int rc = bge_v4_wrap_dek(pwd, pwd_len, m, t, p, dek, &slot);
    if (rc != 0) {
        cfx_memzero_s(dek, sizeof(dek));
        return -1;
    }

    uint8_t data_nonce[24];
    cfx_rand_bytes(data_nonce, sizeof(data_nonce));

    size_t aad_len = BGE_STORE_HDR_LEN + BGE_STORE_SLOT_LEN + 24;
    uint8_t *aad = malloc(aad_len);
    if (!aad) {
        cfx_memzero_s(dek, sizeof(dek));
        return -1;
    }
    bge_v4_build_aad(aad, &hdr, &slot, 1, data_nonce);

    uint8_t *ct = malloc(pt_len > 0 ? pt_len : 1);
    uint8_t tag[BGE_TAG_LEN];
    if (!ct) {
        free(aad);
        cfx_memzero_s(dek, sizeof(dek));
        return -1;
    }

    rc = cfx_xchacha20_poly1305_encrypt(
        ct, tag, pt, pt_len, aad, aad_len, dek, data_nonce);
    cfx_memzero_s(dek, sizeof(dek));
    if (rc != 0) {
        free(ct);
        free(aad);
        return -1;
    }

    /* assemble: fixed_header | slot_buf | data_nonce | ct | tag */
    size_t blob_len = BGE_STORE_HDR_LEN + BGE_STORE_SLOT_LEN + 24 + pt_len + BGE_TAG_LEN;
    uint8_t *blob = malloc(blob_len);
    if (!blob) {
        free(ct);
        free(aad);
        return -1;
    }

    uint8_t *w = blob;
    memcpy(w, &hdr, BGE_STORE_HDR_LEN); w += BGE_STORE_HDR_LEN;
    uint8_t slot_buf[BGE_STORE_SLOT_LEN];
    bge_v4_slot_to_buf(slot_buf, &slot);
    memcpy(w, slot_buf, BGE_STORE_SLOT_LEN); w += BGE_STORE_SLOT_LEN;
    memcpy(w, data_nonce, 24); w += 24;
    if (pt_len > 0) { memcpy(w, ct, pt_len); w += pt_len; }
    memcpy(w, tag, BGE_TAG_LEN);

    cfx_memzero_s(&slot, sizeof(slot));
    free(ct);
    free(aad);

    *out = blob;
    *out_len = blob_len;
    return 0;
}

/* try each slot: Argon2 -> check verifier -> unwrap DEK */
int bge_v4_authenticate_buf(const uint8_t *file_buf, size_t file_len,
                            const char *pwd, size_t pwd_len,
                            bge_v4_store *store) {
    if (file_len < BGE_STORE_HDR_LEN) {
        fprintf(stderr, "error: file too small for v4 header\n");
        return -1;
    }

    /* parse fixed header */
    bge_v4_fixed_header hdr;
    memcpy(&hdr, file_buf, BGE_STORE_HDR_LEN);

    if (memcmp(hdr.magic, BGE_MAGIC, 3) != 0 || hdr.version != BGE_STORE_VERSION) {
        fprintf(stderr, "error: not a valid BGE v4 store\n");
        return -1;
    }

    int slot_count = hdr.slot_count;
    if (slot_count < 1 || slot_count > BGE_STORE_MAX_SLOTS) {
        fprintf(stderr, "error: invalid slot count %d\n", slot_count);
        return -1;
    }

    size_t prefix_len = BGE_STORE_HDR_LEN + (size_t)slot_count * BGE_STORE_SLOT_LEN + 24;
    if (file_len < prefix_len + BGE_TAG_LEN) {
        fprintf(stderr, "error: file too small for v4 data\n");
        return -1;
    }

    /* parse slots */
    bge_v4_slot slots[BGE_STORE_MAX_SLOTS];
    const uint8_t *sp = file_buf + BGE_STORE_HDR_LEN;
    for (int i = 0; i < slot_count; i++) {
        bge_v4_slot_from_buf(&slots[i], sp);
        sp += BGE_STORE_SLOT_LEN;
    }

    /* data nonce follows slots */
    uint8_t data_nonce[24];
    memcpy(data_nonce, sp, 24);

    /* try each slot */
    int matched = -1;
    uint8_t dek[BGE_STORE_DEK_LEN];

    for (int i = 0; i < slot_count; i++) {
        uint32_t m = slots[i].m_cost;
        uint32_t t = slots[i].t_cost;
        uint32_t p = slots[i].p_cost;

        if (m < 8 || t < 1 || p < 1 ||
            m > 4194304 || t > 100 || p > 16) {
            continue;  /* skip unreasonable slot */
        }

        uint8_t kdf_out[48];
    
        int rc = cfx_argon2id(kdf_out, 48,
            (const uint8_t *)pwd, pwd_len,
            slots[i].salt, sizeof(slots[i].salt), m, t, p);
        if (rc != 0) {
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            continue;
        }

        /* constant-time verifier check */
        uint8_t diff = 0;
        for (int j = 0; j < 16; j++)
            diff |= kdf_out[32 + j] ^ slots[i].verifier[j];

        if (diff != 0) {
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            continue;
        }

        /* verifier matched -- unwrap DEK */
        rc = cfx_xchacha20_poly1305_decrypt(
            dek, slots[i].wrapped_dek, BGE_STORE_DEK_LEN,
            NULL, 0,
            kdf_out, slots[i].slot_nonce, slots[i].wrap_tag);

        cfx_memzero_s(kdf_out, sizeof(kdf_out));

        if (rc != 0) {
            /* verifier matched but unwrap failed -- corrupted slot */
            fprintf(stderr, "error: slot %d verifier matched but DEK unwrap failed\n", i);
            cfx_memzero_s(dek, sizeof(dek));
            return -1;
        }

        matched = i;
        break;
    }

    if (matched < 0) {
        fprintf(stderr, "error: wrong passphrase\n");
        return -1;
    }

    /* fill store */
    store->hdr = hdr;
    memcpy(store->slots, slots, sizeof(slots));
    memcpy(store->data_nonce, data_nonce, 24);
    memcpy(store->dek, dek, BGE_STORE_DEK_LEN);
    store->matched_slot = matched;

    store->file_buf = malloc(file_len);
    if (!store->file_buf) {
        cfx_memzero_s(dek, sizeof(dek));
        cfx_memzero_s(store, sizeof(*store));
        return -1;
    }
    memcpy(store->file_buf, file_buf, file_len);
    store->file_len = file_len;

    cfx_memzero_s(dek, sizeof(dek));
    return 0;
}

int bge_v4_decrypt_store(const bge_v4_store *store,
                         uint8_t **pt_out, size_t *pt_len) {
    int slot_count = store->hdr.slot_count;
    size_t prefix_len = BGE_STORE_HDR_LEN + (size_t)slot_count * BGE_STORE_SLOT_LEN + 24;
    size_t ct_len = store->file_len - prefix_len - BGE_TAG_LEN;
    const uint8_t *ct  = store->file_buf + prefix_len;
    const uint8_t *tag = store->file_buf + store->file_len - BGE_TAG_LEN;

    /* rebuild AAD from stored data */
    size_t aad_len = prefix_len;
    uint8_t *aad = malloc(aad_len);
    if (!aad) {
        fprintf(stderr, "error: allocation failed\n");
        return -1;
    }
    bge_v4_build_aad(aad, &store->hdr, store->slots, slot_count, store->data_nonce);

    uint8_t *plaintext = malloc(ct_len + 1);
    if (!plaintext) {
        fprintf(stderr, "error: allocation failed\n");
        free(aad);
        return -1;
    }

    int rc = cfx_xchacha20_poly1305_decrypt(
        plaintext, ct, ct_len,
        aad, aad_len,
        store->dek, store->data_nonce, tag);
    free(aad);

    if (rc != 0) {
        fprintf(stderr, "error: store corrupted or tampered\n");
        cfx_memzero_s(plaintext, ct_len + 1);
        free(plaintext);
        return -1;
    }

    plaintext[ct_len] = '\0';
    *pt_out = plaintext;
    *pt_len = ct_len;
    return 0;
}

/*  unified dispatch  */

int bge_uauthenticate(const char *path, const char *pwd, size_t pwd_len,
                      bge_ustore *us) {
    FILE *f = fopen(path, "rb");
    if (!f) {
        fprintf(stderr, "error: cannot open %s: %s\n", path, strerror(errno));
        return -1;
    }

    uint8_t *file_buf = NULL;
    size_t file_len = 0;
    if (bge_read_all(f, &file_buf, &file_len) != 0) {
        fclose(f);
        fprintf(stderr, "error: cannot read %s\n", path);
        return -1;
    }
    fclose(f);

    if (file_len < 4 || memcmp(file_buf, BGE_MAGIC, 3) != 0) {
        fprintf(stderr, "error: %s is not a valid BGE store\n", path);
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    if (file_buf[3] != BGE_STORE_VERSION) {
        fprintf(stderr, "error: unsupported BGE version %u\n", file_buf[3]);
        cfx_memzero_s(file_buf, file_len);
        free(file_buf);
        return -1;
    }

    int rc = bge_v4_authenticate_buf(file_buf, file_len, pwd, pwd_len, &us->v4);
    cfx_memzero_s(file_buf, file_len);
    free(file_buf);
    return rc;
}

int bge_udecrypt(const bge_ustore *us, uint8_t **pt_out, size_t *pt_len) {
    return bge_v4_decrypt_store(&us->v4, pt_out, pt_len);
}

int bge_uwrite(const char *path, const uint8_t *pt, size_t pt_len,
               const bge_ustore *us) {
    const bge_v4_store *s = &us->v4;
    return bge_v4_encrypt_write(path, pt, pt_len,
                                &s->hdr, s->slots, s->hdr.slot_count,
                                s->dek);
}
