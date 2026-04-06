/* cfx_bge_file.c -- streaming file encryption/decryption for BGE */

#include "cfx_bge_internal.h"

/*  armor helpers  */

int bge_is_armored(const uint8_t *buf, size_t len) {
    return len >= 28 && memcmp(buf, "-----BEGIN BGE MESSAGE-----", 27) == 0;
}

/* wrap binary blob in PEM-style armor with 76-char lines. caller frees *out. */
int bge_armor_encode(const uint8_t *bin, size_t bin_len,
                     uint8_t **out, size_t *out_len) {
    size_t b64_len = 0;
    cfx_base64_encode(NULL, &b64_len, bin, bin_len);

    char *b64 = malloc(b64_len);
    if (!b64) return -1;
    cfx_base64_encode(b64, &b64_len, bin, bin_len);

    /* header + base64 + newlines every 76 chars + footer */
    size_t nlines = (b64_len + 75) / 76;
    size_t total = strlen(BGE_ARMOR_BEGIN) + b64_len + nlines + strlen(BGE_ARMOR_END);
    uint8_t *buf = malloc(total + 1);
    if (!buf) { free(b64); return -1; }

    uint8_t *w = buf;
    size_t hlen = strlen(BGE_ARMOR_BEGIN);
    memcpy(w, BGE_ARMOR_BEGIN, hlen); w += hlen;

    for (size_t i = 0; i < b64_len; i += 76) {
        size_t chunk = b64_len - i;
        if (chunk > 76) chunk = 76;
        memcpy(w, b64 + i, chunk); w += chunk;
        *w++ = '\n';
    }

    size_t flen = strlen(BGE_ARMOR_END);
    memcpy(w, BGE_ARMOR_END, flen); w += flen;

    free(b64);
    *out = buf;
    *out_len = (size_t)(w - buf);
    return 0;
}

/* strip PEM header/footer and base64-decode. caller frees *out. */
int bge_armor_decode(const uint8_t *text, size_t text_len,
                     uint8_t **out, size_t *out_len) {
    /* find start after first newline past header */
    const char *s = (const char *)text;
    const char *body = memchr(s, '\n', text_len);
    if (!body) return -1;
    body++;

    /* find footer */
    const char *footer = strstr(body, "-----END BGE MESSAGE-----");
    if (!footer) return -1;

    size_t b64_len = (size_t)(footer - body);

    /* decode */
    size_t dec_len = cfx_base64_dec_max_len(b64_len);
    uint8_t *dec = malloc(dec_len);
    if (!dec) return -1;

    int rc = cfx_base64_decode(dec, &dec_len, body, b64_len);
    if (rc != 0) { free(dec); return -1; }

    *out = dec;
    *out_len = dec_len;
    return 0;
}

/*  v2 single-shot decrypt from in-memory buffer ─ */
static int bge_decrypt_v2(const uint8_t *file_buf, size_t file_len,
                          const uint8_t key[32], const uint8_t nonce[24],
                          FILE *outf) {
    size_t ct_len = file_len - BGE_AAD_LEN - BGE_TAG_LEN;
    const uint8_t *ct  = file_buf + BGE_AAD_LEN;
    const uint8_t *tag = file_buf + file_len - BGE_TAG_LEN;

    uint8_t *pt = malloc(ct_len + 1);
    if (!pt) {
        fprintf(stderr, "error: allocation failed\n");
        return -1;
    }

    int rc = cfx_xchacha20_poly1305_decrypt(
        pt, ct, ct_len, file_buf, BGE_AAD_LEN, key, nonce, tag);
    if (rc != 0) {
        fprintf(stderr, "error: decryption failed (corrupted or tampered)\n");
        cfx_memzero_s(pt, ct_len + 1); free(pt);
        return -1;
    }

    int ret = 0;
    if (ct_len > 0 && fwrite(pt, 1, ct_len, outf) != ct_len)
        ret = -1;

    cfx_memzero_s(pt, ct_len + 1); free(pt);
    return ret;
}

/*  v3 decrypt from in-memory buffer (armored v3 input)  */

static int bge_decrypt_v3(const uint8_t *file_buf, size_t file_len,
                          const uint8_t key[32], const uint8_t nonce[24],
                          FILE *outf) {
    const uint8_t *p = file_buf + BGE_AAD_LEN;
    const uint8_t *end = file_buf + file_len;
    uint8_t pt_buf[CFX_STREAM_CHUNK_SIZE];
    uint64_t chunk_counter = 0;

    while (p < end) {
        size_t remaining = (size_t)(end - p);
        if (remaining < CFX_STREAM_TAG_SIZE) {
            fprintf(stderr, "error: truncated stream (no tag)\n");
            return -1;
        }

        size_t chunk_plus_tag;
        int is_final;

        if (remaining <= CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE) {
            chunk_plus_tag = remaining;
            is_final = 1;
        } else {
            chunk_plus_tag = CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE;
            is_final = 0;
        }

        size_t ct_len = chunk_plus_tag - CFX_STREAM_TAG_SIZE;
        const uint8_t *ct  = p;
        const uint8_t *tag = p + ct_len;

        int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
            pt_buf, ct, ct_len, tag, chunk_counter, is_final, key, nonce);
        if (rc != 0) {
            fprintf(stderr, "error: chunk %llu authentication failed\n",
                    (unsigned long long)chunk_counter);
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            return -1;
        }

        if (ct_len > 0 && fwrite(pt_buf, 1, ct_len, outf) != ct_len) {
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            return -1;
        }

        cfx_memzero_s(pt_buf, ct_len);
        p += chunk_plus_tag;
        chunk_counter++;
    }

    return 0;
}

/*  v3 streaming decrypt directly from FILE* ─ */
static int bge_decrypt_v3_stream(FILE *inf,
                                  const uint8_t key[32], const uint8_t nonce[24],
                                  FILE *outf) {
    uint8_t chunk_buf[CFX_STREAM_CHUNK_SIZE + CFX_STREAM_TAG_SIZE];
    uint8_t pt_buf[CFX_STREAM_CHUNK_SIZE];
    uint64_t chunk_counter = 0;
    uint8_t lookahead;
    int have_lookahead = 0;

    for (;;) {
        size_t off = 0;
        if (have_lookahead) {
            chunk_buf[0] = lookahead;
            off = 1;
            have_lookahead = 0;
        }

        size_t nread = fread(chunk_buf + off, 1, sizeof(chunk_buf) - off, inf);
        size_t total = off + nread;

        if (nread == 0 && off == 0 && ferror(inf)) {
            fprintf(stderr, "error: read failed\n");
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        if (total == 0) {
            fprintf(stderr, "error: empty stream (missing final chunk)\n");
            return -1;
        }

        if (total < CFX_STREAM_TAG_SIZE) {
            fprintf(stderr, "error: truncated stream (no tag)\n");
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        int is_final;
        if (total < sizeof(chunk_buf)) {
            is_final = 1;
        } else {
            size_t peeked = fread(&lookahead, 1, 1, inf);
            if (peeked == 0) {
                is_final = 1;
            } else {
                is_final = 0;
                have_lookahead = 1;
            }
        }

        size_t ct_len = total - CFX_STREAM_TAG_SIZE;
        const uint8_t *tag = chunk_buf + ct_len;

        int rc = cfx_stream_xchacha20_poly1305_decrypt_chunk(
            pt_buf, chunk_buf, ct_len, tag, chunk_counter, is_final, key, nonce);
        if (rc != 0) {
            fprintf(stderr, "error: chunk %llu authentication failed\n",
                    (unsigned long long)chunk_counter);
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        if (ct_len > 0 && fwrite(pt_buf, 1, ct_len, outf) != ct_len) {
            cfx_memzero_s(pt_buf, sizeof(pt_buf));
            cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
            return -1;
        }

        cfx_memzero_s(pt_buf, ct_len);
        chunk_counter++;

        if (is_final) break;
    }

    cfx_memzero_s(chunk_buf, sizeof(chunk_buf));
    return 0;
}

/*  encrypt file (v3 streaming)  */

int bge_encrypt_file(int argc, char **argv) {
    const char *output = NULL;
    const char *input = NULL;
    int armor = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-e") == 0 || strcmp(argv[i], "--encrypt") == 0) {
            continue;
        } else if (strcmp(argv[i], "-a") == 0 || strcmp(argv[i], "--armor") == 0) {
            armor = 1;
        } else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -o requires a path\n"); return 1; }
            output = argv[++i];
        } else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -i requires a path\n"); return 1; }
            input = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "error: unknown option: %s\n", argv[i]);
            return 1;
        } else if (!input) {
            input = argv[i];
        } else {
            fprintf(stderr, "error: unexpected argument: %s\n", argv[i]);
            return 1;
        }
    }

    FILE *inf = stdin;
    if (input) {
        inf = fopen(input, "rb");
        if (!inf) {
            fprintf(stderr, "error: cannot open %s: %s\n", input, strerror(errno));
            return 1;
        }
    }

    char pwd[256] = {0};
    int pwd_len = prompt_passphrase(pwd, sizeof(pwd));
    if (pwd_len < 0) {
        if (input) fclose(inf);
        return 1;
    }

    /* build v3 header */
    bge_header header;
    uint8_t kdf_out[48], verifier[BGE_VERIFIER_LEN];

    memcpy(header.magic, BGE_MAGIC, 3);
    header.version = BGE_STREAM_VERSION;
    cfx_store32_le(&header.m_cost, BGE_DEFAULT_M);
    cfx_store32_le(&header.t_cost, BGE_DEFAULT_T);
    cfx_store32_le(&header.p_cost, BGE_DEFAULT_P);

    cfx_srand_os();
    cfx_rand_bytes(header.salt,  sizeof(header.salt));
    cfx_rand_bytes(header.nonce, sizeof(header.nonce));

    int rc = cfx_argon2id(kdf_out, 48,
        (const uint8_t *)pwd, (size_t)pwd_len,
        header.salt, sizeof(header.salt),
        BGE_DEFAULT_M, BGE_DEFAULT_T, BGE_DEFAULT_P);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        if (input) fclose(inf);
        return 1;
    }
    memcpy(verifier, kdf_out + 32, BGE_VERIFIER_LEN);

    /* open output */
    FILE *outf = stdout;
    if (output) {
        outf = fopen(output, "wb");
        if (!outf) {
            fprintf(stderr, "error: cannot open %s: %s\n", output, strerror(errno));
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            if (input) fclose(inf);
            return 1;
        }
    }

    /* for armor: accumulate all output in memory, then encode at the end */
    uint8_t *armor_buf = NULL;
    size_t armor_cap = 0, armor_len = 0;
    int ret = 0;

    #define EMIT(data, len) do { \
        if (armor) { \
            size_t _n = (len); \
            if (armor_len + _n > armor_cap) { \
                size_t _nc = armor_cap ? armor_cap * 2 : 4096; \
                while (_nc < armor_len + _n) _nc *= 2; \
                uint8_t *_t = realloc(armor_buf, _nc); \
                if (!_t) { ret = 1; goto done; } \
                armor_buf = _t; armor_cap = _nc; \
            } \
            memcpy(armor_buf + armor_len, (data), _n); \
            armor_len += _n; \
        } else { \
            if (fwrite((data), 1, (len), outf) != (len)) { ret = 1; goto done; } \
        } \
    } while(0)

    /* write header + verifier */
    EMIT(&header, sizeof(header));
    EMIT(verifier, BGE_VERIFIER_LEN);

    /* streaming encrypt loop */
    uint8_t pt_buf[CFX_STREAM_CHUNK_SIZE];
    uint8_t ct_buf[CFX_STREAM_CHUNK_SIZE];
    uint8_t tag[CFX_STREAM_TAG_SIZE];
    uint64_t chunk_counter = 0;

    for (;;) {
        size_t nread = fread(pt_buf, 1, CFX_STREAM_CHUNK_SIZE, inf);
        if (nread == 0 && ferror(inf)) {
            fprintf(stderr, "error: read failed\n");
            ret = 1; goto done;
        }

        int is_final = feof(inf);
        if (!is_final && nread < CFX_STREAM_CHUNK_SIZE) {
            is_final = 1;
        }

        rc = cfx_stream_xchacha20_poly1305_encrypt_chunk(
            ct_buf, tag, pt_buf, nread,
            chunk_counter, is_final, kdf_out, header.nonce);
        cfx_memzero_s(pt_buf, sizeof(pt_buf));
        if (rc != 0) {
            fprintf(stderr, "error: encryption failed at chunk %llu\n",
                    (unsigned long long)chunk_counter);
            ret = 1;
            goto done;
        }

        EMIT(ct_buf, nread);
        EMIT(tag, CFX_STREAM_TAG_SIZE);
        chunk_counter++;

        if (is_final) break;
    }

    #undef EMIT

    /* finalize armor if needed */
    if (armor && ret == 0) {
        uint8_t *armored = NULL;
        size_t armored_len = 0;
        if (bge_armor_encode(armor_buf, armor_len, &armored, &armored_len) != 0) {
            fprintf(stderr, "error: armor encoding failed\n");
            ret = 1;
        } else {
            if (fwrite(armored, 1, armored_len, outf) != armored_len) ret = 1;
            free(armored);
        }
    }

done:
    cfx_memzero_s(kdf_out, sizeof(kdf_out));
    cfx_memzero_s(verifier, sizeof(verifier));
    cfx_memzero_s(pt_buf, sizeof(pt_buf));
    cfx_memzero_s(ct_buf, sizeof(ct_buf));
    if (armor_buf) { free(armor_buf); }
    if (input) fclose(inf);
    if (output && outf) fclose(outf);
    return ret;
}

int bge_decrypt_file(int argc, char **argv) {
    const char *output = NULL;
    const char *input = NULL;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--decrypt") == 0)
            continue;
        else if (strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "--output") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -o requires a path\n"); return 1; }
            output = argv[++i];
        } else if (strcmp(argv[i], "-i") == 0 || strcmp(argv[i], "--input") == 0) {
            if (i + 1 >= argc) { fprintf(stderr, "error: -i requires a path\n"); return 1; }
            input = argv[++i];
        } else if (argv[i][0] == '-') {
            fprintf(stderr, "error: unknown option: %s\n", argv[i]);
            return 1;
        } else if (!input) {
            input = argv[i];
        } else {
            fprintf(stderr, "error: unexpected argument: %s\n", argv[i]);
            return 1;
        }
    }

    FILE *inf = stdin;
    if (input) {
        inf = fopen(input, "rb");
        if (!inf) {
            fprintf(stderr, "error: cannot open %s: %s\n", input, strerror(errno));
            return 1;
        }
    }

    /* peek at the header to decide between streaming and buffered paths */
    uint8_t hdr_peek[BGE_AAD_LEN];
    size_t hdr_n = fread(hdr_peek, 1, BGE_AAD_LEN, inf);

    int is_binary_v3 = (hdr_n == BGE_AAD_LEN &&
                        memcmp(hdr_peek, BGE_MAGIC, 3) == 0 &&
                        hdr_peek[3] == BGE_STREAM_VERSION);

    if (is_binary_v3) {
        /*  non-armored v3: true streaming decrypt from FILE*  */
        bge_header hdr;
        memcpy(&hdr, hdr_peek, sizeof(hdr));

        uint32_t m_cost = cfx_load32_le(&hdr.m_cost);
        uint32_t t_cost = cfx_load32_le(&hdr.t_cost);
        uint32_t p      = cfx_load32_le(&hdr.p_cost);

        if (m_cost < BGE_MIN_M || t_cost < 1 || p < 1 ||
            m_cost > BGE_MAX_M || t_cost > BGE_MAX_T || p > BGE_MAX_P) {
            fprintf(stderr, "error: unreasonable argon2 parameters\n");
            if (input) fclose(inf);
            return 1;
        }

        char pwd[256] = {0};
        int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
        if (pwd_len <= 0) {
            fprintf(stderr, "error: passphrase required\n");
            cfx_memzero_s(pwd, sizeof(pwd));
            if (input) fclose(inf);
            return 1;
        }

        uint8_t kdf_out[48];

        int rc = cfx_argon2id(kdf_out, 48,
            (const uint8_t *)pwd, (size_t)pwd_len,
            hdr.salt, sizeof(hdr.salt), m_cost, t_cost, p);
        cfx_memzero_s(pwd, sizeof(pwd));
        if (rc != 0) {
            fprintf(stderr, "error: argon2 key derivation failed\n");
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            if (input) fclose(inf);
            return 1;
        }

        const uint8_t *stored_verifier = hdr_peek + BGE_HEADER_LEN;
        uint8_t diff = 0;
        for (int i = 0; i < BGE_VERIFIER_LEN; i++)
            diff |= kdf_out[32 + i] ^ stored_verifier[i];

        if (diff != 0) {
            fprintf(stderr, "error: wrong passphrase\n");
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            if (input) fclose(inf);
            return 1;
        }

        FILE *outf = stdout;
        if (output) {
            outf = fopen(output, "wb");
            if (!outf) {
                fprintf(stderr, "error: cannot open %s: %s\n", output, strerror(errno));
                cfx_memzero_s(kdf_out, sizeof(kdf_out));
                if (input) fclose(inf);
                return 1;
            }
        }

        int ret = bge_decrypt_v3_stream(inf, kdf_out, hdr.nonce, outf);
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        cfx_memzero_s(hdr_peek, sizeof(hdr_peek));
        if (input) fclose(inf);
        if (output) fclose(outf);
        return ret != 0 ? 1 : 0;
    }

    /*  armored or v2: read everything into memory  */
    uint8_t *rest = NULL;
    size_t rest_len = 0;
    int rc = bge_read_all(inf, &rest, &rest_len);
    if (input) fclose(inf);
    if (rc != 0) {
        fprintf(stderr, "error: cannot read input\n");
        return 1;
    }

    /* prepend the bytes we already peeked */
    size_t raw_len = hdr_n + rest_len;
    uint8_t *raw = malloc(raw_len);
    if (!raw) {
        fprintf(stderr, "error: allocation failed\n");
        if (rest) { cfx_memzero_s(rest, rest_len); free(rest); }
        return 1;
    }
    memcpy(raw, hdr_peek, hdr_n);
    if (rest_len > 0) memcpy(raw + hdr_n, rest, rest_len);
    if (rest) { cfx_memzero_s(rest, rest_len); free(rest); }
    cfx_memzero_s(hdr_peek, sizeof(hdr_peek));

    /* auto-detect and strip armor */
    uint8_t *file_buf = raw;
    size_t file_len = raw_len;
    int decoded_armor = 0;

    if (bge_is_armored(raw, raw_len)) {
        uint8_t *dec = NULL;
        size_t dec_len = 0;
        if (bge_armor_decode(raw, raw_len, &dec, &dec_len) != 0) {
            fprintf(stderr, "error: invalid armored input\n");
            cfx_memzero_s(raw, raw_len); free(raw);
            return 1;
        }
        file_buf = dec;
        file_len = dec_len;
        decoded_armor = 1;
    }

    /* validate header */
    if (file_len < BGE_AAD_LEN + BGE_TAG_LEN ||
        memcmp(file_buf, BGE_MAGIC, 3) != 0) {
        fprintf(stderr, "error: not a valid BGE file\n");
        goto fail;
    }

    uint8_t version = file_buf[3];
    if (version != BGE_VERSION && version != BGE_STREAM_VERSION) {
        fprintf(stderr, "error: unsupported BGE version %u\n", version);
        goto fail;
    }

    char pwd[256] = {0};
    int pwd_len = bge_read_passphrase("Enter passphrase: ", pwd, sizeof(pwd));
    if (pwd_len <= 0) {
        fprintf(stderr, "error: passphrase required\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        goto fail;
    }

    /* authenticate via argon2 + verifier */
    bge_header hdr;
    memcpy(&hdr, file_buf, sizeof(hdr));

    uint32_t m_cost = cfx_load32_le(&hdr.m_cost);
    uint32_t t_cost = cfx_load32_le(&hdr.t_cost);
    uint32_t p      = cfx_load32_le(&hdr.p_cost);

    if (m_cost < BGE_MIN_M || t_cost < 1 || p < 1 ||
        m_cost > BGE_MAX_M || t_cost > BGE_MAX_T || p > BGE_MAX_P) {
        fprintf(stderr, "error: unreasonable argon2 parameters\n");
        cfx_memzero_s(pwd, sizeof(pwd));
        goto fail;
    }

    uint8_t kdf_out[48];
    rc = cfx_argon2id(kdf_out, 48,
        (const uint8_t *)pwd, (size_t)pwd_len,
        hdr.salt, sizeof(hdr.salt), m_cost, t_cost, p);
    cfx_memzero_s(pwd, sizeof(pwd));
    if (rc != 0) {
        fprintf(stderr, "error: argon2 key derivation failed\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        goto fail;
    }

    const uint8_t *stored_verifier = file_buf + BGE_HEADER_LEN;
    uint8_t diff = 0;
    for (int i = 0; i < BGE_VERIFIER_LEN; i++)
        diff |= kdf_out[32 + i] ^ stored_verifier[i];

    if (diff != 0) {
        fprintf(stderr, "error: wrong passphrase\n");
        cfx_memzero_s(kdf_out, sizeof(kdf_out));
        goto fail;
    }

    /* open output */
    FILE *outf = stdout;
    if (output) {
        outf = fopen(output, "wb");
        if (!outf) {
            fprintf(stderr, "error: cannot open %s: %s\n", output, strerror(errno));
            cfx_memzero_s(kdf_out, sizeof(kdf_out));
            goto fail;
        }
    }

    /* dispatch by version */
    int ret;
    if (version == BGE_VERSION) {
        ret = bge_decrypt_v2(file_buf, file_len, kdf_out, hdr.nonce, outf);
    } else {
        ret = bge_decrypt_v3(file_buf, file_len, kdf_out, hdr.nonce, outf);
    }

    cfx_memzero_s(kdf_out, sizeof(kdf_out));
    if (output) fclose(outf);
    cfx_memzero_s(file_buf, file_len);
    if (decoded_armor) free(file_buf);
    cfx_memzero_s(raw, raw_len); free(raw);
    return ret != 0 ? 1 : 0;

fail:
    cfx_memzero_s(file_buf, file_len);
    if (decoded_armor) free(file_buf);
    cfx_memzero_s(raw, raw_len); free(raw);
    return 1;
}
