/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_cli_store.c - buffer-level regression tests for BGE store encrypt/decrypt */

#include "cfx_store_internal.h"

#define PWD      "testpass"
#define PWD_LEN  8

/* minimal argon2 params for fast tests */
#define T_M 64
#define T_T 1
#define T_P 1

/* encrypt empty store, decrypt, verify zero-length plaintext */
static void test_empty_store_roundtrip(void) {
    uint8_t *blob = NULL;
    size_t blob_len = 0;

    int rc = bge_v4_init_to_buf(&blob, &blob_len, NULL, 0, PWD, PWD_LEN, T_M, T_T, T_P);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(blob != NULL);

    /* verify magic and version */
    CFX_ASSERT(memcmp(blob, BGE_MAGIC, 3) == 0);
    CFX_ASSERT(blob[3] == BGE_STORE_VERSION);

    /* authenticate and decrypt */
    bge_v4_store store;
    rc = bge_v4_authenticate_buf(blob, blob_len, PWD, PWD_LEN, &store);
    CFX_ASSERT(rc == 0);

    uint8_t *pt = NULL;
    size_t pt_len = 0;
    rc = bge_v4_decrypt_store(&store, &pt, &pt_len);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(pt_len == 0);

    free(pt);
    free(blob);
    bge_v4_store_wipe(&store);
    printf("[OK] test_empty_store_roundtrip\n");
}

/* encrypt a KV store, decrypt, verify plaintext matches */
static void test_kv_roundtrip(void) {
    /* build plaintext: "token" -> "abc123" */
    uint8_t pt_buf[256];
    uint8_t *w = pt_buf;
    cfx_store16_le(w, 5);     w += 2;
    memcpy(w, "token", 5);    w += 5;
    cfx_store32_le(w, 6);     w += 4;
    memcpy(w, "abc123", 6);   w += 6;
    size_t pt_len = (size_t)(w - pt_buf);

    uint8_t *blob = NULL;
    size_t blob_len = 0;

    int rc = bge_v4_init_to_buf(&blob, &blob_len, pt_buf, pt_len, PWD, PWD_LEN, T_M, T_T, T_P);
    CFX_ASSERT(rc == 0);

    bge_v4_store store;
    rc = bge_v4_authenticate_buf(blob, blob_len, PWD, PWD_LEN, &store);
    CFX_ASSERT(rc == 0);

    uint8_t *pt_out = NULL;
    size_t pt_out_len = 0;
    rc = bge_v4_decrypt_store(&store, &pt_out, &pt_out_len);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(pt_out_len == pt_len);
    CFX_ASSERT(memcmp(pt_out, pt_buf, pt_len) == 0);

    /* verify we can read back the value via store_get */
    size_t vlen = 0;
    const uint8_t *val = store_get(pt_out, pt_out_len, "token", &vlen);
    CFX_ASSERT(val != NULL);
    CFX_ASSERT(vlen == 6);
    CFX_ASSERT(memcmp(val, "abc123", 6) == 0);

    free(pt_out);
    free(blob);
    bge_v4_store_wipe(&store);
    printf("[OK] test_kv_roundtrip\n");
}

/* wrong password must fail authentication */
static void test_wrong_password(void) {
    uint8_t pt[] = "secret";
    size_t pt_len = 6;

    uint8_t *blob = NULL;
    size_t blob_len = 0;

    int rc = bge_v4_init_to_buf(&blob, &blob_len, pt, pt_len, PWD, PWD_LEN, T_M, T_T, T_P);
    CFX_ASSERT(rc == 0);

    bge_v4_store store;
    rc = bge_v4_authenticate_buf(blob, blob_len, "wrongpwd", 8, &store);
    CFX_ASSERT(rc != 0);

    free(blob);
    printf("[OK] test_wrong_password\n");
}

/* tampered ciphertext must fail decryption */
static void test_tampered_ciphertext(void) {
    uint8_t pt[] = "secret";
    size_t pt_len = 6;

    uint8_t *blob = NULL;
    size_t blob_len = 0;

    int rc = bge_v4_init_to_buf(&blob, &blob_len, pt, pt_len, PWD, PWD_LEN, T_M, T_T, T_P);
    CFX_ASSERT(rc == 0);

    /* flip a byte in the ciphertext region (after header + slot + nonce) */
    size_t ct_offset = BGE_STORE_HDR_LEN + BGE_STORE_SLOT_LEN + 24;
    blob[ct_offset] ^= 0xff;

    bge_v4_store store;
    rc = bge_v4_authenticate_buf(blob, blob_len, PWD, PWD_LEN, &store);
    CFX_ASSERT(rc == 0);  /* auth passes (verifier is fine) */

    uint8_t *pt_out = NULL;
    size_t pt_out_len = 0;
    rc = bge_v4_decrypt_store(&store, &pt_out, &pt_out_len);
    CFX_ASSERT(rc != 0);  /* decryption must fail (AEAD tag mismatch) */

    free(blob);
    bge_v4_store_wipe(&store);
    printf("[OK] test_tampered_ciphertext\n");
}

/* multiple entries: set, get, rm via store_* helpers */
static void test_store_set_get_rm(void) {
    /* start with empty plaintext */
    size_t pt_len = 0;
    uint8_t *pt = NULL;

    /* set "user" -> "alice" */
    size_t new_len;
    uint8_t *pt2 = store_set(pt, pt_len, "user",
        (const uint8_t *)"alice", 5, &new_len);
    CFX_ASSERT(pt2 != NULL);
    free(pt);
    pt = pt2; pt_len = new_len;

    /* set "pass" -> "hunter2" */
    pt2 = store_set(pt, pt_len, "pass",
        (const uint8_t *)"hunter2", 7, &new_len);
    CFX_ASSERT(pt2 != NULL);
    free(pt);
    pt = pt2; pt_len = new_len;

    /* get "user" */
    size_t vlen;
    const uint8_t *v = store_get(pt, pt_len, "user", &vlen);
    CFX_ASSERT(v && vlen == 5 && memcmp(v, "alice", 5) == 0);

    /* get "pass" */
    v = store_get(pt, pt_len, "pass", &vlen);
    CFX_ASSERT(v && vlen == 7 && memcmp(v, "hunter2", 7) == 0);

    /* encrypt roundtrip the whole store */
    uint8_t *blob = NULL;
    size_t blob_len = 0;
    int rc = bge_v4_init_to_buf(&blob, &blob_len, pt, pt_len, PWD, PWD_LEN, T_M, T_T, T_P);
    CFX_ASSERT(rc == 0);

    bge_v4_store store;
    rc = bge_v4_authenticate_buf(blob, blob_len, PWD, PWD_LEN, &store);
    CFX_ASSERT(rc == 0);

    uint8_t *dec = NULL;
    size_t dec_len = 0;
    rc = bge_v4_decrypt_store(&store, &dec, &dec_len);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(dec_len == pt_len);
    CFX_ASSERT(memcmp(dec, pt, pt_len) == 0);

    /* rm "user", verify only "pass" remains */
    pt2 = store_rm(dec, dec_len, "user", &new_len);
    CFX_ASSERT(pt2 != NULL);
    CFX_ASSERT(store_get(pt2, new_len, "user", &vlen) == NULL);
    v = store_get(pt2, new_len, "pass", &vlen);
    CFX_ASSERT(v && vlen == 7);

    free(pt2);
    free(dec);
    free(blob);
    free(pt);
    bge_v4_store_wipe(&store);
    printf("[OK] test_store_set_get_rm\n");
}

/* overwrite existing key */
static void test_store_overwrite(void) {
    size_t len = 0;
    uint8_t *pt = NULL;
    size_t new_len;

    pt = store_set(pt, len, "key", (const uint8_t *)"old", 3, &new_len);
    CFX_ASSERT(pt);
    len = new_len;

    uint8_t *pt2 = store_set(pt, len, "key", (const uint8_t *)"new_value", 9, &new_len);
    CFX_ASSERT(pt2);
    free(pt);
    pt = pt2; len = new_len;

    size_t vlen;
    const uint8_t *v = store_get(pt, len, "key", &vlen);
    CFX_ASSERT(v && vlen == 9 && memcmp(v, "new_value", 9) == 0);

    free(pt);
    printf("[OK] test_store_overwrite\n");
}

int main(void) {
    test_empty_store_roundtrip();
    test_kv_roundtrip();
    test_wrong_password();
    test_tampered_ciphertext();
    test_store_set_get_rm();
    test_store_overwrite();
    printf("All store tests passed.\n");
    return 0;
}
