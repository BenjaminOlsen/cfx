/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/aead_chacha20_poly1271.h"
#include "cfx/macros.h"
#include "cfx/rand.h"

#include <string.h>
#include <stdint.h>
#include <stdio.h>

static void print_hex(const char *label, const uint8_t *buf, size_t len) {
    printf("%s: ", label);
    for (size_t i = 0; i < len && i < 32; ++i) printf("%02x", buf[i]);
    if (len > 32) printf("...");
    printf("\n");
}

static void test_basic_roundtrip(void) {
    printf("== test_basic_roundtrip ==\n");

    const uint8_t key[32] = {
        0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
        0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
        0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
        0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f
    };
    const uint8_t nonce[12] = {0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x4a,0x00,0x00,0x00,0x00};
    const uint8_t pt[] = "Ladies and Gentlemen of the class of '99: "
                         "If I could offer you only one tip for the future, "
                         "sunscreen would be it.";
    const uint8_t aad[] = "Additional authenticated data";
    size_t pt_len = strlen((const char *)pt);
    size_t aad_len = strlen((const char *)aad);

    uint8_t ct[256], tag[16], pt_out[256];

    int rc = cfx_chacha20_poly1271_encrypt(ct, tag, pt, pt_len, aad, aad_len, key, nonce);
    CFX_ASSERT(rc == 0);
    print_hex("ct", ct, pt_len);
    print_hex("tag", tag, 16);

    rc = cfx_chacha20_poly1271_decrypt(pt_out, ct, pt_len, aad, aad_len, key, nonce, tag);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(pt_out, pt, pt_len) == 0);
    printf("   roundtrip OK\n");
}

static void test_bad_tag(void) {
    printf("== test_bad_tag ==\n");

    const uint8_t key[32] = {0};
    const uint8_t nonce[12] = {0};
    const uint8_t pt[] = "Test message for authentication";
    size_t pt_len = strlen((const char *)pt);

    uint8_t ct[64], tag[16], bad_tag[16], pt_out[64];

    int rc = cfx_chacha20_poly1271_encrypt(ct, tag, pt, pt_len, NULL, 0, key, nonce);
    CFX_ASSERT(rc == 0);

    memcpy(bad_tag, tag, 16);
    bad_tag[0] ^= 0x01;

    rc = cfx_chacha20_poly1271_decrypt(pt_out, ct, pt_len, NULL, 0, key, nonce, bad_tag);
    CFX_ASSERT(rc != 0);
    printf("   bad tag rejected OK\n");
}

static void test_empty_plaintext(void) {
    printf("== test_empty_plaintext ==\n");

    const uint8_t key[32] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
                             17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32};
    const uint8_t nonce[12] = {0};
    const uint8_t aad[] = "Only AAD, no plaintext";
    size_t aad_len = strlen((const char *)aad);
    uint8_t tag[16], dummy;

    int rc = cfx_chacha20_poly1271_encrypt(&dummy, tag, NULL, 0, aad, aad_len, key, nonce);
    CFX_ASSERT(rc == 0);
    rc = cfx_chacha20_poly1271_decrypt(&dummy, NULL, 0, aad, aad_len, key, nonce, tag);
    CFX_ASSERT(rc == 0);
    printf("   empty plaintext OK\n");
}

static void test_empty_aad(void) {
    printf("== test_empty_aad ==\n");

    const uint8_t key[32] = {0};
    const uint8_t nonce[12] = {0};
    const uint8_t pt[] = "Message without AAD";
    size_t pt_len = strlen((const char *)pt);
    uint8_t ct[64], tag[16], pt_out[64];

    int rc = cfx_chacha20_poly1271_encrypt(ct, tag, pt, pt_len, NULL, 0, key, nonce);
    CFX_ASSERT(rc == 0);
    rc = cfx_chacha20_poly1271_decrypt(pt_out, ct, pt_len, NULL, 0, key, nonce, tag);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(pt_out, pt, pt_len) == 0);
    printf("   empty AAD OK\n");
}

static void test_length_coverage(void) {
    printf("== test_length_coverage ==\n");

    const uint8_t key[32] = {0x42};
    const uint8_t nonce[12] = {0x43};
    size_t test_lengths[] = {0, 1, 14, 15, 16, 29, 30, 31, 44, 45, 46, 100, 255};
    size_t num_tests = sizeof(test_lengths) / sizeof(test_lengths[0]);
    uint8_t pt[256], ct[256], pt_out[256], tag[16], aad[32];

    for (size_t i = 0; i < num_tests; i++) {
        size_t len = test_lengths[i];
        for (size_t j = 0; j < len; j++) pt[j] = (uint8_t)(j ^ i);
        for (size_t j = 0; j < sizeof(aad); j++) aad[j] = (uint8_t)(j + i);

        int rc = cfx_chacha20_poly1271_encrypt(ct, tag, pt, len, aad, sizeof(aad), key, nonce);
        CFX_ASSERT(rc == 0);
        rc = cfx_chacha20_poly1271_decrypt(pt_out, ct, len, aad, sizeof(aad), key, nonce, tag);
        CFX_ASSERT(rc == 0);
        if (len > 0) CFX_ASSERT(memcmp(pt_out, pt, len) == 0);
        printf("   len=%zu OK\n", len);
    }
}

static void test_xchacha_roundtrip(void) {
    printf("== test_xchacha_roundtrip ==\n");

    const uint8_t key[32] = {
        0x80,0x81,0x82,0x83,0x84,0x85,0x86,0x87,
        0x88,0x89,0x8a,0x8b,0x8c,0x8d,0x8e,0x8f,
        0x90,0x91,0x92,0x93,0x94,0x95,0x96,0x97,
        0x98,0x99,0x9a,0x9b,0x9c,0x9d,0x9e,0x9f
    };
    const uint8_t nonce[24] = {
        0x40,0x41,0x42,0x43,0x44,0x45,0x46,0x47,
        0x48,0x49,0x4a,0x4b,0x4c,0x4d,0x4e,0x4f,
        0x50,0x51,0x52,0x53,0x54,0x55,0x56,0x57
    };
    const uint8_t pt[] = "XChaCha20-Poly1271 test message with 24-byte nonce";
    const uint8_t aad[] = "XChaCha AAD";
    size_t pt_len = strlen((const char *)pt);
    size_t aad_len = strlen((const char *)aad);
    uint8_t ct[128], tag[16], pt_out[128], bad_tag[16];

    int rc = cfx_xchacha20_poly1271_encrypt(ct, tag, pt, pt_len, aad, aad_len, key, nonce);
    CFX_ASSERT(rc == 0);
    rc = cfx_xchacha20_poly1271_decrypt(pt_out, ct, pt_len, aad, aad_len, key, nonce, tag);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(memcmp(pt_out, pt, pt_len) == 0);

    memcpy(bad_tag, tag, 16);
    bad_tag[8] ^= 0x80;
    rc = cfx_xchacha20_poly1271_decrypt(pt_out, ct, pt_len, aad, aad_len, key, nonce, bad_tag);
    CFX_ASSERT(rc != 0);
    printf("   xchacha roundtrip OK\n");
}

static void test_fuzz(void) {
    printf("== test_fuzz ==\n");

    #define MAX_PT 512
    #define MAX_AAD 128
    #define FUZZ_ITERS 5000

    uint8_t key[32], nonce[12], aad[MAX_AAD], pt[MAX_PT], ct[MAX_PT], pt_out[MAX_PT], tag[16], bad_tag[16];
    cfx_srand(0xDEADBEEF);

    for (int iter = 0; iter < FUZZ_ITERS; iter++) {
        for (size_t i = 0; i < sizeof key; i++) key[i] = (uint8_t)cfx_rand();
        for (size_t i = 0; i < sizeof nonce; i++) nonce[i] = (uint8_t)cfx_rand();

        size_t pt_len = cfx_rand() % (MAX_PT + 1);
        size_t aad_len = cfx_rand() % (MAX_AAD + 1);
        for (size_t i = 0; i < pt_len; i++) pt[i] = (uint8_t)cfx_rand();
        for (size_t i = 0; i < aad_len; i++) aad[i] = (uint8_t)cfx_rand();

        const uint8_t *aad_ptr = aad_len ? aad : NULL;

        int rc = cfx_chacha20_poly1271_encrypt(ct, tag, pt, pt_len, aad_ptr, aad_len, key, nonce);
        CFX_ASSERT(rc == 0);
        rc = cfx_chacha20_poly1271_decrypt(pt_out, ct, pt_len, aad_ptr, aad_len, key, nonce, tag);
        CFX_ASSERT(rc == 0);
        if (pt_len > 0) CFX_ASSERT(memcmp(pt_out, pt, pt_len) == 0);

        memcpy(bad_tag, tag, 16);
        bad_tag[cfx_rand() % 16] ^= (uint8_t)(1u << (cfx_rand() & 7u));
        rc = cfx_chacha20_poly1271_decrypt(pt_out, ct, pt_len, aad_ptr, aad_len, key, nonce, bad_tag);
        CFX_ASSERT(rc != 0);

        if (iter % 500 == 0) printf("   fuzz iter %d OK (pt=%zu, aad=%zu)\n", iter, pt_len, aad_len);
    }
    printf("   fuzz %d iterations OK\n", FUZZ_ITERS);
}

static void test_xchacha_fuzz(void) {
    printf("== test_xchacha_fuzz ==\n");

    #define XMAX_PT 256
    #define XMAX_AAD 64
    #define XFUZZ_ITERS 3000

    uint8_t key[32], nonce[24], aad[XMAX_AAD], pt[XMAX_PT], ct[XMAX_PT], pt_out[XMAX_PT], tag[16], bad_tag[16];
    cfx_srand(0xCAFEBABE);

    for (int iter = 0; iter < XFUZZ_ITERS; iter++) {
        for (size_t i = 0; i < sizeof key; i++) key[i] = (uint8_t)cfx_rand();
        for (size_t i = 0; i < sizeof nonce; i++) nonce[i] = (uint8_t)cfx_rand();

        size_t pt_len = cfx_rand() % (XMAX_PT + 1);
        size_t aad_len = cfx_rand() % (XMAX_AAD + 1);
        for (size_t i = 0; i < pt_len; i++) pt[i] = (uint8_t)cfx_rand();
        for (size_t i = 0; i < aad_len; i++) aad[i] = (uint8_t)cfx_rand();

        const uint8_t *aad_ptr = aad_len ? aad : NULL;

        int rc = cfx_xchacha20_poly1271_encrypt(ct, tag, pt, pt_len, aad_ptr, aad_len, key, nonce);
        CFX_ASSERT(rc == 0);
        rc = cfx_xchacha20_poly1271_decrypt(pt_out, ct, pt_len, aad_ptr, aad_len, key, nonce, tag);
        CFX_ASSERT(rc == 0);
        if (pt_len > 0) CFX_ASSERT(memcmp(pt_out, pt, pt_len) == 0);

        memcpy(bad_tag, tag, 16);
        bad_tag[cfx_rand() % 16] ^= (uint8_t)(1u << (cfx_rand() & 7u));
        rc = cfx_xchacha20_poly1271_decrypt(pt_out, ct, pt_len, aad_ptr, aad_len, key, nonce, bad_tag);
        CFX_ASSERT(rc != 0);

        if (iter % 500 == 0) printf("   xchacha fuzz iter %d OK\n", iter);
    }
    printf("   xchacha fuzz %d iterations OK\n", XFUZZ_ITERS);
}

int main(void) {
    CFX_TEST(test_basic_roundtrip);
    CFX_TEST(test_bad_tag);
    CFX_TEST(test_empty_plaintext);
    CFX_TEST(test_empty_aad);
    CFX_TEST(test_length_coverage);
    CFX_TEST(test_xchacha_roundtrip);
    CFX_TEST(test_fuzz);
    CFX_TEST(test_xchacha_fuzz);
    puts("ok");
    return 0;
}
