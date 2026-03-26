/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/poly1271.h"
#include "cfx/macros.h"
#include "cfx/rand.h"
#include "cfx/arch.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

static void print_hex(const char *label, const uint8_t *buf, size_t len) {
    printf("%s: ", label);
    for (size_t i = 0; i < len; i++) printf("%02x", buf[i]);
    printf("\n");
}

static void test_empty_message(void) {
    uint8_t key[32] = {
        0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08,
        0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10,
        0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18,
        0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f, 0x20
    };
    uint8_t tag1[16], tag2[16];

    cfx_poly1271(tag1, NULL, 0, key);
    cfx_poly1271(tag2, NULL, 0, key);

    print_hex("empty tag", tag1, 16);
    CFX_ASSERT(memcmp(tag1, tag2, 16) == 0);
}

static void test_single_block(void) {
    uint8_t key[32] = {0};
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)i;

    uint8_t msg[15] = "Hello, Poly127";
    uint8_t tag1[16], tag2[16];

    cfx_poly1271(tag1, msg, 15, key);
    cfx_poly1271(tag2, msg, 15, key);

    print_hex("single block tag", tag1, 16);
    CFX_ASSERT(memcmp(tag1, tag2, 16) == 0);
}

static void test_multiple_blocks(void) {
    uint8_t key[32] = {0};
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)(i * 3 + 7);

    uint8_t msg[100];
    for (int i = 0; i < 100; i++) msg[i] = (uint8_t)i;

    uint8_t tag1[16], tag2[16];

    cfx_poly1271(tag1, msg, 100, key);
    cfx_poly1271(tag2, msg, 100, key);

    print_hex("multi block tag", tag1, 16);
    CFX_ASSERT(memcmp(tag1, tag2, 16) == 0);
}

static void test_streaming(void) {
    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)(i ^ 0xAA);

    uint8_t msg[100];
    for (int i = 0; i < 100; i++) msg[i] = (uint8_t)(i * 7);

    uint8_t tag_oneshot[16];
    cfx_poly1271(tag_oneshot, msg, 100, key);

    size_t chunk_sizes[] = {1, 3, 7, 13, 15, 16, 17, 31, 50};
    for (size_t c = 0; c < sizeof(chunk_sizes)/sizeof(chunk_sizes[0]); c++) {
        size_t chunk = chunk_sizes[c];
        cfx_poly1271_ctx_t ctx;
        cfx_poly1271_init(&ctx, key);

        size_t pos = 0;
        while (pos < 100) {
            size_t n = (pos + chunk > 100) ? (100 - pos) : chunk;
            cfx_poly1271_update(&ctx, msg + pos, n);
            pos += n;
        }

        uint8_t tag_stream[16];
        cfx_poly1271_finish(&ctx, tag_stream);

        CFX_ASSERT(memcmp(tag_oneshot, tag_stream, 16) == 0);
    }
    printf("streaming test passed for all chunk sizes\n");
}

static void test_different_messages(void) {
    uint8_t key[32] = {0};
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)i;

    uint8_t msg1[20] = "message one here!!!";
    uint8_t msg2[20] = "message two here!!!";

    uint8_t tag1[16], tag2[16];
    cfx_poly1271(tag1, msg1, 20, key);
    cfx_poly1271(tag2, msg2, 20, key);

    print_hex("msg1 tag", tag1, 16);
    print_hex("msg2 tag", tag2, 16);
    CFX_ASSERT(memcmp(tag1, tag2, 16) != 0);
}

static void test_different_keys(void) {
    uint8_t key1[32] = {0};
    uint8_t key2[32] = {0};
    for (int i = 0; i < 32; i++) {
        key1[i] = (uint8_t)i;
        key2[i] = (uint8_t)(i + 1);
    }

    uint8_t msg[30] = "test message for key diff!!!";
    uint8_t tag1[16], tag2[16];

    cfx_poly1271(tag1, msg, 30, key1);
    cfx_poly1271(tag2, msg, 30, key2);

    print_hex("key1 tag", tag1, 16);
    print_hex("key2 tag", tag2, 16);
    CFX_ASSERT(memcmp(tag1, tag2, 16) != 0);
}

static void test_all_lengths(void) {
    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)(i * 11);

    uint8_t msg[200];
    for (int i = 0; i < 200; i++) msg[i] = (uint8_t)(i * 13);

    uint8_t prev_tag[16] = {0};

    for (size_t len = 0; len <= 60; len++) {
        uint8_t tag[16];
        cfx_poly1271(tag, msg, len, key);

        if (len > 0) {
            CFX_ASSERT(memcmp(tag, prev_tag, 16) != 0);
        }
        memcpy(prev_tag, tag, 16);
    }
    printf("all lengths 0-60 produce unique tags\n");
}

static void test_fuzz(void) {
    #define FUZZ_ITERS 5000
    #define FUZZ_MAX_LEN 500

    cfx_srand(0xDEADBEEF);

    for (int iter = 0; iter < FUZZ_ITERS; iter++) {
        uint8_t key[32];
        for (int i = 0; i < 32; i++) key[i] = (uint8_t)cfx_rand();

        size_t len = cfx_rand() % (FUZZ_MAX_LEN + 1);
        uint8_t *msg = (uint8_t *)malloc(len + 1);
        CFX_ASSERT(msg != NULL);
        for (size_t i = 0; i < len; i++) msg[i] = (uint8_t)cfx_rand();

        uint8_t tag1[16], tag2[16];
        cfx_poly1271(tag1, msg, len, key);

        cfx_poly1271_ctx_t ctx;
        cfx_poly1271_init(&ctx, key);
        size_t pos = 0;
        while (pos < len) {
            size_t chunk = (cfx_rand() % 50) + 1;
            if (pos + chunk > len) chunk = len - pos;
            cfx_poly1271_update(&ctx, msg + pos, chunk);
            pos += chunk;
        }
        cfx_poly1271_finish(&ctx, tag2);

        CFX_ASSERT(memcmp(tag1, tag2, 16) == 0);

        free(msg);

        if (iter % (FUZZ_ITERS / 10) == 0) {
            printf("fuzz iter %d ok\n", iter);
        }
    }
    printf("fuzz test passed (%d iterations)\n", FUZZ_ITERS);
}

static void test_zeros(void) {
    uint8_t key[32] = {0};
    key[0] = 1;

    uint8_t msg[100] = {0};
    uint8_t tag[16];

    cfx_poly1271(tag, msg, 100, key);
    print_hex("zeros tag", tag, 16);

    int all_zero = 1;
    for (int i = 0; i < 16; i++) if (tag[i] != 0) all_zero = 0;
    CFX_ASSERT(!all_zero);
}

static void test_ones(void) {
    uint8_t key[32];
    memset(key, 0xFF, 32);

    uint8_t msg[100];
    memset(msg, 0xFF, 100);
    uint8_t tag[16];

    cfx_poly1271(tag, msg, 100, key);
    print_hex("ones tag", tag, 16);
}

/* test vectors generated by Python reference implementation */
static void hex_to_bytes(const char *hex, uint8_t *out, size_t len) {
    for (size_t i = 0; i < len; i++) {
        unsigned int byte;
        sscanf(hex + 2*i, "%02x", &byte);
        out[i] = (uint8_t)byte;
    }
}

static int check_vector(const char *name, const char *key_hex,
    const uint8_t *msg, size_t msg_len,
    const char *expected_hex) {
    uint8_t key[32], expected[16], tag[16];

    hex_to_bytes(key_hex, key, 32);
    hex_to_bytes(expected_hex, expected, 16);

    cfx_poly1271(tag, msg, msg_len, key);

    if (memcmp(tag, expected, 16) != 0) {
        printf("FAIL %s\n", name);
        print_hex("  expected", expected, 16);
        print_hex("  got     ", tag, 16);
        return 0;
    }
    printf("  %s: ok\n", name);
    return 1;
}

static void test_verify(void) {
    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)(i * 3);

    uint8_t msg[] = "Test message for verification";
    uint8_t tag[16];

    cfx_poly1271(tag, msg, sizeof(msg) - 1, key);

    /* valid tag should verify */
    CFX_ASSERT(cfx_poly1271_verify(tag, msg, sizeof(msg) - 1, key) == 0);

    /* modified tag should fail */
    tag[5] ^= 0x01;
    CFX_ASSERT(cfx_poly1271_verify(tag, msg, sizeof(msg) - 1, key) != 0);

    /* restore tag, modify message */
    tag[5] ^= 0x01;
    msg[0] ^= 0x01;
    CFX_ASSERT(cfx_poly1271_verify(tag, msg, sizeof(msg) - 1, key) != 0);

    printf("verify test passed\n");
}

static void test_reference_vectors(void) {
    /* vectors from Python reference implementation (poly1271_verify.py) */
    int passed = 1;

    /* length 0 */
    passed &= check_vector("len_0",
        "b5739948a249856c49e54909ebb2f31d497377aea932d57ae80e81139e4bb6dd",
        NULL, 0,
        "497377aea932d57ae80e81139e4bb6dd");

    /* length 1: msg = 0x0d */
    {
        uint8_t msg[] = {0x0d};
        passed &= check_vector("len_1",
            "f298f36369b075bf9168ceda5e9500704ac3c1475720556168946cb49dfcf6d3",
            msg, 1,
            "9479b96ea37dff9fc873500f55ee93d4");
    }

    /* length 15 (exactly one block) */
    {
        uint8_t msg[15];
        for (int i = 0; i < 15; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        passed &= check_vector("len_15",
            "f6118ce0bcdae1277b9b82025e5aeb46c9c4015d0f8b0ea4cdb292950dbca174",
            msg, 15,
            "00bd5684bbb9ec0fdb6812fec2b683da");
    }

    /* length 16 (one block + 1 byte) */
    {
        uint8_t msg[16];
        for (int i = 0; i < 16; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        passed &= check_vector("len_16",
            "bdce2ecf9a9f4388b2e3c055eb3940f7f6ecfcd6405806aa01b11984937b4a5f",
            msg, 16,
            "1a7e08badf6a311a1326e180535fe5b0");
    }

    /* length 30 (exactly two blocks) */
    {
        uint8_t msg[30];
        for (int i = 0; i < 30; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        passed &= check_vector("len_30",
            "cef0ca2975a63992cd201d7a474f4f37d284bc4c4c402ccecafe52f8eac671d3",
            msg, 30,
            "d64d14c694fc3ade68fa7e6ee9862108");
    }

    /* length 31 (two blocks + 1 byte) */
    {
        uint8_t msg[31];
        for (int i = 0; i < 31; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        passed &= check_vector("len_31",
            "ecd5a699516a991bd618eedacbb7f7ba9815587976e7a8a7ddfd7790c077b09e",
            msg, 31,
            "4bf292fab6af8bff3202c45b9e9a8701");
    }

    /* length 100 */
    {
        uint8_t msg[100];
        for (int i = 0; i < 100; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        passed &= check_vector("len_100",
            "b0fdf754bb59965e413060dedf8ed9b0d96e8c23544792dca17e89998da75a3a",
            msg, 100,
            "68977a391e9bdf4d904fdb851cde1599");
    }

    /* length 256 */
    {
        uint8_t msg[256];
        for (int i = 0; i < 256; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        passed &= check_vector("len_256",
            "595f063a7941f6cb8b74934e31a7c0f477b7240be9f726a321ee09b53eadf174",
            msg, 256,
            "1fff8a650754ea5eb55963a614335ab0");
    }

    /* length 1000 */
    {
        uint8_t *msg = (uint8_t *)malloc(1000);
        CFX_ASSERT(msg != NULL);
        for (int i = 0; i < 1000; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        passed &= check_vector("len_1000",
            "5476a825f45840087522474e83e1ec8d70a34954735d043659b18aa9244668cf",
            msg, 1000,
            "53d0bba62df4206e277f411f099fa60c");
        free(msg);
    }

    /*
     * Walt Whitman - "O Me! O Life!" (1892) - last stanza
     * "That the powerful play goes on, and you may contribute a verse."
     */
    {
        const char *whitman =
            "Answer.\n"
            "That you are here-that life exists and identity,\n"
            "That the powerful play goes on, and you may contribute a verse.";
        passed &= check_vector("uncle_walt",
            "c20ae296d65cb1ef04ef2cd32acdeb835f9bfacd77a2ff311b1e573f9f949677",
            (const uint8_t *)whitman, 120,
            "e0840ca53a6744daa78e6d9c65c3bbe9");
    }

    CFX_ASSERT(passed);
}

/* AVX2 tests */

#if CFX_HAVE_AVX2

static void test_avx2_vs_scalar(void) {
    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)(i * 3 + 7);

    /* test various lengths */
    size_t lens[] = {0, 1, 14, 15, 16, 29, 30, 31, 45, 59, 60, 61, 100, 256, 1000};

    for (size_t i = 0; i < sizeof(lens)/sizeof(lens[0]); i++) {
        size_t len = lens[i];
        uint8_t *msg = len > 0 ? (uint8_t *)malloc(len) : NULL;
        if (len > 0) {
            CFX_ASSERT(msg != NULL);
            for (size_t j = 0; j < len; j++) msg[j] = (uint8_t)((j * 7 + 13) & 0xff);
        }

        uint8_t tag_scalar[16], tag_avx2[16];
        cfx_poly1271(tag_scalar, msg, len, key);
        cfx_poly1271_avx2(tag_avx2, msg, len, key);

        if (memcmp(tag_scalar, tag_avx2, 16) != 0) {
            printf("FAIL avx2 vs scalar len_%zu\n", len);
            print_hex("  scalar", tag_scalar, 16);
            print_hex("  avx2  ", tag_avx2, 16);
            CFX_ASSERT(0);
        }

        if (len > 0) free(msg);
    }
    printf("avx2 vs scalar: all lengths match\n");
}

static void test_avx2_streaming(void) {
    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)i;

    uint8_t msg[200];
    for (int i = 0; i < 200; i++) msg[i] = (uint8_t)i;

    uint8_t tag_oneshot[16];
    cfx_poly1271_avx2(tag_oneshot, msg, 200, key);

    /* test various chunk sizes */
    size_t chunks[] = {1, 7, 15, 16, 30, 60, 61};
    for (size_t c = 0; c < sizeof(chunks)/sizeof(chunks[0]); c++) {
        size_t chunk = chunks[c];
        cfx_poly1271_avx2_ctx_t ctx;
        cfx_poly1271_avx2_init(&ctx, key);

        size_t pos = 0;
        while (pos < 200) {
            size_t n = (pos + chunk > 200) ? (200 - pos) : chunk;
            cfx_poly1271_avx2_update(&ctx, msg + pos, n);
            pos += n;
        }

        uint8_t tag[16];
        cfx_poly1271_avx2_finish(&ctx, tag);

        CFX_ASSERT(memcmp(tag_oneshot, tag, 16) == 0);
    }
    printf("avx2 streaming: all chunk sizes match\n");
}

static void test_avx2_reference_vectors(void) {
    int passed = 1;

    /* reuse same vectors as scalar */
    {
        uint8_t key[32], expected[16], tag[16];
        hex_to_bytes("b5739948a249856c49e54909ebb2f31d497377aea932d57ae80e81139e4bb6dd", key, 32);
        hex_to_bytes("497377aea932d57ae80e81139e4bb6dd", expected, 16);
        cfx_poly1271_avx2(tag, NULL, 0, key);
        if (memcmp(tag, expected, 16) != 0) {
            printf("FAIL avx2 len_0\n");
            passed = 0;
        }
    }

    {
        uint8_t key[32], expected[16], tag[16], msg[] = {0x0d};
        hex_to_bytes("f298f36369b075bf9168ceda5e9500704ac3c1475720556168946cb49dfcf6d3", key, 32);
        hex_to_bytes("9479b96ea37dff9fc873500f55ee93d4", expected, 16);
        cfx_poly1271_avx2(tag, msg, 1, key);
        if (memcmp(tag, expected, 16) != 0) {
            printf("FAIL avx2 len_1\n");
            passed = 0;
        }
    }

    {
        uint8_t key[32], expected[16], tag[16], msg[256];
        for (int i = 0; i < 256; i++) msg[i] = (uint8_t)((i * 7 + 13) & 0xff);
        hex_to_bytes("595f063a7941f6cb8b74934e31a7c0f477b7240be9f726a321ee09b53eadf174", key, 32);
        hex_to_bytes("1fff8a650754ea5eb55963a614335ab0", expected, 16);
        cfx_poly1271_avx2(tag, msg, 256, key);
        if (memcmp(tag, expected, 16) != 0) {
            printf("FAIL avx2 len_256\n");
            passed = 0;
        }
    }

    CFX_ASSERT(passed);
    printf("avx2 reference vectors: ok\n");
}

static void test_avx2_verify(void) {
    uint8_t key[32];
    for (int i = 0; i < 32; i++) key[i] = (uint8_t)(i * 3);

    uint8_t msg[] = "Test message for AVX2 verification";
    uint8_t tag[16];

    cfx_poly1271_avx2(tag, msg, sizeof(msg) - 1, key);

    CFX_ASSERT(cfx_poly1271_avx2_verify(tag, msg, sizeof(msg) - 1, key) == 0);

    tag[5] ^= 0x01;
    CFX_ASSERT(cfx_poly1271_avx2_verify(tag, msg, sizeof(msg) - 1, key) != 0);

    printf("avx2 verify: ok\n");
}

#endif /* CFX_HAVE_AVX2 */

int main(void) {
    CFX_TEST(test_empty_message);
    CFX_TEST(test_single_block);
    CFX_TEST(test_multiple_blocks);
    CFX_TEST(test_streaming);
    CFX_TEST(test_different_messages);
    CFX_TEST(test_different_keys);
    CFX_TEST(test_all_lengths);
    CFX_TEST(test_zeros);
    CFX_TEST(test_ones);
    CFX_TEST(test_verify);
    CFX_TEST(test_reference_vectors);
    CFX_TEST(test_fuzz);

#if CFX_HAVE_AVX2
    printf("\nAVX2 Tests\n");
    printf("----------\n");
    CFX_TEST(test_avx2_vs_scalar);
    CFX_TEST(test_avx2_streaming);
    CFX_TEST(test_avx2_reference_vectors);
    CFX_TEST(test_avx2_verify);
#else
    printf("\nAVX2 not available, skipping AVX2 tests\n");
#endif

    puts("ok");
    return 0;
}
