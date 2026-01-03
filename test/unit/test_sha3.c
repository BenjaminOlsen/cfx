/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/sha3.h"
#include "cfx/macros.h"
#include <stdio.h>
#include <string.h>

/* NIST FIPS 202 test vectors */

/* SHA3-256("") */
static void test_sha3_256_empty(void) {
    uint8_t out[32];
    const uint8_t expected[32] = {
        0xa7, 0xff, 0xc6, 0xf8, 0xbf, 0x1e, 0xd7, 0x66,
        0x51, 0xc1, 0x47, 0x56, 0xa0, 0x61, 0xd6, 0x62,
        0xf5, 0x80, 0xff, 0x4d, 0xe4, 0x3b, 0x49, 0xfa,
        0x82, 0xd8, 0x0a, 0x4b, 0x80, 0xf8, 0x43, 0x4a
    };

    cfx_sha3_256(out, NULL, 0);
    CFX_ASSERT(memcmp(out, expected, 32) == 0);
}

/* SHA3-256("abc") */
static void test_sha3_256_abc(void) {
    uint8_t out[32];
    const uint8_t expected[32] = {
        0x3a, 0x98, 0x5d, 0xa7, 0x4f, 0xe2, 0x25, 0xb2,
        0x04, 0x5c, 0x17, 0x2d, 0x6b, 0xd3, 0x90, 0xbd,
        0x85, 0x5f, 0x08, 0x6e, 0x3e, 0x9d, 0x52, 0x5b,
        0x46, 0xbf, 0xe2, 0x45, 0x11, 0x43, 0x15, 0x32
    };

    cfx_sha3_256(out, "abc", 3);
    CFX_ASSERT(memcmp(out, expected, 32) == 0);
}

/* SHA3-224("") */
static void test_sha3_224_empty(void) {
    uint8_t out[28];
    const uint8_t expected[28] = {
        0x6b, 0x4e, 0x03, 0x42, 0x36, 0x67, 0xdb, 0xb7,
        0x3b, 0x6e, 0x15, 0x45, 0x4f, 0x0e, 0xb1, 0xab,
        0xd4, 0x59, 0x7f, 0x9a, 0x1b, 0x07, 0x8e, 0x3f,
        0x5b, 0x5a, 0x6b, 0xc7
    };

    cfx_sha3_224(out, NULL, 0);
    CFX_ASSERT(memcmp(out, expected, 28) == 0);
}

/* SHA3-384("") */
static void test_sha3_384_empty(void) {
    uint8_t out[48];
    const uint8_t expected[48] = {
        0x0c, 0x63, 0xa7, 0x5b, 0x84, 0x5e, 0x4f, 0x7d,
        0x01, 0x10, 0x7d, 0x85, 0x2e, 0x4c, 0x24, 0x85,
        0xc5, 0x1a, 0x50, 0xaa, 0xaa, 0x94, 0xfc, 0x61,
        0x99, 0x5e, 0x71, 0xbb, 0xee, 0x98, 0x3a, 0x2a,
        0xc3, 0x71, 0x38, 0x31, 0x26, 0x4a, 0xdb, 0x47,
        0xfb, 0x6b, 0xd1, 0xe0, 0x58, 0xd5, 0xf0, 0x04
    };

    cfx_sha3_384(out, NULL, 0);
    CFX_ASSERT(memcmp(out, expected, 48) == 0);
}

/* SHA3-512("") */
static void test_sha3_512_empty(void) {
    uint8_t out[64];
    const uint8_t expected[64] = {
        0xa6, 0x9f, 0x73, 0xcc, 0xa2, 0x3a, 0x9a, 0xc5,
        0xc8, 0xb5, 0x67, 0xdc, 0x18, 0x5a, 0x75, 0x6e,
        0x97, 0xc9, 0x82, 0x16, 0x4f, 0xe2, 0x58, 0x59,
        0xe0, 0xd1, 0xdc, 0xc1, 0x47, 0x5c, 0x80, 0xa6,
        0x15, 0xb2, 0x12, 0x3a, 0xf1, 0xf5, 0xf9, 0x4c,
        0x11, 0xe3, 0xe9, 0x40, 0x2c, 0x3a, 0xc5, 0x58,
        0xf5, 0x00, 0x19, 0x9d, 0x95, 0xb6, 0xd3, 0xe3,
        0x01, 0x75, 0x85, 0x86, 0x28, 0x1d, 0xcd, 0x26
    };

    cfx_sha3_512(out, NULL, 0);
    CFX_ASSERT(memcmp(out, expected, 64) == 0);
}

/* SHA3-512("abc") */
static void test_sha3_512_abc(void) {
    uint8_t out[64];
    const uint8_t expected[64] = {
        0xb7, 0x51, 0x85, 0x0b, 0x1a, 0x57, 0x16, 0x8a,
        0x56, 0x93, 0xcd, 0x92, 0x4b, 0x6b, 0x09, 0x6e,
        0x08, 0xf6, 0x21, 0x82, 0x74, 0x44, 0xf7, 0x0d,
        0x88, 0x4f, 0x5d, 0x02, 0x40, 0xd2, 0x71, 0x2e,
        0x10, 0xe1, 0x16, 0xe9, 0x19, 0x2a, 0xf3, 0xc9,
        0x1a, 0x7e, 0xc5, 0x76, 0x47, 0xe3, 0x93, 0x40,
        0x57, 0x34, 0x0b, 0x4c, 0xf4, 0x08, 0xd5, 0xa5,
        0x65, 0x92, 0xf8, 0x27, 0x4e, 0xec, 0x53, 0xf0
    };

    cfx_sha3_512(out, "abc", 3);
    CFX_ASSERT(memcmp(out, expected, 64) == 0);
}

/* SHAKE128("", 32) */
static void test_shake128_empty(void) {
    uint8_t out[32];
    const uint8_t expected[32] = {
        0x7f, 0x9c, 0x2b, 0xa4, 0xe8, 0x8f, 0x82, 0x7d,
        0x61, 0x60, 0x45, 0x50, 0x76, 0x05, 0x85, 0x3e,
        0xd7, 0x3b, 0x80, 0x93, 0xf6, 0xef, 0xbc, 0x88,
        0xeb, 0x1a, 0x6e, 0xac, 0xfa, 0x66, 0xef, 0x26
    };

    cfx_shake128(out, 32, NULL, 0);
    CFX_ASSERT(memcmp(out, expected, 32) == 0);
}

/* SHAKE256("", 64) */
static void test_shake256_empty(void) {
    uint8_t out[64];
    const uint8_t expected[64] = {
        0x46, 0xb9, 0xdd, 0x2b, 0x0b, 0xa8, 0x8d, 0x13,
        0x23, 0x3b, 0x3f, 0xeb, 0x74, 0x3e, 0xeb, 0x24,
        0x3f, 0xcd, 0x52, 0xea, 0x62, 0xb8, 0x1b, 0x82,
        0xb5, 0x0c, 0x27, 0x64, 0x6e, 0xd5, 0x76, 0x2f,
        0xd7, 0x5d, 0xc4, 0xdd, 0xd8, 0xc0, 0xf2, 0x00,
        0xcb, 0x05, 0x01, 0x9d, 0x67, 0xb5, 0x92, 0xf6,
        0xfc, 0x82, 0x1c, 0x49, 0x47, 0x9a, 0xb4, 0x86,
        0x40, 0x29, 0x2e, 0xac, 0xb3, 0xb7, 0xc4, 0xbe
    };

    cfx_shake256(out, 64, NULL, 0);
    CFX_ASSERT(memcmp(out, expected, 64) == 0);
}

/* streaming API */
static void test_sha3_256_streaming(void) {
    cfx_sha3_ctx_t ctx;
    uint8_t out_oneshot[32], out_stream[32];
    const char* msg = "The quick brown fox jumps over the lazy dog";

    cfx_sha3_256(out_oneshot, msg, strlen(msg));

    cfx_sha3_256_init(&ctx);
    cfx_sha3_256_update(&ctx, msg, 10);
    cfx_sha3_256_update(&ctx, msg + 10, 20);
    cfx_sha3_256_update(&ctx, msg + 30, strlen(msg) - 30);
    cfx_sha3_256_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 32) == 0);
}

/* SHAKE extendable output - multiple squeezes */
static void test_shake128_multiple_squeeze(void) {
    cfx_sha3_ctx_t ctx1, ctx2;
    uint8_t out1[64], out2a[32], out2b[32];

    /* one-shot 64 bytes */
    cfx_shake128(out1, 64, "test", 4);

    /* streaming: squeeze 32 + 32 */
    cfx_shake128_init(&ctx2);
    cfx_shake128_update(&ctx2, "test", 4);
    cfx_shake128_final(&ctx2);
    cfx_shake128_squeeze(&ctx2, out2a, 32);
    cfx_shake128_squeeze(&ctx2, out2b, 32);

    /* should match */
    CFX_ASSERT(memcmp(out1, out2a, 32) == 0);
    CFX_ASSERT(memcmp(out1 + 32, out2b, 32) == 0);
}

/* SHAKE variable output */
static void test_shake256_variable_output(void) {
    uint8_t out16[16], out32[32], out64[64];

    cfx_shake256(out16, 16, "test", 4);
    cfx_shake256(out32, 32, "test", 4);
    cfx_shake256(out64, 64, "test", 4);

    /* first 16 bytes should match */
    CFX_ASSERT(memcmp(out16, out32, 16) == 0);
    CFX_ASSERT(memcmp(out16, out64, 16) == 0);
    CFX_ASSERT(memcmp(out32, out64, 32) == 0);
}

/* longer message (NIST: 200 bytes of 0xa3) */
static void test_sha3_256_200_bytes(void) {
    uint8_t msg[200];
    uint8_t out[32];
    const uint8_t expected[32] = {
        0x79, 0xf3, 0x8a, 0xde, 0xc5, 0xc2, 0x03, 0x07,
        0xa9, 0x8e, 0xf7, 0x6e, 0x83, 0x24, 0xaf, 0xbf,
        0xd4, 0x6c, 0xfd, 0x81, 0xb2, 0x2e, 0x39, 0x73,
        0xc6, 0x5f, 0xa1, 0xbd, 0x9d, 0xe3, 0x17, 0x87
    };

    memset(msg, 0xa3, 200);
    cfx_sha3_256(out, msg, 200);
    CFX_ASSERT(memcmp(out, expected, 32) == 0);
}

/* different inputs produce different outputs */
static void test_sha3_256_different_inputs(void) {
    uint8_t out1[32], out2[32];

    cfx_sha3_256(out1, "test1", 5);
    cfx_sha3_256(out2, "test2", 5);

    CFX_ASSERT(memcmp(out1, out2, 32) != 0);
}

/* byte-by-byte */
static void test_sha3_256_byte_by_byte(void) {
    cfx_sha3_ctx_t ctx;
    uint8_t out_oneshot[32], out_stream[32];
    const char* msg = "hello world";
    size_t len = strlen(msg);

    cfx_sha3_256(out_oneshot, msg, len);

    cfx_sha3_256_init(&ctx);
    for (size_t i = 0; i < len; i++) {
        cfx_sha3_256_update(&ctx, msg + i, 1);
    }
    cfx_sha3_256_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 32) == 0);
}

/* all variants produce different outputs for same input */
static void test_sha3_variants_differ(void) {
    uint8_t out224[28], out256[32], out384[48], out512[64];

    cfx_sha3_224(out224, "test", 4);
    cfx_sha3_256(out256, "test", 4);
    cfx_sha3_384(out384, "test", 4);
    cfx_sha3_512(out512, "test", 4);

    /* different sizes, but even truncated they should differ */
    CFX_ASSERT(memcmp(out224, out256, 28) != 0);
    CFX_ASSERT(memcmp(out256, out384, 32) != 0);
    CFX_ASSERT(memcmp(out384, out512, 48) != 0);
}

/* SHAKE vs SHA3 differ (different domain separator) */
static void test_shake_vs_sha3_differ(void) {
    uint8_t sha3_out[32], shake_out[32];

    cfx_sha3_256(sha3_out, "test", 4);
    cfx_shake256(shake_out, 32, "test", 4);

    CFX_ASSERT(memcmp(sha3_out, shake_out, 32) != 0);
}

int main(void) {
    CFX_TEST(test_sha3_256_empty);
    CFX_TEST(test_sha3_256_abc);
    CFX_TEST(test_sha3_224_empty);
    CFX_TEST(test_sha3_384_empty);
    CFX_TEST(test_sha3_512_empty);
    CFX_TEST(test_sha3_512_abc);
    CFX_TEST(test_shake128_empty);
    CFX_TEST(test_shake256_empty);
    CFX_TEST(test_sha3_256_streaming);
    CFX_TEST(test_shake128_multiple_squeeze);
    CFX_TEST(test_shake256_variable_output);
    CFX_TEST(test_sha3_256_200_bytes);
    CFX_TEST(test_sha3_256_different_inputs);
    CFX_TEST(test_sha3_256_byte_by_byte);
    CFX_TEST(test_sha3_variants_differ);
    CFX_TEST(test_shake_vs_sha3_differ);

    printf("\nAll SHA-3 tests passed!\n");
    return 0;
}
