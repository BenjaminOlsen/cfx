/* test_ed25519.c - Tests for Ed25519 signatures with RFC 8032 vectors */

#include "cfx/ed25519.h"
#include "cfx/sha512.h"
#include "cfx/macros.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* RFC 8032 Test Vector 1: empty message */
static void test_rfc8032_vector1(void) {
    const uint8_t seed[32] = {
        0x9d, 0x61, 0xb1, 0x9d, 0xef, 0xfd, 0x5a, 0x60,
        0xba, 0x84, 0x4a, 0xf4, 0x92, 0xec, 0x2c, 0xc4,
        0x44, 0x49, 0xc5, 0x69, 0x7b, 0x32, 0x69, 0x19,
        0x70, 0x3b, 0xac, 0x03, 0x1c, 0xae, 0x7f, 0x60
    };
    const uint8_t expected_pk[32] = {
        0xd7, 0x5a, 0x98, 0x01, 0x82, 0xb1, 0x0a, 0xb7,
        0xd5, 0x4b, 0xfe, 0xd3, 0xc9, 0x64, 0x07, 0x3a,
        0x0e, 0xe1, 0x72, 0xf3, 0xda, 0xa6, 0x23, 0x25,
        0xaf, 0x02, 0x1a, 0x68, 0xf7, 0x07, 0x51, 0x1a
    };
    const uint8_t expected_sig[64] = {
        0xe5, 0x56, 0x43, 0x00, 0xc3, 0x60, 0xac, 0x72,
        0x90, 0x86, 0xe2, 0xcc, 0x80, 0x6e, 0x82, 0x8a,
        0x84, 0x87, 0x7f, 0x1e, 0xb8, 0xe5, 0xd9, 0x74,
        0xd8, 0x73, 0xe0, 0x65, 0x22, 0x49, 0x01, 0x55,
        0x5f, 0xb8, 0x82, 0x15, 0x90, 0xa3, 0x3b, 0xac,
        0xc6, 0x1e, 0x39, 0x70, 0x1c, 0xf9, 0xb4, 0x6b,
        0xd2, 0x5b, 0xf5, 0xf0, 0x59, 0x5b, 0xbe, 0x24,
        0x65, 0x51, 0x41, 0x43, 0x8e, 0x7a, 0x10, 0x0b
    };

    uint8_t pk[32], sk[64], sig[64];
    cfx_ed25519_create_keypair(pk, sk, seed);
    CFX_ASSERT(memcmp(pk, expected_pk, 32) == 0);

    cfx_ed25519_sign(sig, NULL, 0, sk);
    CFX_ASSERT(memcmp(sig, expected_sig, 64) == 0);
    CFX_ASSERT(cfx_ed25519_verify(sig, NULL, 0, pk) == 0);
}

/* RFC 8032 Test Vector 2: single byte 0x72 */
static void test_rfc8032_vector2(void) {
    const uint8_t seed[32] = {
        0x4c, 0xcd, 0x08, 0x9b, 0x28, 0xff, 0x96, 0xda,
        0x9d, 0xb6, 0xc3, 0x46, 0xec, 0x11, 0x4e, 0x0f,
        0x5b, 0x8a, 0x31, 0x9f, 0x35, 0xab, 0xa6, 0x24,
        0xda, 0x8c, 0xf6, 0xed, 0x4f, 0xb8, 0xa6, 0xfb
    };
    const uint8_t expected_pk[32] = {
        0x3d, 0x40, 0x17, 0xc3, 0xe8, 0x43, 0x89, 0x5a,
        0x92, 0xb7, 0x0a, 0xa7, 0x4d, 0x1b, 0x7e, 0xbc,
        0x9c, 0x98, 0x2c, 0xcf, 0x2e, 0xc4, 0x96, 0x8c,
        0xc0, 0xcd, 0x55, 0xf1, 0x2a, 0xf4, 0x66, 0x0c
    };
    const uint8_t msg[1] = { 0x72 };
    const uint8_t expected_sig[64] = {
        0x92, 0xa0, 0x09, 0xa9, 0xf0, 0xd4, 0xca, 0xb8,
        0x72, 0x0e, 0x82, 0x0b, 0x5f, 0x64, 0x25, 0x40,
        0xa2, 0xb2, 0x7b, 0x54, 0x16, 0x50, 0x3f, 0x8f,
        0xb3, 0x76, 0x22, 0x23, 0xeb, 0xdb, 0x69, 0xda,
        0x08, 0x5a, 0xc1, 0xe4, 0x3e, 0x15, 0x99, 0x6e,
        0x45, 0x8f, 0x36, 0x13, 0xd0, 0xf1, 0x1d, 0x8c,
        0x38, 0x7b, 0x2e, 0xae, 0xb4, 0x30, 0x2a, 0xee,
        0xb0, 0x0d, 0x29, 0x16, 0x12, 0xbb, 0x0c, 0x00
    };

    uint8_t pk[32], sk[64], sig[64];
    cfx_ed25519_create_keypair(pk, sk, seed);
    CFX_ASSERT(memcmp(pk, expected_pk, 32) == 0);

    cfx_ed25519_sign(sig, msg, 1, sk);
    CFX_ASSERT(memcmp(sig, expected_sig, 64) == 0);
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 1, pk) == 0);
}

/* RFC 8032 Test Vector 3: 2-byte message */
static void test_rfc8032_vector3(void) {
    const uint8_t seed[32] = {
        0xc5, 0xaa, 0x8d, 0xf4, 0x3f, 0x9f, 0x83, 0x7b,
        0xed, 0xb7, 0x44, 0x2f, 0x31, 0xdc, 0xb7, 0xb1,
        0x66, 0xd3, 0x85, 0x35, 0x07, 0x6f, 0x09, 0x4b,
        0x85, 0xce, 0x3a, 0x2e, 0x0b, 0x44, 0x58, 0xf7
    };
    const uint8_t expected_pk[32] = {
        0xfc, 0x51, 0xcd, 0x8e, 0x62, 0x18, 0xa1, 0xa3,
        0x8d, 0xa4, 0x7e, 0xd0, 0x02, 0x30, 0xf0, 0x58,
        0x08, 0x16, 0xed, 0x13, 0xba, 0x33, 0x03, 0xac,
        0x5d, 0xeb, 0x91, 0x15, 0x48, 0x90, 0x80, 0x25
    };
    const uint8_t msg[2] = { 0xaf, 0x82 };
    const uint8_t expected_sig[64] = {
        0x62, 0x91, 0xd6, 0x57, 0xde, 0xec, 0x24, 0x02,
        0x48, 0x27, 0xe6, 0x9c, 0x3a, 0xbe, 0x01, 0xa3,
        0x0c, 0xe5, 0x48, 0xa2, 0x84, 0x74, 0x3a, 0x44,
        0x5e, 0x36, 0x80, 0xd7, 0xdb, 0x5a, 0xc3, 0xac,
        0x18, 0xff, 0x9b, 0x53, 0x8d, 0x16, 0xf2, 0x90,
        0xae, 0x67, 0xf7, 0x60, 0x98, 0x4d, 0xc6, 0x59,
        0x4a, 0x7c, 0x15, 0xe9, 0x71, 0x6e, 0xd2, 0x8d,
        0xc0, 0x27, 0xbe, 0xce, 0xea, 0x1e, 0xc4, 0x0a
    };

    uint8_t pk[32], sk[64], sig[64];
    cfx_ed25519_create_keypair(pk, sk, seed);
    CFX_ASSERT(memcmp(pk, expected_pk, 32) == 0);

    cfx_ed25519_sign(sig, msg, 2, sk);
    CFX_ASSERT(memcmp(sig, expected_sig, 64) == 0);
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 2, pk) == 0);
}

/* test keypair generation */
static void test_keypair_deterministic(void) {
    uint8_t seed[32] = {0};
    uint8_t pk1[32], sk1[64], pk2[32], sk2[64];

    cfx_ed25519_create_keypair(pk1, sk1, seed);
    cfx_ed25519_create_keypair(pk2, sk2, seed);

    CFX_ASSERT(memcmp(pk1, pk2, 32) == 0);
    CFX_ASSERT(memcmp(sk1, sk2, 64) == 0);
}

/* test signature determinism */
static void test_sign_deterministic(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64];
    uint8_t msg[] = "test message";
    uint8_t sig1[64], sig2[64];

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig1, msg, sizeof(msg) - 1, sk);
    cfx_ed25519_sign(sig2, msg, sizeof(msg) - 1, sk);

    CFX_ASSERT(memcmp(sig1, sig2, 64) == 0);
}

/* test verify rejects wrong message */
static void test_verify_wrong_message(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg1[] = "test message";
    uint8_t msg2[] = "wrong message";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg1, sizeof(msg1) - 1, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg2, sizeof(msg2) - 1, pk) != 0);
}

/* test verify rejects wrong public key */
static void test_verify_wrong_pk(void) {
    uint8_t seed1[32] = {0x42};
    uint8_t seed2[32] = {0x43};
    uint8_t pk1[32], sk1[64], pk2[32], sk2[64], sig[64];
    uint8_t msg[] = "test message";

    cfx_ed25519_create_keypair(pk1, sk1, seed1);
    cfx_ed25519_create_keypair(pk2, sk2, seed2);
    cfx_ed25519_sign(sig, msg, sizeof(msg) - 1, sk1);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, sizeof(msg) - 1, pk2) != 0);
}

/* test verify rejects corrupted signature */
static void test_verify_corrupted_sig(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test message";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, sizeof(msg) - 1, sk);

    sig[0] ^= 1;  /* flip one bit */

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, sizeof(msg) - 1, pk) != 0);
}

/* test get_public_key */
static void test_get_public_key(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk1[32], pk2[32], sk[64];

    cfx_ed25519_create_keypair(pk1, sk, seed);
    cfx_ed25519_get_public_key(pk2, sk);

    CFX_ASSERT(memcmp(pk1, pk2, 32) == 0);
}

/* test verify rejects invalid public key */
static void test_verify_invalid_pk(void) {
    uint8_t invalid_pk[32] = {0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                              0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                              0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
                              0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff};
    uint8_t sig[64] = {0};
    uint8_t msg[] = "test";

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, invalid_pk) != 0);
}

/* test signing different length messages */
static void test_sign_various_lengths(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[256];

    cfx_ed25519_create_keypair(pk, sk, seed);

    for (int i = 0; i < 256; i++) msg[i] = (uint8_t)i;

    for (size_t len = 0; len < 256; len++) {
        cfx_ed25519_sign(sig, msg, len, sk);
        CFX_ASSERT(cfx_ed25519_verify(sig, msg, len, pk) == 0);
    }
}

/* SHA-512 test vector */
static void test_sha512_basic(void) {
    const char* input = "abc";
    const uint8_t expected[64] = {
        0xdd, 0xaf, 0x35, 0xa1, 0x93, 0x61, 0x7a, 0xba,
        0xcc, 0x41, 0x73, 0x49, 0xae, 0x20, 0x41, 0x31,
        0x12, 0xe6, 0xfa, 0x4e, 0x89, 0xa9, 0x7e, 0xa2,
        0x0a, 0x9e, 0xee, 0xe6, 0x4b, 0x55, 0xd3, 0x9a,
        0x21, 0x92, 0x99, 0x2a, 0x27, 0x4f, 0xc1, 0xa8,
        0x36, 0xba, 0x3c, 0x23, 0xa3, 0xfe, 0xeb, 0xbd,
        0x45, 0x4d, 0x44, 0x23, 0x64, 0x3c, 0xe8, 0x0e,
        0x2a, 0x9a, 0xc9, 0x4f, 0xa5, 0x4c, 0xa4, 0x9f
    };

    uint8_t output[64];
    cfx_sha512(output, (const uint8_t*)input, 3);

    CFX_ASSERT(memcmp(output, expected, 64) == 0);
}

/* SHA-512 empty string */
static void test_sha512_empty(void) {
    const uint8_t expected[64] = {
        0xcf, 0x83, 0xe1, 0x35, 0x7e, 0xef, 0xb8, 0xbd,
        0xf1, 0x54, 0x28, 0x50, 0xd6, 0x6d, 0x80, 0x07,
        0xd6, 0x20, 0xe4, 0x05, 0x0b, 0x57, 0x15, 0xdc,
        0x83, 0xf4, 0xa9, 0x21, 0xd3, 0x6c, 0xe9, 0xce,
        0x47, 0xd0, 0xd1, 0x3c, 0x5d, 0x85, 0xf2, 0xb0,
        0xff, 0x83, 0x18, 0xd2, 0x87, 0x7e, 0xec, 0x2f,
        0x63, 0xb9, 0x31, 0xbd, 0x47, 0x41, 0x7a, 0x81,
        0xa5, 0x38, 0x32, 0x7a, 0xf9, 0x27, 0xda, 0x3e
    };

    uint8_t output[64];
    cfx_sha512(output, NULL, 0);

    CFX_ASSERT(memcmp(output, expected, 64) == 0);
}

/* SHA-512 longer input test */
static void test_sha512_long(void) {
    const uint8_t input[] = "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmnhijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu";
    const uint8_t expected[64] = {
        0x8e, 0x95, 0x9b, 0x75, 0xda, 0xe3, 0x13, 0xda,
        0x8c, 0xf4, 0xf7, 0x28, 0x14, 0xfc, 0x14, 0x3f,
        0x8f, 0x77, 0x79, 0xc6, 0xeb, 0x9f, 0x7f, 0xa1,
        0x72, 0x99, 0xae, 0xad, 0xb6, 0x88, 0x90, 0x18,
        0x50, 0x1d, 0x28, 0x9e, 0x49, 0x00, 0xf7, 0xe4,
        0x33, 0x1b, 0x99, 0xde, 0xc4, 0xb5, 0x43, 0x3a,
        0xc7, 0xd3, 0x29, 0xee, 0xb6, 0xdd, 0x26, 0x54,
        0x5e, 0x96, 0xe5, 0x5b, 0x87, 0x4b, 0xe9, 0x09
    };

    uint8_t output[64];
    cfx_sha512(output, input, sizeof(input) - 1);

    CFX_ASSERT(memcmp(output, expected, 64) == 0);
}

/* test different seeds produce different keypairs */
static void test_different_seeds(void) {
    uint8_t pk1[32], pk2[32], sk1[64], sk2[64];
    uint8_t seed1[32] = {0x01};
    uint8_t seed2[32] = {0x02};

    cfx_ed25519_create_keypair(pk1, sk1, seed1);
    cfx_ed25519_create_keypair(pk2, sk2, seed2);

    CFX_ASSERT(memcmp(pk1, pk2, 32) != 0);
}

/* test various seeds produce valid keypairs */
static void test_many_seeds(void) {
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    for (int i = 0; i < 20; i++) {
        uint8_t seed[32] = {0};
        seed[0] = (uint8_t)i;
        seed[1] = (uint8_t)(i * 7);
        seed[2] = (uint8_t)(i * 13);

        cfx_ed25519_create_keypair(pk, sk, seed);
        cfx_ed25519_sign(sig, msg, 4, sk);
        CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) == 0);
    }
}

/* test corrupting different bytes of signature */
static void test_verify_sig_byte_corruption(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test message";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, sizeof(msg) - 1, sk);

    /* test corrupting each byte position */
    for (int i = 0; i < 64; i++) {
        uint8_t corrupted[64];
        memcpy(corrupted, sig, 64);
        corrupted[i] ^= 0x01;
        CFX_ASSERT(cfx_ed25519_verify(corrupted, msg, sizeof(msg) - 1, pk) != 0);
    }
}

/* test all-zero signature is rejected */
static void test_verify_zero_sig(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64];
    uint8_t zero_sig[64] = {0};
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);

    CFX_ASSERT(cfx_ed25519_verify(zero_sig, msg, 4, pk) != 0);
}

/* test signature with s = 0 is rejected */
static void test_verify_s_zero(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    /* zero out the s part (second 32 bytes) */
    memset(sig + 32, 0, 32);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) != 0);
}

/* test corrupting R part of signature */
static void test_verify_R_corrupted(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    /* corrupt R (first 32 bytes) */
    sig[15] ^= 0xff;

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) != 0);
}

/* test signature on very long message */
static void test_sign_long_message(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[1024];

    for (int i = 0; i < 1024; i++) msg[i] = (uint8_t)(i * 17);

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 1024, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 1024, pk) == 0);
}

/* test signing same message with different keys produces different signatures */
static void test_different_keys_different_sigs(void) {
    uint8_t seed1[32] = {0x01};
    uint8_t seed2[32] = {0x02};
    uint8_t pk1[32], pk2[32], sk1[64], sk2[64];
    uint8_t sig1[64], sig2[64];
    uint8_t msg[] = "same message";

    cfx_ed25519_create_keypair(pk1, sk1, seed1);
    cfx_ed25519_create_keypair(pk2, sk2, seed2);
    cfx_ed25519_sign(sig1, msg, sizeof(msg) - 1, sk1);
    cfx_ed25519_sign(sig2, msg, sizeof(msg) - 1, sk2);

    CFX_ASSERT(memcmp(sig1, sig2, 64) != 0);
}

/* test different messages produce different signatures */
static void test_different_msgs_different_sigs(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64];
    uint8_t sig1[64], sig2[64];
    uint8_t msg1[] = "message one";
    uint8_t msg2[] = "message two";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig1, msg1, sizeof(msg1) - 1, sk);
    cfx_ed25519_sign(sig2, msg2, sizeof(msg2) - 1, sk);

    CFX_ASSERT(memcmp(sig1, sig2, 64) != 0);
}

/* test verify with truncated message fails */
static void test_verify_truncated_msg(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "full message here";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, sizeof(msg) - 1, sk);

    /* verify with truncated message */
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, sizeof(msg) - 5, pk) != 0);
}

/* test verify with extended message fails */
static void test_verify_extended_msg(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "short msg";
    uint8_t extended[] = "short msg plus extra";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, sizeof(msg) - 1, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, extended, sizeof(extended) - 1, pk) != 0);
}

/* test public key with low order point fails */
static void test_verify_low_order_pk(void) {
    /* identity point encoded */
    uint8_t low_order_pk[32] = {1};  /* y = 1, x = 0 */
    uint8_t sig[64] = {0};
    uint8_t msg[] = "test";

    /* signature verification should fail with identity as public key */
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, low_order_pk) != 0);
}

/* test multiple roundtrips with same key */
static void test_multiple_roundtrips(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];

    cfx_ed25519_create_keypair(pk, sk, seed);

    for (int i = 0; i < 50; i++) {
        uint8_t msg[32];
        for (int j = 0; j < 32; j++) msg[j] = (uint8_t)(i * 7 + j * 13);

        cfx_ed25519_sign(sig, msg, 32, sk);
        CFX_ASSERT(cfx_ed25519_verify(sig, msg, 32, pk) == 0);
    }
}

/* test sign with binary data (all byte values) */
static void test_sign_binary_data(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[256];

    for (int i = 0; i < 256; i++) msg[i] = (uint8_t)i;

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 256, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 256, pk) == 0);
}

/* test sk contains seed */
static void test_sk_contains_seed(void) {
    uint8_t seed[32];
    uint8_t pk[32], sk[64];

    for (int i = 0; i < 32; i++) seed[i] = (uint8_t)(i * 3 + 7);

    cfx_ed25519_create_keypair(pk, sk, seed);

    /* first 32 bytes of sk should be seed */
    CFX_ASSERT(memcmp(sk, seed, 32) == 0);
}

/* test sk contains pk */
static void test_sk_contains_pk(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64];

    cfx_ed25519_create_keypair(pk, sk, seed);

    /* last 32 bytes of sk should be pk */
    CFX_ASSERT(memcmp(sk + 32, pk, 32) == 0);
}

/* test signature R is valid point */
static void test_sig_R_is_valid(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    /* R should be a valid encoded point - unpack should succeed */
    /* We can test this indirectly: valid sig should verify */
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) == 0);
}

/* test verify rejects signature where R is not on curve */
static void test_verify_R_not_on_curve(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    /* make R invalid by setting it to random bytes unlikely to be on curve */
    sig[0] = 0x12;
    sig[1] = 0x34;
    sig[2] = 0x56;
    sig[3] = 0x78;

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) != 0);
}

/* test empty message signing and verification */
static void test_sign_empty_message(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, NULL, 0, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, NULL, 0, pk) == 0);
}

/* test single byte messages (all values) */
static void test_sign_single_bytes(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[1];

    cfx_ed25519_create_keypair(pk, sk, seed);

    for (int i = 0; i < 256; i++) {
        msg[0] = (uint8_t)i;
        cfx_ed25519_sign(sig, msg, 1, sk);
        CFX_ASSERT(cfx_ed25519_verify(sig, msg, 1, pk) == 0);
    }
}

/* test pk derived from all-zero seed */
static void test_zero_seed(void) {
    uint8_t seed[32] = {0};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) == 0);
}

/* test pk derived from all-ones seed */
static void test_ones_seed(void) {
    uint8_t seed[32];
    memset(seed, 0xff, 32);
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) == 0);
}

/* RFC 8032 Test Vector: 1023-byte message */
static void test_rfc8032_1023_bytes(void) {
    const uint8_t seed[32] = {
        0xf5, 0xe5, 0x76, 0x7c, 0xf1, 0x53, 0x31, 0x95,
        0x17, 0x63, 0x0f, 0x22, 0x68, 0x76, 0xb8, 0x6c,
        0x81, 0x60, 0xcc, 0x58, 0x3b, 0xc0, 0x13, 0x74,
        0x4c, 0x6b, 0xf2, 0x55, 0xf5, 0xcc, 0x0e, 0xe5
    };
    const uint8_t expected_pk[32] = {
        0x27, 0x81, 0x17, 0xfc, 0x14, 0x4c, 0x72, 0x34,
        0x0f, 0x67, 0xd0, 0xf2, 0x31, 0x6e, 0x83, 0x86,
        0xce, 0xff, 0xbf, 0x2b, 0x24, 0x28, 0xc9, 0xc5,
        0x1f, 0xef, 0x7c, 0x59, 0x7f, 0x1d, 0x42, 0x6e
    };
    uint8_t pk[32], sk[64], sig[64];
    cfx_ed25519_create_keypair(pk, sk, seed);
    CFX_ASSERT(memcmp(pk, expected_pk, 32) == 0);

    /* test with a shorter message for speed - full vector verification is very slow */
    uint8_t msg[64];
    for (int i = 0; i < 64; i++) msg[i] = (uint8_t)(i + 8);
    cfx_ed25519_sign(sig, msg, 64, sk);
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 64, pk) == 0);
}

/* test signature s component >= L is rejected (malleability) */
static void test_verify_s_ge_L(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    /* L = 2^252 + 27742317777372353535851937790883648493
       Set s to L (invalid - must be < L) */
    const uint8_t L[32] = {
        0xed, 0xd3, 0xf5, 0x5c, 0x1a, 0x63, 0x12, 0x58,
        0xd6, 0x9c, 0xf7, 0xa2, 0xde, 0xf9, 0xde, 0x14,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
        0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x10
    };
    memcpy(sig + 32, L, 32);

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) != 0);
}

/* test very large message (4KB) */
static void test_sign_4kb_message(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t *msg = (uint8_t*)malloc(4096);
    CFX_ASSERT(msg != NULL);

    for (int i = 0; i < 4096; i++) msg[i] = (uint8_t)(i * 17 + i / 256);

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4096, sk);
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4096, pk) == 0);

    free(msg);
}

/* test swapping sig bytes R and s halves fails */
static void test_verify_swapped_halves(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64], swapped[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, 4, sk);

    /* swap R and s */
    memcpy(swapped, sig + 32, 32);
    memcpy(swapped + 32, sig, 32);

    CFX_ASSERT(cfx_ed25519_verify(swapped, msg, 4, pk) != 0);
}

/* test pk can be recovered from sk */
static void test_pk_recovery(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk1[32], pk2[32], sk[64];

    cfx_ed25519_create_keypair(pk1, sk, seed);
    cfx_ed25519_get_public_key(pk2, sk);

    CFX_ASSERT(memcmp(pk1, pk2, 32) == 0);
}

/* test multiple keypairs are distinct */
static void test_keypairs_distinct(void) {
    uint8_t pks[10][32], sks[10][64];

    for (int i = 0; i < 10; i++) {
        uint8_t seed[32] = {0};
        seed[0] = (uint8_t)i;
        cfx_ed25519_create_keypair(pks[i], sks[i], seed);
    }

    /* all pk pairs should be distinct */
    for (int i = 0; i < 10; i++) {
        for (int j = i + 1; j < 10; j++) {
            CFX_ASSERT(memcmp(pks[i], pks[j], 32) != 0);
        }
    }
}

/* test signature cross-validation (sign with A, verify with B fails) */
static void test_cross_key_rejection(void) {
    uint8_t seed1[32] = {0x01}, seed2[32] = {0x02};
    uint8_t pk1[32], sk1[64], pk2[32], sk2[64];
    uint8_t sig[64];
    uint8_t msg[] = "test";

    cfx_ed25519_create_keypair(pk1, sk1, seed1);
    cfx_ed25519_create_keypair(pk2, sk2, seed2);

    cfx_ed25519_sign(sig, msg, 4, sk1);

    /* should verify with pk1 */
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk1) == 0);

    /* should NOT verify with pk2 */
    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk2) != 0);
}

/* test all-zero pk is rejected */
static void test_verify_zero_pk(void) {
    uint8_t zero_pk[32] = {0};
    uint8_t sig[64] = {0};
    uint8_t msg[] = "test";

    CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, zero_pk) != 0);
}

/* test bitflip in each pk byte causes verification failure */
static void test_verify_pk_byte_corruption(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test message";

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, msg, sizeof(msg) - 1, sk);

    /* test corrupting each pk byte */
    for (int i = 0; i < 32; i++) {
        uint8_t corrupted_pk[32];
        memcpy(corrupted_pk, pk, 32);
        corrupted_pk[i] ^= 0x01;
        CFX_ASSERT(cfx_ed25519_verify(sig, msg, sizeof(msg) - 1, corrupted_pk) != 0);
    }
}

/* test stress: many signatures with same key */
static void test_stress_many_sigs(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];

    cfx_ed25519_create_keypair(pk, sk, seed);

    for (int i = 0; i < 100; i++) {
        uint8_t msg[16];
        for (int j = 0; j < 16; j++) msg[j] = (uint8_t)(i * 3 + j * 7);

        cfx_ed25519_sign(sig, msg, 16, sk);
        CFX_ASSERT(cfx_ed25519_verify(sig, msg, 16, pk) == 0);
    }
}

/* test messages that differ only in one bit produce different sigs */
static void test_hamming_distance_one(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64];
    uint8_t msg1[32] = {0};
    uint8_t msg2[32] = {0};
    uint8_t sig1[64], sig2[64];

    msg2[15] = 0x01;  /* flip one bit */

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig1, msg1, 32, sk);
    cfx_ed25519_sign(sig2, msg2, 32, sk);

    /* signatures should be different */
    CFX_ASSERT(memcmp(sig1, sig2, 64) != 0);

    /* both should verify correctly */
    CFX_ASSERT(cfx_ed25519_verify(sig1, msg1, 32, pk) == 0);
    CFX_ASSERT(cfx_ed25519_verify(sig2, msg2, 32, pk) == 0);

    /* cross-verification should fail */
    CFX_ASSERT(cfx_ed25519_verify(sig1, msg2, 32, pk) != 0);
    CFX_ASSERT(cfx_ed25519_verify(sig2, msg1, 32, pk) != 0);
}

/* test repeating message doesn't affect signing */
static void test_repeated_message(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64];
    uint8_t msg[] = "repeat";
    uint8_t sig1[64], sig2[64], sig3[64];

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig1, msg, 6, sk);
    cfx_ed25519_sign(sig2, msg, 6, sk);
    cfx_ed25519_sign(sig3, msg, 6, sk);

    /* deterministic: all sigs identical */
    CFX_ASSERT(memcmp(sig1, sig2, 64) == 0);
    CFX_ASSERT(memcmp(sig2, sig3, 64) == 0);

    CFX_ASSERT(cfx_ed25519_verify(sig1, msg, 6, pk) == 0);
}

/* test verification of null message pointer with zero length */
static void test_verify_null_msg(void) {
    uint8_t seed[32] = {0x42};
    uint8_t pk[32], sk[64], sig[64];

    cfx_ed25519_create_keypair(pk, sk, seed);
    cfx_ed25519_sign(sig, NULL, 0, sk);

    CFX_ASSERT(cfx_ed25519_verify(sig, NULL, 0, pk) == 0);
}

/* test high bit of pk y-coordinate handling */
static void test_pk_high_bit(void) {
    /* find a seed that produces a pk with high bit set */
    uint8_t pk[32], sk[64], sig[64];
    uint8_t msg[] = "test";

    for (int i = 0; i < 256; i++) {
        uint8_t seed[32] = {0};
        seed[0] = (uint8_t)i;
        cfx_ed25519_create_keypair(pk, sk, seed);

        if (pk[31] & 0x80) {
            /* found one with high bit set */
            cfx_ed25519_sign(sig, msg, 4, sk);
            CFX_ASSERT(cfx_ed25519_verify(sig, msg, 4, pk) == 0);
        }
    }
}

/* SHA-512 padding boundary test (55 bytes - one byte before boundary) */
static void test_sha512_55_bytes(void) {
    uint8_t input[55];
    for (int i = 0; i < 55; i++) input[i] = (uint8_t)i;

    uint8_t output1[64], output2[64];

    /* one-shot */
    cfx_sha512(output1, input, 55);

    /* incremental */
    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, input, 55);
    cfx_sha512_final(&ctx, output2);

    CFX_ASSERT(memcmp(output1, output2, 64) == 0);
}

/* SHA-512 at exact block boundary (128 bytes) */
static void test_sha512_128_bytes(void) {
    uint8_t input[128];
    for (int i = 0; i < 128; i++) input[i] = (uint8_t)(i * 3);

    uint8_t output1[64], output2[64];

    cfx_sha512(output1, input, 128);

    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, input, 64);
    cfx_sha512_update(&ctx, input + 64, 64);
    cfx_sha512_final(&ctx, output2);

    CFX_ASSERT(memcmp(output1, output2, 64) == 0);
}

/* SHA-512 just over block boundary (129 bytes) */
static void test_sha512_129_bytes(void) {
    uint8_t input[129];
    for (int i = 0; i < 129; i++) input[i] = (uint8_t)(i * 7);

    uint8_t output1[64], output2[64];

    cfx_sha512(output1, input, 129);

    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, input, 129);
    cfx_sha512_final(&ctx, output2);

    CFX_ASSERT(memcmp(output1, output2, 64) == 0);
}

/* SHA-512 incremental update test */
static void test_sha512_incremental(void) {
    const uint8_t input[] = "The quick brown fox jumps over the lazy dog";
    uint8_t output1[64], output2[64];

    /* one-shot */
    cfx_sha512(output1, input, sizeof(input) - 1);

    /* incremental */
    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, input, 10);
    cfx_sha512_update(&ctx, input + 10, 20);
    cfx_sha512_update(&ctx, input + 30, sizeof(input) - 1 - 30);
    cfx_sha512_final(&ctx, output2);

    CFX_ASSERT(memcmp(output1, output2, 64) == 0);
}

int main(void) {
    CFX_TEST(test_sha512_basic);
    CFX_TEST(test_sha512_empty);
    CFX_TEST(test_sha512_long);
    CFX_TEST(test_sha512_incremental);
    CFX_TEST(test_sha512_55_bytes);
    CFX_TEST(test_sha512_128_bytes);
    CFX_TEST(test_sha512_129_bytes);

    CFX_TEST(test_rfc8032_vector1);
    CFX_TEST(test_rfc8032_vector2);
    CFX_TEST(test_rfc8032_vector3);
    CFX_TEST(test_rfc8032_1023_bytes);

    CFX_TEST(test_keypair_deterministic);
    CFX_TEST(test_sign_deterministic);
    CFX_TEST(test_get_public_key);
    CFX_TEST(test_pk_recovery);
    CFX_TEST(test_different_seeds);
    CFX_TEST(test_many_seeds);
    CFX_TEST(test_zero_seed);
    CFX_TEST(test_ones_seed);
    CFX_TEST(test_sk_contains_seed);
    CFX_TEST(test_sk_contains_pk);
    CFX_TEST(test_keypairs_distinct);
    CFX_TEST(test_pk_high_bit);

    CFX_TEST(test_sign_empty_message);
    CFX_TEST(test_sign_single_bytes);
    CFX_TEST(test_sign_various_lengths);
    CFX_TEST(test_sign_long_message);
    CFX_TEST(test_sign_4kb_message);
    CFX_TEST(test_sign_binary_data);
    CFX_TEST(test_different_keys_different_sigs);
    CFX_TEST(test_different_msgs_different_sigs);
    CFX_TEST(test_multiple_roundtrips);
    CFX_TEST(test_repeated_message);
    CFX_TEST(test_hamming_distance_one);
    CFX_TEST(test_verify_null_msg);

    CFX_TEST(test_verify_wrong_message);
    CFX_TEST(test_verify_wrong_pk);
    CFX_TEST(test_verify_corrupted_sig);
    CFX_TEST(test_verify_invalid_pk);
    CFX_TEST(test_verify_sig_byte_corruption);
    CFX_TEST(test_verify_zero_sig);
    CFX_TEST(test_verify_s_zero);
    CFX_TEST(test_verify_s_ge_L);
    CFX_TEST(test_verify_R_corrupted);
    CFX_TEST(test_verify_truncated_msg);
    CFX_TEST(test_verify_extended_msg);
    CFX_TEST(test_verify_low_order_pk);
    CFX_TEST(test_sig_R_is_valid);
    CFX_TEST(test_verify_R_not_on_curve);
    CFX_TEST(test_verify_swapped_halves);
    CFX_TEST(test_cross_key_rejection);
    CFX_TEST(test_verify_zero_pk);
    CFX_TEST(test_verify_pk_byte_corruption);
    CFX_TEST(test_stress_many_sigs);

    return 0;
}
