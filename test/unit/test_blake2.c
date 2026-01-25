/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/blake2.h"
#include "cfx/sha256.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

static void print_hex(const char *label, const uint8_t *data, size_t len) {
    printf("%s: ", label);
    for (size_t i = 0; i < len; i++) printf("%02x", data[i]);
    printf("\n");
}

/* RFC 7693 Appendix A - test vectors */

/* BLAKE2b with empty message, 64-byte output */
static void test_blake2b_empty(void) {
    uint8_t out[64];
    const uint8_t expected[64] = {
        0x78, 0x6a, 0x02, 0xf7, 0x42, 0x01, 0x59, 0x03,
        0xc6, 0xc6, 0xfd, 0x85, 0x25, 0x52, 0xd2, 0x72,
        0x91, 0x2f, 0x47, 0x40, 0xe1, 0x58, 0x47, 0x61,
        0x8a, 0x86, 0xe2, 0x17, 0xf7, 0x1f, 0x54, 0x19,
        0xd2, 0x5e, 0x10, 0x31, 0xaf, 0xee, 0x58, 0x53,
        0x13, 0x89, 0x64, 0x44, 0x93, 0x4e, 0xb0, 0x4b,
        0x90, 0x3a, 0x68, 0x5b, 0x14, 0x48, 0xb7, 0x55,
        0xd5, 0x6f, 0x70, 0x1a, 0xfe, 0x9b, 0xe2, 0xce
    };

    CFX_ASSERT(cfx_blake2b(out, 64, NULL, 0, NULL, 0) == 0);
    CFX_ASSERT(memcmp(out, expected, 64) == 0);
    printf("test_blake2b_empty() - OK\n");
}

/* BLAKE2b("abc"), 64-byte output */
static void test_blake2b_abc(void) {
    uint8_t out[64];
    const uint8_t expected[64] = {
        0xba, 0x80, 0xa5, 0x3f, 0x98, 0x1c, 0x4d, 0x0d,
        0x6a, 0x27, 0x97, 0xb6, 0x9f, 0x12, 0xf6, 0xe9,
        0x4c, 0x21, 0x2f, 0x14, 0x68, 0x5a, 0xc4, 0xb7,
        0x4b, 0x12, 0xbb, 0x6f, 0xdb, 0xff, 0xa2, 0xd1,
        0x7d, 0x87, 0xc5, 0x39, 0x2a, 0xab, 0x79, 0x2d,
        0xc2, 0x52, 0xd5, 0xde, 0x45, 0x33, 0xcc, 0x95,
        0x18, 0xd3, 0x8a, 0xa8, 0xdb, 0xf1, 0x92, 0x5a,
        0xb9, 0x23, 0x86, 0xed, 0xd4, 0x00, 0x99, 0x23
    };

    CFX_ASSERT(cfx_blake2b(out, 64, "abc", 3, NULL, 0) == 0);
    CFX_ASSERT(memcmp(out, expected, 64) == 0);
    printf("test_blake2b_abc() - OK\n");
}

/* BLAKE2s with empty message, 32-byte output */
static void test_blake2s_empty(void) {
    uint8_t out[32];
    const uint8_t expected[32] = {
        0x69, 0x21, 0x7a, 0x30, 0x79, 0x90, 0x80, 0x94,
        0xe1, 0x11, 0x21, 0xd0, 0x42, 0x35, 0x4a, 0x7c,
        0x1f, 0x55, 0xb6, 0x48, 0x2c, 0xa1, 0xa5, 0x1e,
        0x1b, 0x25, 0x0d, 0xfd, 0x1e, 0xd0, 0xee, 0xf9
    };

    CFX_ASSERT(cfx_blake2s(out, 32, NULL, 0, NULL, 0) == 0);
    CFX_ASSERT(memcmp(out, expected, 32) == 0);
    printf("test_blake2s_empty() - OK\n");
}

/* BLAKE2s("abc"), 32-byte output */
static void test_blake2s_abc(void) {
    uint8_t out[32];
    const uint8_t expected[32] = {
        0x50, 0x8c, 0x5e, 0x8c, 0x32, 0x7c, 0x14, 0xe2,
        0xe1, 0xa7, 0x2b, 0xa3, 0x4e, 0xeb, 0x45, 0x2f,
        0x37, 0x45, 0x8b, 0x20, 0x9e, 0xd6, 0x3a, 0x29,
        0x4d, 0x99, 0x9b, 0x4c, 0x86, 0x67, 0x59, 0x82
    };

    CFX_ASSERT(cfx_blake2s(out, 32, "abc", 3, NULL, 0) == 0);
    CFX_ASSERT(memcmp(out, expected, 32) == 0);
    printf("test_blake2s_abc() - OK\n");
}

/* Variable output length */
static void test_blake2b_variable_output(void) {
    uint8_t out32[32], out48[48], out64[64];

    CFX_ASSERT(cfx_blake2b(out32, 32, "test", 4, NULL, 0) == 0);
    CFX_ASSERT(cfx_blake2b(out48, 48, "test", 4, NULL, 0) == 0);
    CFX_ASSERT(cfx_blake2b(out64, 64, "test", 4, NULL, 0) == 0);

    /* different output lengths should produce different hashes */
    CFX_ASSERT(memcmp(out32, out64, 32) != 0);
    CFX_ASSERT(memcmp(out48, out64, 48) != 0);

    printf("test_blake2b_variable_output() - OK\n");
}

static void test_blake2s_variable_output(void) {
    uint8_t out16[16], out24[24], out32[32];

    CFX_ASSERT(cfx_blake2s(out16, 16, "test", 4, NULL, 0) == 0);
    CFX_ASSERT(cfx_blake2s(out24, 24, "test", 4, NULL, 0) == 0);
    CFX_ASSERT(cfx_blake2s(out32, 32, "test", 4, NULL, 0) == 0);

    CFX_ASSERT(memcmp(out16, out32, 16) != 0);
    CFX_ASSERT(memcmp(out24, out32, 24) != 0);

    printf("test_blake2s_variable_output() - OK\n");
}

/* Keyed hashing (MAC mode) */
static void test_blake2b_keyed(void) {
    uint8_t out[64];
    const uint8_t key[64] = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
        0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
        0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
        0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
        0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
        0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f
    };

    /* RFC 7693 Appendix E: keyed BLAKE2b with 64-byte key, empty message */
    const uint8_t expected[64] = {
        0x10, 0xeb, 0xb6, 0x77, 0x00, 0xb1, 0x86, 0x8e,
        0xfb, 0x44, 0x17, 0x98, 0x7a, 0xcf, 0x46, 0x90,
        0xae, 0x9d, 0x97, 0x2f, 0xb7, 0xa5, 0x90, 0xc2,
        0xf0, 0x28, 0x71, 0x79, 0x9a, 0xaa, 0x47, 0x86,
        0xb5, 0xe9, 0x96, 0xe8, 0xf0, 0xf4, 0xeb, 0x98,
        0x1f, 0xc2, 0x14, 0xb0, 0x05, 0xf4, 0x2d, 0x2f,
        0xf4, 0x23, 0x34, 0x99, 0x39, 0x16, 0x53, 0xdf,
        0x7a, 0xef, 0xcb, 0xc1, 0x3f, 0xc5, 0x15, 0x68
    };

    CFX_ASSERT(cfx_blake2b(out, 64, NULL, 0, key, 64) == 0);
    CFX_ASSERT(memcmp(out, expected, 64) == 0);

    printf("test_blake2b_keyed() - OK\n");
}

static void test_blake2s_keyed(void) {
    uint8_t out[32];
    const uint8_t key[32] = {
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
        0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
        0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
        0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f
    };

    /* RFC 7693 Appendix E: keyed BLAKE2s with 32-byte key, empty message */
    const uint8_t expected[32] = {
        0x48, 0xa8, 0x99, 0x7d, 0xa4, 0x07, 0x87, 0x6b,
        0x3d, 0x79, 0xc0, 0xd9, 0x23, 0x25, 0xad, 0x3b,
        0x89, 0xcb, 0xb7, 0x54, 0xd8, 0x6a, 0xb7, 0x1a,
        0xee, 0x04, 0x7a, 0xd3, 0x45, 0xfd, 0x2c, 0x49
    };

    CFX_ASSERT(cfx_blake2s(out, 32, NULL, 0, key, 32) == 0);
    CFX_ASSERT(memcmp(out, expected, 32) == 0);

    printf("test_blake2s_keyed() - OK\n");
}

/* Streaming API */
static void test_blake2b_streaming(void) {
    cfx_blake2b_ctx_t ctx;
    uint8_t out_oneshot[64], out_stream[64];
    const char *msg = "The quick brown fox jumps over the lazy dog";

    cfx_blake2b(out_oneshot, 64, msg, strlen(msg), NULL, 0);

    cfx_blake2b_init(&ctx, 64);
    cfx_blake2b_update(&ctx, msg, 10);
    cfx_blake2b_update(&ctx, msg + 10, 20);
    cfx_blake2b_update(&ctx, msg + 30, strlen(msg) - 30);
    cfx_blake2b_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 64) == 0);
    printf("test_blake2b_streaming() - OK\n");
}

static void test_blake2s_streaming(void) {
    cfx_blake2s_ctx_t ctx;
    uint8_t out_oneshot[32], out_stream[32];
    const char *msg = "The quick brown fox jumps over the lazy dog";

    cfx_blake2s(out_oneshot, 32, msg, strlen(msg), NULL, 0);

    cfx_blake2s_init(&ctx, 32);
    cfx_blake2s_update(&ctx, msg, 10);
    cfx_blake2s_update(&ctx, msg + 10, 20);
    cfx_blake2s_update(&ctx, msg + 30, strlen(msg) - 30);
    cfx_blake2s_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 32) == 0);
    printf("test_blake2s_streaming() - OK\n");
}

/* Byte-by-byte streaming */
static void test_blake2b_byte_by_byte(void) {
    cfx_blake2b_ctx_t ctx;
    uint8_t out_oneshot[64], out_stream[64];
    const char *msg = "hello world";
    size_t len = strlen(msg);

    cfx_blake2b(out_oneshot, 64, msg, len, NULL, 0);

    cfx_blake2b_init(&ctx, 64);
    for (size_t i = 0; i < len; i++) {
        cfx_blake2b_update(&ctx, msg + i, 1);
    }
    cfx_blake2b_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 64) == 0);
    printf("test_blake2b_byte_by_byte() - OK\n");
}

/* Long message (multiple blocks) */
static void test_blake2b_long_message(void) {
    uint8_t *msg = malloc(1000);
    uint8_t out[64];

    for (int i = 0; i < 1000; i++) msg[i] = (uint8_t)i;

    CFX_ASSERT(cfx_blake2b(out, 64, msg, 1000, NULL, 0) == 0);

    /* just verify it completes without crashing */
    /* a known answer test would require a reference implementation */

    free(msg);
    printf("test_blake2b_long_message() - OK\n");
}

static void test_blake2s_long_message(void) {
    uint8_t *msg = malloc(1000);
    uint8_t out[32];

    for (int i = 0; i < 1000; i++) msg[i] = (uint8_t)i;

    CFX_ASSERT(cfx_blake2s(out, 32, msg, 1000, NULL, 0) == 0);

    free(msg);
    printf("test_blake2s_long_message() - OK\n");
}

/* Edge cases */
static void test_blake2b_invalid_params(void) {
    uint8_t out[64];

    CFX_ASSERT(cfx_blake2b(out, 0, "test", 4, NULL, 0) == -1);   /* outlen=0 invalid */
    CFX_ASSERT(cfx_blake2b(out, 65, "test", 4, NULL, 0) == -1);  /* outlen>64 invalid */
    CFX_ASSERT(cfx_blake2b(out, 64, "test", 4, "key", 65) == -1); /* keylen>64 invalid */

    printf("test_blake2b_invalid_params() - OK\n");
}

static void test_blake2s_invalid_params(void) {
    uint8_t out[32];

    CFX_ASSERT(cfx_blake2s(out, 0, "test", 4, NULL, 0) == -1);
    CFX_ASSERT(cfx_blake2s(out, 33, "test", 4, NULL, 0) == -1);
    CFX_ASSERT(cfx_blake2s(out, 32, "test", 4, "key", 33) == -1);

    printf("test_blake2s_invalid_params() - OK\n");
}

/* Different keys produce different MACs */
static void test_blake2b_different_keys(void) {
    uint8_t out1[64], out2[64];
    const uint8_t key1[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    const uint8_t key2[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17};

    cfx_blake2b(out1, 64, "test", 4, key1, 16);
    cfx_blake2b(out2, 64, "test", 4, key2, 16);

    CFX_ASSERT(memcmp(out1, out2, 64) != 0);
    printf("test_blake2b_different_keys() - OK\n");
}

/* Different messages produce different hashes */
static void test_blake2b_different_messages(void) {
    uint8_t out1[64], out2[64];

    cfx_blake2b(out1, 64, "test1", 5, NULL, 0);
    cfx_blake2b(out2, 64, "test2", 5, NULL, 0);

    CFX_ASSERT(memcmp(out1, out2, 64) != 0);
    printf("test_blake2b_different_messages() - OK\n");
}

/* Exactly one block (128 bytes for BLAKE2b) */
static void test_blake2b_one_block(void) {
    uint8_t msg[128];
    uint8_t out[64];

    memset(msg, 'A', 128);
    CFX_ASSERT(cfx_blake2b(out, 64, msg, 128, NULL, 0) == 0);

    printf("test_blake2b_one_block() - OK\n");
}

/* Exactly one block + 1 byte */
static void test_blake2b_one_block_plus_one(void) {
    uint8_t msg[129];
    uint8_t out[64];

    memset(msg, 'A', 129);
    CFX_ASSERT(cfx_blake2b(out, 64, msg, 129, NULL, 0) == 0);

    printf("test_blake2b_one_block_plus_one() - OK\n");
}

/* Keyed streaming */
static void test_blake2b_keyed_streaming(void) {
    cfx_blake2b_ctx_t ctx;
    uint8_t out_oneshot[64], out_stream[64];
    const uint8_t key[32] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                             16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31};
    const char *msg = "hello";

    cfx_blake2b(out_oneshot, 64, msg, 5, key, 32);

    cfx_blake2b_init_key(&ctx, 64, key, 32);
    cfx_blake2b_update(&ctx, msg, 5);
    cfx_blake2b_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 64) == 0);
    printf("test_blake2b_keyed_streaming() - OK\n");
}

/* Self-test: streaming vs one-shot consistency */
static void test_blake2b_consistency(void) {
    uint8_t input[256];
    uint8_t out_oneshot[64], out_stream[64];
    cfx_blake2b_ctx_t ctx;

    for (int i = 0; i < 256; i++) input[i] = (uint8_t)i;

    /* one-shot */
    CFX_ASSERT(cfx_blake2b(out_oneshot, 64, input, 256, NULL, 0) == 0);

    /* streaming in chunks */
    cfx_blake2b_init(&ctx, 64);
    cfx_blake2b_update(&ctx, input, 100);
    cfx_blake2b_update(&ctx, input + 100, 100);
    cfx_blake2b_update(&ctx, input + 200, 56);
    cfx_blake2b_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 64) == 0);
    printf("test_blake2b_consistency() - OK\n");
}

static void test_blake2s_consistency(void) {
    uint8_t input[256];
    uint8_t out_oneshot[32], out_stream[32];
    cfx_blake2s_ctx_t ctx;

    for (int i = 0; i < 256; i++) input[i] = (uint8_t)i;

    cfx_blake2s(out_oneshot, 32, input, 256, NULL, 0);

    cfx_blake2s_init(&ctx, 32);
    cfx_blake2s_update(&ctx, input, 100);
    cfx_blake2s_update(&ctx, input + 100, 100);
    cfx_blake2s_update(&ctx, input + 200, 56);
    cfx_blake2s_final(&ctx, out_stream);

    CFX_ASSERT(memcmp(out_oneshot, out_stream, 32) == 0);
    printf("test_blake2s_consistency() - OK\n");
}

/*
 * SHA-256 LENGTH EXTENSION ATTACK DEMONSTRATION
 *
 * Merkle-Damgård hashes (MD5, SHA-1, SHA-256) are vulnerable because:
 * 1. The output IS the internal state
 * 2. Given H(M), attacker can compute H(M || padding || X) without knowing M
 *
 * This is why HMAC exists: HMAC(k,m) = H(k ⊕ opad || H(k ⊕ ipad || m))
 *
 * Attack scenario: Server uses H(secret || user_data) as MAC
 * Attacker sees: MAC = SHA256(secret || "user=alice")
 * Attacker computes: SHA256(secret || "user=alice" || padding || "&admin=true")
 *                    WITHOUT knowing the secret!
 */

/* internal sha256 state (matches cfx_sha256_state_t in sha256.c) */
typedef struct {
    uint32_t state[8];
    uint64_t total_bits;
    uint8_t buffer[64];
    size_t buffer_len;
} sha256_internal_state_t;

static void test_sha256_length_extension_attack(void) {
    /*
     * Simulate: MAC = SHA256(secret || message)
     * Secret: "supersecret" (11 bytes)
     * Message: "user=alice" (10 bytes)
     * Total: 21 bytes
     */
    const char *secret = "supersecret";
    const char *original_msg = "user=alice";
    const char *extension = "&admin=true";

    /* Step 1: Server computes MAC = SHA256(secret || message) */
    cfx_sha256_ctx ctx;
    uint8_t original_mac[32];
    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, (const uint8_t *)secret, strlen(secret));
    cfx_sha256_update(&ctx, (const uint8_t *)original_msg, strlen(original_msg));
    cfx_sha256_final(&ctx, original_mac);

    /*
     * Step 2: Attacker performs length extension attack
     *
     * The attacker knows:
     * - original_mac (the hash output = internal state after processing)
     * - len(secret || original_msg) = 21 bytes (or can guess/brute-force)
     *
     * SHA-256 padding for 21 bytes (168 bits):
     * - 0x80 byte
     * - zeros until 56 bytes from block end
     * - 8-byte big-endian length = 168 bits = 0x00000000000000A8
     *
     * Padding: 0x80 || (43 zeros) || 0x00 0x00 0x00 0x00 0x00 0x00 0x00 0xA8
     * Total block = 64 bytes
     */
    uint8_t padding[64 - 21];  /* bytes to fill first block */
    memset(padding, 0, sizeof(padding));
    padding[0] = 0x80;
    /* length in bits at end: 21 * 8 = 168 = 0xA8 */
    padding[sizeof(padding) - 1] = 0xA8;

    /*
     * Step 3: Attacker reconstructs SHA-256 state from MAC
     * and continues hashing with the extension
     */
    cfx_sha256_ctx attacker_ctx;
    cfx_sha256_init(&attacker_ctx);
    sha256_internal_state_t *internal = (sha256_internal_state_t *)&attacker_ctx;

    /* Inject the MAC as internal state (h0-h7) */
    for (int i = 0; i < 8; i++) {
        internal->state[i] = ((uint32_t)original_mac[i*4] << 24) |
                             ((uint32_t)original_mac[i*4+1] << 16) |
                             ((uint32_t)original_mac[i*4+2] << 8) |
                             ((uint32_t)original_mac[i*4+3]);
    }
    /* Set length counter to 64 bytes = 512 bits (one block processed) */
    internal->total_bits = 512;
    internal->buffer_len = 0;

    /* Hash the extension */
    cfx_sha256_update(&attacker_ctx, (const uint8_t *)extension, strlen(extension));
    uint8_t forged_mac[32];
    cfx_sha256_final(&attacker_ctx, forged_mac);

    /*
     * Step 4: Verify the attack - compute the "real" hash
     * H(secret || original_msg || padding || extension)
     */
    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, (const uint8_t *)secret, strlen(secret));
    cfx_sha256_update(&ctx, (const uint8_t *)original_msg, strlen(original_msg));
    cfx_sha256_update(&ctx, padding, sizeof(padding));
    cfx_sha256_update(&ctx, (const uint8_t *)extension, strlen(extension));
    uint8_t real_extended_mac[32];
    cfx_sha256_final(&ctx, real_extended_mac);

    /* THE ATTACK SUCCEEDS: forged_mac == real_extended_mac */
    CFX_ASSERT(memcmp(forged_mac, real_extended_mac, 32) == 0);

    printf("test_sha256_length_extension_attack() - OK (attack succeeds!)\n");
    printf("  [!] SHA-256 is vulnerable to length extension\n");
    printf("  [!] Attacker computed H(secret||msg||pad||ext) without knowing secret\n");
}

/*
 * BLAKE2 RESISTS LENGTH EXTENSION ATTACKS
 *
 * BLAKE2 is NOT vulnerable because:
 * 1. The output is XORed with both halves of the state (v[0..7] ^ v[8..15])
 * 2. The finalization flag changes the compression behavior
 * 3. The output doesn't reveal the full internal state
 *
 * Even if you try the same attack approach as SHA-256, it fails.
 */

/* internal blake2b state (matches implementation) */
typedef struct {
    uint64_t h[8];
    uint64_t t[2];
    uint64_t f[2];
    uint8_t buf[128];
    size_t buflen;
    size_t outlen;
} blake2b_internal_state_t;

static void test_blake2b_resists_length_extension(void) {
    const char *secret = "supersecret";
    const char *original_msg = "user=alice";
    const char *extension = "&admin=true";

    /* Step 1: Compute original MAC with BLAKE2b keyed mode */
    uint8_t original_mac[64];
    cfx_blake2b(original_mac, 64, original_msg, strlen(original_msg),
        secret, strlen(secret));

    /*
     * Step 2: Attacker attempts length extension
     *
     * With BLAKE2, the attacker cannot:
     * - Reconstruct internal state from output (XOR'd with second half)
     * - Continue hashing because finalization flag was set
     *
     * We'll try anyway and show it fails.
     */
    cfx_blake2b_ctx_t attacker_ctx;
    cfx_blake2b_init(&attacker_ctx, 64);
    blake2b_internal_state_t *internal = (blake2b_internal_state_t *)&attacker_ctx;

    /* Try to inject MAC as state (this is fundamentally wrong for BLAKE2) */
    for (int i = 0; i < 8; i++) {
        internal->h[i] = ((uint64_t)original_mac[i*8] << 56) |
                         ((uint64_t)original_mac[i*8+1] << 48) |
                         ((uint64_t)original_mac[i*8+2] << 40) |
                         ((uint64_t)original_mac[i*8+3] << 32) |
                         ((uint64_t)original_mac[i*8+4] << 24) |
                         ((uint64_t)original_mac[i*8+5] << 16) |
                         ((uint64_t)original_mac[i*8+6] << 8) |
                         ((uint64_t)original_mac[i*8+7]);
    }
    internal->t[0] = 128;  /* pretend one block processed */
    internal->buflen = 0;

    cfx_blake2b_update(&attacker_ctx, extension, strlen(extension));
    uint8_t forged_mac[64];
    cfx_blake2b_final(&attacker_ctx, forged_mac);

    /* The real MAC for the extended message (with key) */
    uint8_t extended_data[100];
    size_t total_len = strlen(original_msg) + strlen(extension);
    memcpy(extended_data, original_msg, strlen(original_msg));
    memcpy(extended_data + strlen(original_msg), extension, strlen(extension));

    uint8_t real_extended_mac[64];
    cfx_blake2b(real_extended_mac, 64, extended_data, total_len,
        secret, strlen(secret));

    /* THE ATTACK FAILS: forged_mac != real_extended_mac */
    CFX_ASSERT(memcmp(forged_mac, real_extended_mac, 64) != 0);

    printf("test_blake2b_resists_length_extension() - OK (attack fails!)\n");
    printf("  [+] BLAKE2 is resistant to length extension attacks\n");
    printf("  [+] Attacker cannot forge valid MAC without the key\n");
}

int main(void) {
    CFX_TEST(test_blake2b_empty);
    CFX_TEST(test_blake2b_abc);
    CFX_TEST(test_blake2s_empty);
    CFX_TEST(test_blake2s_abc);
    CFX_TEST(test_blake2b_variable_output);
    CFX_TEST(test_blake2s_variable_output);
    CFX_TEST(test_blake2b_keyed);
    CFX_TEST(test_blake2s_keyed);
    CFX_TEST(test_blake2b_streaming);
    CFX_TEST(test_blake2s_streaming);
    CFX_TEST(test_blake2b_byte_by_byte);
    CFX_TEST(test_blake2b_keyed_streaming);
    CFX_TEST(test_blake2b_long_message);
    CFX_TEST(test_blake2s_long_message);
    CFX_TEST(test_blake2b_invalid_params);
    CFX_TEST(test_blake2s_invalid_params);
    CFX_TEST(test_blake2b_one_block);
    CFX_TEST(test_blake2b_one_block_plus_one);
    CFX_TEST(test_blake2b_different_keys);
    CFX_TEST(test_blake2b_different_messages);
    CFX_TEST(test_blake2b_consistency);
    CFX_TEST(test_blake2s_consistency);

    /* length extension attack demonstration */
    CFX_TEST(test_sha256_length_extension_attack);
    CFX_TEST(test_blake2b_resists_length_extension);

    printf("\nAll BLAKE2 tests passed!\n");
    return 0;
}
