/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx/sha512.h"
#include "cfx/macros.h"

static int hex_nibble(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

static int hex_to_bytes(const char *hex, uint8_t *out, size_t out_len) {
    size_t hex_len = strlen(hex);
    if (hex_len % 2 != 0) return -1;
    if (out_len * 2 != hex_len) return -1;

    for (size_t i = 0; i < out_len; ++i) {
        int hi = hex_nibble(hex[2 * i]);
        int lo = hex_nibble(hex[2 * i + 1]);
        if (hi < 0 || lo < 0) return -1;
        out[i] = (uint8_t)((hi << 4) | lo);
    }
    return 0;
}

static void print_hex(const uint8_t *buf, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        printf("%02x", buf[i]);
    }
}

static void cfx_sha512_once(const uint8_t *msg, size_t len, uint8_t out[64]) {
    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, msg, len);
    cfx_sha512_final(&ctx, out);
}

static void check_vector(const char *name, const uint8_t *msg, size_t msg_len, const char *expected_hex) {
    uint8_t got[64];
    uint8_t expect[64];

    if (hex_to_bytes(expected_hex, expect, sizeof(expect)) != 0) {
        fprintf(stderr, "[%s] invalid expected hex string\n", name);
        abort();
    }

    cfx_sha512_once(msg, msg_len, got);

    if (memcmp(got, expect, sizeof(got)) != 0) {
        fprintf(stderr, "[%s] SHA-512 mismatch!\n", name);
        fprintf(stderr, "  expected: ");
        print_hex(expect, sizeof(expect));
        fprintf(stderr, "\n");
        fprintf(stderr, "       got: ");
        print_hex(got, sizeof(got));
        fprintf(stderr, "\n");
        abort();
    } else {
        printf("[OK] %s\n", name);
    }
}

static void check_incremental(const char *name, const uint8_t *msg, size_t msg_len) {
    uint8_t dig_one[64];
    uint8_t dig_inc[64];

    cfx_sha512_once(msg, msg_len, dig_one);

    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);

    size_t remaining = msg_len;
    const uint8_t *p = msg;
    size_t chunk_sizes[] = {1, 2, 3, 7, 128, 13, 5, 8, 64};
    size_t num_chunks = sizeof(chunk_sizes) / sizeof(chunk_sizes[0]);
    size_t idx = 0;

    while (remaining > 0) {
        size_t chunk = chunk_sizes[idx++ % num_chunks];
        if (chunk > remaining) chunk = remaining;
        cfx_sha512_update(&ctx, p, chunk);
        p += chunk;
        remaining -= chunk;
    }

    cfx_sha512_final(&ctx, dig_inc);

    if (memcmp(dig_one, dig_inc, sizeof(dig_one)) != 0) {
        fprintf(stderr, "[%s] incremental vs one-shot mismatch!\n", name);
        fprintf(stderr, "  one-shot: ");
        print_hex(dig_one, sizeof(dig_one));
        fprintf(stderr, "\n  incremental: ");
        print_hex(dig_inc, sizeof(dig_inc));
        fprintf(stderr, "\n");
        abort();
    } else {
        printf("[OK] %s (incremental)\n", name);
    }
}

static void test_nist_vectors(void) {
    /* FIPS 180-4 test vectors for SHA-512 */

    /* empty string */
    check_vector("SHA512(\"\")", (const uint8_t *)"", 0,
        "cf83e1357eefb8bdf1542850d66d8007d620e4050b5715dc83f4a921d36ce9ce"
        "47d0d13c5d85f2b0ff8318d2877eec2f63b931bd47417a81a538327af927da3e");

    /* "abc" */
    check_vector("SHA512(\"abc\")", (const uint8_t *)"abc", 3,
        "ddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a"
        "2192992a274fc1a836ba3c23a3feebbd454d4423643ce80e2a9ac94fa54ca49f");

    /* two-block message (896 bits) */
    const char *msg3 =
        "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn"
        "hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu";
    check_vector("SHA512(112-byte)", (const uint8_t *)msg3, strlen(msg3),
        "8e959b75dae313da8cf4f72814fc143f8f7779c6eb9f7fa17299aeadb6889018"
        "501d289e4900f7e4331b99dec4b5433ac7d329eeb6dd26545e96e55b874be909");

    /* 1,000,000 x 'a' */
    size_t N = 1000000;
    uint8_t *million = malloc(N);
    memset(million, 'a', N);
    check_vector("SHA512(\"a\" x 1000000)", million, N,
        "e718483d0ce769644e2e42c7bc15b4638e1f98b13b2044285632a803afa973eb"
        "de0ff244877ea60a4cb0432ce577c31beb009c5c2c49aa2e4eadb217ad8cc09b");
    free(million);
}

static void test_incremental_hashing(void) {
    const uint8_t *msg = (const uint8_t *)"abc";
    check_incremental("SHA512(\"abc\")", msg, 3);

    /* Longer message */
    const char *msg2 =
        "abcdefghbcdefghicdefghijdefghijkefghijklfghijklmghijklmn"
        "hijklmnoijklmnopjklmnopqklmnopqrlmnopqrsmnopqrstnopqrstu";
    check_incremental("SHA512(112-byte)", (const uint8_t *)msg2, strlen(msg2));

    /* Million a's */
    size_t N = 1000000;
    uint8_t *million = malloc(N);
    memset(million, 'a', N);
    check_incremental("SHA512(\"a\" x 1000000)", million, N);
    free(million);
}

static void test_convenience_function(void) {
    /* Test cfx_sha512() one-shot convenience function */
    const uint8_t msg[] = "abc";
    uint8_t got[64];
    uint8_t expect[64];

    const char *expected_hex =
        "ddaf35a193617abacc417349ae20413112e6fa4e89a97ea20a9eeee64b55d39a"
        "2192992a274fc1a836ba3c23a3feebbd454d4423643ce80e2a9ac94fa54ca49f";

    hex_to_bytes(expected_hex, expect, 64);
    cfx_sha512(got, msg, 3);

    if (memcmp(got, expect, 64) != 0) {
        fprintf(stderr, "cfx_sha512() one-shot function failed!\n");
        abort();
    }
    printf("[OK] cfx_sha512() one-shot\n");
}

int main(void) {
    CFX_TEST(test_nist_vectors);
    CFX_TEST(test_incremental_hashing);
    CFX_TEST(test_convenience_function);
    printf("All SHA-512 tests passed.\n");
    return 0;
}
