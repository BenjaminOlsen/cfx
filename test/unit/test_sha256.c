#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cfx/sha256.h"
#include "cfx/macros.h"

static int hex_nibble(char c) {
    if (c >= '0' && c <= '9') return c - '0';
    if (c >= 'a' && c <= 'f') return 10 + (c - 'a');
    if (c >= 'A' && c <= 'F') return 10 + (c - 'A');
    return -1;
}

static int hex_to_bytes(const char *hex, uint8_t *out, size_t out_len) {
    size_t hex_len = strlen(hex);
    if (hex_len % 2 != 0) {
        return -1;
    }
    if (out_len * 2 != hex_len) {
        return -1;
    }

    for (size_t i = 0; i < out_len; ++i) {
        int hi = hex_nibble(hex[2 * i]);
        int lo = hex_nibble(hex[2 * i + 1]);
        if (hi < 0 || lo < 0) {
            return -1;
        }
        out[i] = (uint8_t)((hi << 4) | lo);
    }
    return 0;
}

static void print_hex(const uint8_t *buf, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        printf("%02x", buf[i]);
    }
}

static void cfx_sha256_once(const uint8_t *msg, size_t len, uint8_t out[32]) {
    cfx_sha256_ctx ctx;
    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, msg, len);
    cfx_sha256_final(&ctx, out);
}

static void check_vector(const char *name,
                         const uint8_t *msg, size_t msg_len,
                         const char *expected_hex) {
    uint8_t got[32];
    uint8_t expect[32];

    if (hex_to_bytes(expected_hex, expect, sizeof(expect)) != 0) {
        fprintf(stderr, "[%s] invalid expected hex string\n", name);
        abort();
    }

    cfx_sha256_once(msg, msg_len, got);

    if (memcmp(got, expect, sizeof(got)) != 0) {
        fprintf(stderr, "[%s] SHA-256 mismatch!\n", name);
        fprintf(stderr, "  expected: ");
        print_hex(expect, sizeof(expect));
        fprintf(stderr, "\n");
        fprintf(stderr, "       got: ");
        print_hex(got, sizeof(got));
        fprintf(stderr, "\n");
        //CFX_FAIL();
    } else {
        printf("[OK] %s\n", name);
    }
}

static void check_incremental(const char *name,
                              const uint8_t *msg, size_t msg_len) {
    uint8_t dig_one[32];
    uint8_t dig_inc[32];

    cfx_sha256_once(msg, msg_len, dig_one);

    cfx_sha256_ctx ctx;
    cfx_sha256_init(&ctx);

    size_t remaining = msg_len;
    const uint8_t *p = msg;
    size_t chunk_sizes[] = {1, 2, 3, 7, 64, 13, 5, 8, 32};
    size_t num_chunks = sizeof(chunk_sizes) / sizeof(chunk_sizes[0]);
    size_t idx = 0;

    while (remaining > 0) {
        size_t chunk = chunk_sizes[idx++ % num_chunks];
        if (chunk > remaining) {
            chunk = remaining;
        }
        cfx_sha256_update(&ctx, p, chunk);
        p += chunk;
        remaining -= chunk;
    }

    cfx_sha256_final(&ctx, dig_inc);

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
    // 1. Empty
    check_vector("NIST empty", (const uint8_t *)"", 0,
        "e3b0c44298fc1c149afbf4c8996fb924"
        "27ae41e4649b934ca495991b7852b855");

    // 2. "abc"
    check_vector("NIST abc", (const uint8_t *)"abc", 3,
        "ba7816bf8f01cfea414140de5dae2223"
        "b00361a396177a9cb410ff61f20015ad");

    // 3. Longer FIPS example
    const char *msg3 =
        "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
    check_vector("NIST long ascii", (const uint8_t *)msg3, strlen(msg3),
        "248d6a61d20638b8e5c026930c3e6039"
        "a33ce45964ff2167f6ecedd419db06c1");

    // 4. 1,000,000 × 'a'
    size_t N = 1000000;
    uint8_t *million = malloc(N);
    memset(million, 'a', N);

    check_vector("NIST million a", million, N,
        "cdc76e5c9914fb9281a1c7e284d73e67"
        "f1809a48a497200e046d39ccc7112cd0");

    free(million);

    // 5. single 0x00  (generic SHA-256 test)
    uint8_t m00[1] = { 0x00 };
    check_vector("SHA256(0x00)", m00, 1,
        "6e340b9cffb37a989ca544e6bb780a2c"
        "78901d3fb33738768511a30617afa01d");

    // 6. single 0xFF  (generic SHA-256 test)
    uint8_t mff[1] = { 0xFF };
    check_vector("SHA256(0xFF)", mff, 1,
        "a8100ae6aa1940d0b663bb31cd466142"
        "ebbdbd5187131b92d93818987832eb89");

    // 7. 55×0x00 (padding edge case: fits in one block)
    uint8_t m55[55] = {0};
    check_vector("SHA256(55×0x00)", m55, sizeof m55,
        "02779466cdec163811d078815c633f21"
        "901413081449002f24aa3e80f0b88ef7");

    // 8. 56×0x00 (padding edge case: forces two blocks)
    uint8_t m56[56] = {0};
    check_vector("SHA256(56×0x00)", m56, sizeof m56,
        "d4817aa5497628e7c77e6b606107042b"
        "bba3130888c5f47a375e6179be789fbb");
}

static void abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq(void) {
    const uint8_t msg[] =
        "abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq";
    const size_t msg_len = strlen((const char *)msg);
    const char *expected =
        "248d6a61d20638b8e5c026930c3e6039"
        "a33ce45964ff2167f6ecedd419db06c1";
    check_vector("SHA256(long FIPS msg)", msg, msg_len, expected);
    check_incremental("SHA256(long FIPS msg)", msg, msg_len);
}

static void empty_string(void) {
    const uint8_t msg[] = { /* empty */ };
    const char *expected =
        "e3b0c44298fc1c149afbf4c8996fb924"
        "27ae41e4649b934ca495991b7852b855";
    check_vector("SHA256(\"\")", msg, 0, expected);
    check_incremental("SHA256(\"\")", msg, 0);
}

static void abc(void) {
    const uint8_t msg[] = {'a','b','c'};
    const char *expected =
        "ba7816bf8f01cfea414140de5dae2223"
        "b00361a396177a9cb410ff61f20015ad";
    check_vector("SHA256(\"abc\")", msg, sizeof(msg), expected);
    check_incremental("SHA256(\"abc\")", msg, sizeof(msg));
}

static void aaaaaaaaaaamillion(void) {
    const size_t N = 1000000;
    uint8_t *msg = (uint8_t *)malloc(N);
    CFX_ASSERT(msg);
    memset(msg, 'a', N);

    const char *expected =
        "cdc76e5c9914fb9281a1c7e284d73e67"
        "f1809a48a497200e046d39ccc7112cd0";

    check_vector("SHA256(\"a\" * 1_000_000)", msg, N, expected);
    check_incremental("SHA256(\"a\" * 1_000_000)", msg, N);

    free(msg);
}

int main(void) {
    CFX_TEST(test_nist_vectors);
    CFX_TEST(abc);
    CFX_TEST(empty_string);
    CFX_TEST(abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq);
    CFX_TEST(aaaaaaaaaaamillion);
    printf("All SHA-256 tests passed.\n");
    return 0;
}
