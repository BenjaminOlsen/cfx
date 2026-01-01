/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_bytes.c - Tests for big integer byte conversion */

#include "test_common.h"

static void test_big_to_bytes_be_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    cfx_big_from_limb(&b, 0);

    /* size query */
    size_t need = 123;
    int rc = cfx_big_to_bytes_be(NULL, &need, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(need == 0);

    /* write into buffer */
    uint8_t out[8];
    memset(out, 0xAA, sizeof(out));
    size_t out_len = sizeof(out);
    rc = cfx_big_to_bytes_be(out, &out_len, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(out_len == 0);

    cfx_big_free(&b);
}

static void test_big_to_bytes_be(void) {

    struct {
        const char*hex_in;
        const char *hex_be;
    } cases[] = {
        { "01", "01" },
        { "ff", "FF" },
        { "0100", "0100" },
        { "123456", "123456" },
        { "1234567890abcdef", "1234567890ABCDEF" },
        { "0102030405060708090a0b0c0d0e0f10", "0102030405060708090A0B0C0D0E0F10" },
        { "8000000000000000", "8000000000000000" },         /* top bit set */
        { "ffffffffffffffffffffffffffffffff", "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF" },
        { "010000000000000000", "010000000000000000" },
        { "010203040506070809090909090909", "010203040506070809090909090909"}
    };

    cfx_big_t b;
    cfx_big_init(&b);

    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); i++) {
        int ok = cfx_big_from_hex(&b, cases[i].hex_in);
        CFX_ASSERT(ok == 0);

        size_t need = 0;
        int rc = cfx_big_to_bytes_be(NULL, &need, &b);
        CFX_ASSERT(rc == 0);

        if (strcmp(cases[i].hex_in, "0") == 0) {
            CFX_ASSERT(need == 0);
            continue;
        }

        size_t exp_len = strlen(cases[i].hex_be)/2;
        CFX_ASSERT(need == exp_len);

        uint8_t exp[64];
        hex_to_bytes(exp, exp_len, cases[i].hex_be);

        uint8_t out[64];
        size_t out_len = sizeof(out);
        rc = cfx_big_to_bytes_be(out, &out_len, &b);
        CFX_ASSERT(rc == 0);
        CFX_ASSERT(out_len == exp_len);
        CFX_ASSERT(memcmp(out, exp, exp_len) == 0);
    }

    cfx_big_free(&b);
}

static void test_big_to_bytes_be_buffer_too_small(void) {
    cfx_big_t b;
    cfx_big_init(&b);

    cfx_big_from_limb(&b, 0x010203u);

    size_t need = 0;
    int rc = cfx_big_to_bytes_be(NULL, &need, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(need == 3);

    uint8_t out[3];
    size_t out_len = 2; /* too small on purpose */
    rc = cfx_big_to_bytes_be(out, &out_len, &b);
    CFX_ASSERT(rc == -1);

    cfx_big_free(&b);
}

static void test_to_bin_zero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    size_t sz = 0;
    char* s = cfx_big_bin_alloc(&a, &sz);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strcmp(s, "0") == 0);
    CFX_ASSERT(sz == 1);
    free(s);
    cfx_big_free(&a);
}

static void test_to_bin_nonzero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 5);  /* 101 in binary */
    size_t sz = 0;
    char* s = cfx_big_bin_alloc(&a, &sz);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strcmp(s, "101") == 0);
    free(s);
    cfx_big_free(&a);
}

static void test_to_b64_zero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    size_t sz = 0;
    char* s = cfx_big_b64_alloc(&a, &sz);
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strcmp(s, "0") == 0);
    free(s);
    cfx_big_free(&a);
}

static void test_from_bytes_be_all_zeros(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    uint8_t zeros[] = {0, 0, 0, 0};
    int rc = cfx_big_from_bytes_be(&a, zeros, 4);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

static void test_from_bytes_be_null(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    int rc = cfx_big_from_bytes_be(&a, NULL, 0);
    CFX_ASSERT(rc == 0);
    rc = cfx_big_from_bytes_be(&a, NULL, 5);  /* NULL with nonzero len */
    CFX_ASSERT(rc != 0);
    cfx_big_free(&a);
}

static void test_to_bytes_be_query_zero(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    size_t len = 100;
    int rc = cfx_big_to_bytes_be(NULL, &len, &a);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(len == 0);
    cfx_big_free(&a);
}

int main(void) {
    CFX_TEST(test_big_to_bytes_be_zero);
    CFX_TEST(test_big_to_bytes_be);
    CFX_TEST(test_big_to_bytes_be_buffer_too_small);
    CFX_TEST(test_to_bin_zero);
    CFX_TEST(test_to_bin_nonzero);
    CFX_TEST(test_to_b64_zero);
    CFX_TEST(test_from_bytes_be_all_zeros);
    CFX_TEST(test_from_bytes_be_null);
    CFX_TEST(test_to_bytes_be_query_zero);
    puts("OK");
    return 0;
}
