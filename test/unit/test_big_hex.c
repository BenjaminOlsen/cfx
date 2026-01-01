/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_hex.c - Tests for big integer hex conversion */

#include "test_common.h"

static void test_hex_zero_empty_n(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    size_t len = 12345;
    char* s = cfx_big_hex_alloc(&b, &len);
    CFX_ASSERT(s);
    CFX_ASSERT(strcmp(s, "0") == 0);
    CFX_ASSERT(len == 1);
    CFX_ASSERT(s[len] == '\0');
    free(s);
    cfx_big_free(&b);
    PRINT_TEST(1);
}

static void test_hex_zero_explicit_limb_zero(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs[] = {0};
    cfx_big_from_limbs(&b, limbs, 1);
    size_t len = 0;
    char* s = cfx_big_hex_alloc(&b, &len);
    CFX_ASSERT(s);
    CFX_ASSERT(strcmp(s, "0") == 0);
    CFX_ASSERT(len == 1);
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_single_limb_basic(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs1[] = {UINT64_C(0x1)};
    cfx_big_from_limbs(&b, limbs1, 1);
    size_t len1 = 0;
    char* s1 = cfx_big_hex_alloc(&b, &len1);
    CFX_ASSERT(strcmp(s1, "1") == 0);
    CFX_ASSERT(len1 == strlen("1"));
    CFX_ASSERT(s1[len1] == '\0');
    free(s1);

    b.limb[0] = UINT64_C(0xabcdef);
    size_t len2 = 0;
    char* s2 = cfx_big_hex_alloc(&b, &len2);
    CFX_ASSERT(strcmp(s2, "abcdef") == 0);
    CFX_ASSERT(len2 == strlen("abcdef"));
    CFX_ASSERT(s2[len2] == '\0');
    cfx_big_free(&b);
    free(s2);
    PRINT_TEST(1);
}

static void test_hex_single_limb_hex_digit_count(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    uint64_t v = UINT64_C(0x1000000000000000);
    cfx_big_from_u64(&b, v);
    size_t len = 0;
    char* s = cfx_big_hex_alloc(&b, &len);
    CFX_ASSERT(strcmp(s, "1000000000000000") == 0);
    CFX_ASSERT(len == 16);
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_two_limbs_padding(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs[2];
    limbs[0] = (cfx_limb_t)1;
    limbs[1] = (cfx_limb_t)1;
    cfx_big_from_limbs(&b, limbs, 2);
    size_t len = 0;
    char* s = cfx_big_hex_alloc(&b, &len);
#if (CFX_LIMB_BITS == 64)
    CFX_ASSERT(strcmp(s, "10000000000000001") == 0);
    CFX_ASSERT(len == 17);
#elif (CFX_LIMB_BITS == 32)
    CFX_ASSERT(strcmp(s, "100000001") == 0);
    CFX_ASSERT(len == 9);
#endif
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_two_limbs_mixed_digits(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs[2];
    char* s = NULL;
    size_t szout = 0;
#if (CFX_LIMB_BITS==64)
    limbs[0] = UINT64_C(0x0011223344556677);
    limbs[1] = UINT64_C(0x0000000000000ABC);
    cfx_big_from_limbs(&b, limbs, 2);
    s = cfx_big_hex_alloc(&b, &szout);
    CFX_ASSERT(strcmp(s, "abc0011223344556677") == 0);
#else
    limbs[0] = UINT32_C(0x00112233);
    limbs[1] = UINT32_C(0x00000ABC);
    cfx_big_from_limbs(&b, limbs, 2);
    s = cfx_big_hex_alloc(&b, &szout);
    CFX_ASSERT(strcmp(s, "abc00112233") == 0);
#endif
    CFX_ASSERT(s != NULL);
    CFX_ASSERT(strlen(s) == szout);
    CFX_ASSERT(s[szout] == '\0');
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_leading_zero_limb_skipped(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs3[3];
#if CFX_LIMB_BITS == 64
    limbs3[0] = UINT64_C(0xDEADBEEFCAFEBABE);
    limbs3[1] = UINT64_C(0x0000000000000123);
    limbs3[2] = UINT64_C(0x0000000000000000);
    const char* expect_low = "deadbeefcafebabe";
    char expect_buf[3 + 16 + 1];
#elif CFX_LIMB_BITS == 32
    limbs3[0] = UINT64_C(0xDEADBEEF);
    limbs3[1] = UINT64_C(0x00000123);
    limbs3[2] = UINT64_C(0x00000000);
    const char* expect_low = "deadbeef";
    char expect_buf[3 + 8 + 1];
#endif
    cfx_big_from_limbs(&b, limbs3, 3);
    size_t len = 0;
    char* s = cfx_big_hex_alloc(&b, &len);
    snprintf(expect_buf, sizeof(expect_buf), "%s%s", "123", expect_low);
    CFX_ASSERT(strcmp(s, expect_buf) == 0);
    CFX_ASSERT(len == strlen(expect_buf));
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}

static void test_hex_no_leading_zeros_on_msl(void) {
    cfx_big_t b;
    cfx_big_init(&b);
    cfx_limb_t limbs[] = {
        (cfx_limb_t)0,
        (cfx_limb_t)0xAB
    };
    cfx_big_from_limbs(&b, limbs, 2);
    size_t len = 0;
    char* s = cfx_big_hex_alloc(&b, &len);
#if CFX_LIMB_BITS == 64
    CFX_ASSERT(strcmp(s, "ab0000000000000000") == 0);
    CFX_ASSERT(len == 2 + 16);
#else
    CFX_ASSERT(strcmp(s, "ab00000000") == 0);
    CFX_ASSERT(len == 2 + 8);
#endif
    CFX_ASSERT(s[len] == '\0');
    cfx_big_free(&b);
    free(s);
    PRINT_TEST(1);
}

int main(void) {
    CFX_TEST(test_hex_zero_empty_n);
    CFX_TEST(test_hex_zero_explicit_limb_zero);
    CFX_TEST(test_hex_single_limb_basic);
    CFX_TEST(test_hex_single_limb_hex_digit_count);
    CFX_TEST(test_hex_two_limbs_padding);
    CFX_TEST(test_hex_two_limbs_mixed_digits);
    CFX_TEST(test_hex_leading_zero_limb_skipped);
    CFX_TEST(test_hex_no_leading_zeros_on_msl);
    puts("OK");
    return 0;
}
