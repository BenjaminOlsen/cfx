/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_shift.c - Tests for big integer bit shift operations */

#include "test_common.h"

static void test_shifts(void) {
    cfx_big_t a, b;
    cfx_big_init(&a);
    cfx_big_init(&b);

    cfx_limb_t lsb = 123, msb = 283912;
    cfx_limb_t limbs[] = {lsb, msb};
    cfx_big_from_limbs(&a, limbs, sizeof(limbs)/sizeof(limbs[0]));
    cfx_big_copy(&b, &a);

    cfx_big_shl_bits(&b, &b, 1);
    CFX_ASSERT(b.limb[1] == 2 * msb);
    CFX_ASSERT(b.limb[0] == 2 * lsb);

    cfx_big_shl_bits(&b, &b, 1);
    CFX_ASSERT(b.limb[1] == 4 * msb);
    CFX_ASSERT(b.limb[0] == 4 * lsb);

    cfx_limb_t l0 = b.limb[0];
    cfx_limb_t l1 = b.limb[1];

    cfx_big_shl_bits(&b, &b, CFX_LIMB_BITS+1);

    CFX_ASSERT(b.limb[2] == l1*2);
    CFX_ASSERT(b.limb[1] == l0*2);
    CFX_ASSERT(b.limb[0] == 0);

    l0 = b.limb[0];
    l1 = b.limb[1];

    cfx_big_shr_bits(&b, &b, 1);
    CFX_ASSERT(b.limb[1] == l1/2);
    CFX_ASSERT(b.limb[0] == l0/2);

    cfx_big_free(&a);
    cfx_big_free(&b);
}

static void test_shl_basic_identity(void) {
    CHECK_SHL_CASE("0", 0, "0");
    CHECK_SHL_CASE("1", 0, "1");
    CHECK_SHL_CASE("deadbeef", 0, "deadbeef");
}

static void test_shl_create_top_limb(void) {
#if CFX_LIMB_BITS == 64
    CHECK_SHL_CASE("8000000000000000", 1, "10000000000000000");
    CHECK_SHL_CASE("1", 64, "10000000000000000");
    CHECK_SHL_CASE("1", 68, "100000000000000000");
#else
    CHECK_SHL_CASE("80000000", 1, "100000000");
    CHECK_SHL_CASE("1", 32, "100000000");
    CHECK_SHL_CASE("1", 36, "1000000000");
#endif
}

static void test_shl_cross_limb_1bit(void) {
#if CFX_LIMB_BITS == 64
    CHECK_SHL_CASE("100000000000000008000000000000001", 1, "200000000000000010000000000000002");
#else
    CHECK_SHL_CASE("180000001", 1, "300000002");
#endif
}

static void test_shr_basic_identity(void) {
    CHECK_SHR_CASE("0", 0, "0");
    CHECK_SHR_CASE("1", 0, "1");
    CHECK_SHR_CASE("deadbeef", 0, "deadbeef");
}

static void test_shr_drop_whole_limb(void) {
#if CFX_LIMB_BITS == 64
    CHECK_SHR_CASE("10000000000000000", 64, "1");
    CHECK_SHR_CASE("100000000000000000", 68, "1");
#else
    CHECK_SHR_CASE("100000000", 32, "1");
    CHECK_SHR_CASE("1000000000", 36, "1");
#endif
}

static void test_shr_to_zero(void) {
#if CFX_LIMB_BITS == 64
    CHECK_SHR_CASE("1", 65, "0");
    CHECK_SHR_CASE("ffffffffffffffff", 128, "0");
#else
    CHECK_SHR_CASE("1", 33, "0");
    CHECK_SHR_CASE("ffffffff", 64, "0");
#endif
}

static void test_shr_cross_limb_carry_6bits(void) {
#if CFX_LIMB_BITS == 64
    CHECK_SHR_CASE("3e0000000000000c40", 6, "f800000000000031");
#else
    CHECK_SHR_CASE("3e00000c40", 6, "f8000031");
#endif
}

static void test_shr_mixed_cases(void) {
#if CFX_LIMB_BITS == 64
    CHECK_SHR_CASE("100000000000000000000", 4, "10000000000000000000");
    CHECK_SHR_CASE("abcdef0123456789", 4, "abcdef012345678");
    CHECK_SHR_CASE("abcdef0123456789", 60, "a");
#else
    CHECK_SHR_CASE("10000000000", 4, "1000000000");
    CHECK_SHR_CASE("abcdef01", 4, "abcdef0");
    CHECK_SHR_CASE("abcdef01", 28, "a");
#endif
}

static void test_shl_zero_shift(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 42);
    cfx_big_shl_bits_eq(&a, 0);
    CFX_ASSERT(cfx_big_eq_u64(&a, 42));
    cfx_big_free(&a);
}

static void test_shr_zero_shift(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 42);
    cfx_big_shr_bits_eq(&a, 0);
    CFX_ASSERT(cfx_big_eq_u64(&a, 42));
    cfx_big_free(&a);
}

static void test_shl_zero_value(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_shl_bits_eq(&a, 10);
    CFX_ASSERT(cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

static void test_shr_zero_value(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_shr_bits_eq(&a, 10);
    CFX_ASSERT(cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

static void test_shr_more_than_size(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 0xFF);
    cfx_big_shr_bits_eq(&a, 1000);
    CFX_ASSERT(cfx_big_is_zero(&a));
    cfx_big_free(&a);
}

static void test_shl_pure_limb_shift(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 1);
    cfx_big_shl_bits_eq(&a, CFX_LIMB_BITS);  /* pure limb shift, bit_shift==0 */
    CFX_ASSERT(a.n == 2);
    CFX_ASSERT(a.limb[0] == 0);
    CFX_ASSERT(a.limb[1] == 1);
    cfx_big_free(&a);
}

static void test_shr_pure_limb_shift(void) {
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_limb(&a, 1);
    cfx_big_shl_bits_eq(&a, CFX_LIMB_BITS);
    CFX_ASSERT(a.n == 2);
    cfx_big_shr_bits_eq(&a, CFX_LIMB_BITS);  /* pure limb shift back */
    CFX_ASSERT(cfx_big_eq_u64(&a, 1));
    cfx_big_free(&a);
}

int main(void) {
    CFX_TEST(test_shifts);
    CFX_TEST(test_shl_basic_identity);
    CFX_TEST(test_shl_create_top_limb);
    CFX_TEST(test_shl_cross_limb_1bit);
    CFX_TEST(test_shr_basic_identity);
    CFX_TEST(test_shr_drop_whole_limb);
    CFX_TEST(test_shr_to_zero);
    CFX_TEST(test_shr_cross_limb_carry_6bits);
    CFX_TEST(test_shr_mixed_cases);
    CFX_TEST(test_shl_zero_shift);
    CFX_TEST(test_shr_zero_shift);
    CFX_TEST(test_shl_zero_value);
    CFX_TEST(test_shr_zero_value);
    CFX_TEST(test_shr_more_than_size);
    CFX_TEST(test_shl_pure_limb_shift);
    CFX_TEST(test_shr_pure_limb_shift);
    puts("OK");
    return 0;
}
