/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_div.c - Tests for big integer division */

#include "test_common.h"

static void test_big_div_divide_by_zero(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_limb(&n, 123);
    cfx_big_from_limb(&d, 0);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == -1);

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_zero_dividend(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_limb(&n, 0);
    cfx_big_from_limb(&d, 42);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "0");
    expect_dec_eq(&r, "0");

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_n_less_than_d(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_limb(&n, 123456);
    cfx_big_from_limb(&d, 123456789);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "0");
    CFX_ASSERT(cfx_big_eq(&r, &n));

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_equal_numbers(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_dec(&n, "1234567890123456789012345678901234567890");
    cfx_big_copy(&d, &n);

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    expect_dec_eq(&q, "1");
    expect_dec_eq(&r, "0");

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_single_limb_divisor_property(void) {
    cfx_big_t n, d, q, r;
    cfx_big_init(&n);
    cfx_big_init(&d);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_from_dec(&n, "340282366920938463463374607431768211455");
    cfx_big_from_limb(&d, UINT64_C(123456789));

    int rc = cfx_big_divrem(&q, &r, &n, &d);
    CFX_ASSERT(rc == 0);
    assert_n_eq_qd_plus_r(&n, &q, &d, &r);

    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&q);
    cfx_big_free(&r);
}

static void test_big_div_multi_limb_divisor_exact_and_remainder(void) {
    cfx_big_t a, b, r, n, q, rem;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_init(&rem);
    cfx_big_from_dec(&a, "123456789012345678901234567890123456789");
    cfx_big_from_dec(&b, "987654321098765432109876543210987654321");
    cfx_big_from_dec(&r, "12345678901234567890");

    cfx_big_copy(&n, &a);
    cfx_big_mul_eq(&n, &b);
    cfx_big_add_eq(&n, &r);

    int rc = cfx_big_divrem(&q, &rem, &n, &b);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(cfx_big_eq(&q, &a));
    CFX_ASSERT(cfx_big_eq(&rem, &r));
    assert_n_eq_qd_plus_r(&n, &q, &b, &rem);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&r);
    cfx_big_free(&n);
    cfx_big_free(&q);
    cfx_big_free(&rem);
}

static void test_big_div_in_place_eq_with_remainder(void) {
    cfx_big_t a, b, n, rem, forty_two;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&rem);
    cfx_big_init(&forty_two);
    cfx_big_from_dec(&a, "1122334455667788990011223344556677889900");
    cfx_big_from_dec(&b, "18446744073709551616");
    cfx_big_copy(&n, &a);
    cfx_big_mul_eq(&n, &b);
    cfx_big_from_limb(&forty_two, 42);
    cfx_big_add_eq(&n, &forty_two);

    int rc = cfx_big_divrem_eq(&n, &b, &rem);
    CFX_ASSERT(rc == 0);
    CFX_ASSERT(cfx_big_eq(&n, &a));
    expect_dec_eq(&rem, "42");

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&n);
    cfx_big_free(&rem);
    cfx_big_free(&forty_two);
}

static void test_big_div_quotient_only_and_remainder_only(void) {
    cfx_big_t a, b, n, q, r, b_minus_1;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&n);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_init(&b_minus_1);
    cfx_big_from_dec(&a, "3141592653589793238462643383279502884197");
    cfx_big_from_dec(&b, "2718281828459045235360287471352662497757");

    cfx_big_copy(&n, &a);
    cfx_big_mul_eq(&n, &b);

    cfx_big_copy(&b_minus_1, &b);
    big_dec1(&b_minus_1);

    cfx_big_add_eq(&n, &b_minus_1);

    CFX_ASSERT(cfx_big_div(&q, &n, &b) == 0);
    CFX_ASSERT(cfx_big_eq(&q, &a));

    CFX_ASSERT(cfx_big_mod(&r, &n, &b) == 0);
    CFX_ASSERT(cfx_big_eq(&r, &b_minus_1));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&n);
    cfx_big_free(&q);
    cfx_big_free(&r);
    cfx_big_free(&b_minus_1);
}

static void test_big_div_alias_remainder_eq_src(void) {
    CFX_ASSERT(1);
    /* TODO: test aliasing cases */
}

int main(void) {
    CFX_TEST(test_big_div_divide_by_zero);
    CFX_TEST(test_big_div_zero_dividend);
    CFX_TEST(test_big_div_n_less_than_d);
    CFX_TEST(test_big_div_equal_numbers);
    CFX_TEST(test_big_div_single_limb_divisor_property);
    CFX_TEST(test_big_div_multi_limb_divisor_exact_and_remainder);
    CFX_TEST(test_big_div_in_place_eq_with_remainder);
    CFX_TEST(test_big_div_quotient_only_and_remainder_only);
    CFX_TEST(test_big_div_alias_remainder_eq_src);
    puts("OK");
    return 0;
}
