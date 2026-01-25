/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "test_common.h"
#include "cfx/sbig.h"

static void test_init_free(void) {
    cfx_sbig_t s;
    cfx_sbig_init(&s);
    CFX_ASSERT(cfx_sbig_is_zero(&s));
    CFX_ASSERT(cfx_sbig_sign(&s) == 0);
    cfx_sbig_free(&s);
    PRINT_TEST(1);
}

static void test_from_i64(void) {
    cfx_sbig_t s;
    cfx_sbig_init(&s);

    /* positive */
    cfx_sbig_from_i64(&s, 12345);
    CFX_ASSERT(cfx_sbig_sign(&s) == 1);
    CFX_ASSERT(cfx_sbig_is_positive(&s));
    CFX_ASSERT(cfx_big_eq_u64(&s.mag, 12345));

    /* zero */
    cfx_sbig_from_i64(&s, 0);
    CFX_ASSERT(cfx_sbig_sign(&s) == 0);
    CFX_ASSERT(cfx_sbig_is_zero(&s));

    /* negative */
    cfx_sbig_from_i64(&s, -42);
    CFX_ASSERT(cfx_sbig_sign(&s) == -1);
    CFX_ASSERT(cfx_sbig_is_negative(&s));
    CFX_ASSERT(cfx_big_eq_u64(&s.mag, 42));

    /* INT64_MIN edge case: |INT64_MIN| = 2^63 = 0x8000000000000000 */
    cfx_sbig_from_i64(&s, INT64_MIN);
    CFX_ASSERT(cfx_sbig_sign(&s) == -1);
    /* Can't use cfx_big_eq_u64 for 64-bit value on 32-bit limbs */
    cfx_big_t expected;
    cfx_big_init(&expected);
    cfx_big_from_hex(&expected, "8000000000000000");
    CFX_ASSERT(cfx_big_eq(&s.mag, &expected));
    cfx_big_free(&expected);

    cfx_sbig_free(&s);
    PRINT_TEST(1);
}

static void test_neg(void) {
    cfx_sbig_t s;
    cfx_sbig_init(&s);

    cfx_sbig_from_i64(&s, 100);
    CFX_ASSERT(cfx_sbig_sign(&s) == 1);

    cfx_sbig_neg_eq(&s);
    CFX_ASSERT(cfx_sbig_sign(&s) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&s.mag, 100));

    cfx_sbig_neg_eq(&s);
    CFX_ASSERT(cfx_sbig_sign(&s) == 1);

    /* negating zero stays zero */
    cfx_sbig_from_i64(&s, 0);
    cfx_sbig_neg_eq(&s);
    CFX_ASSERT(cfx_sbig_is_zero(&s));

    cfx_sbig_free(&s);
    PRINT_TEST(1);
}

static void test_abs(void) {
    cfx_sbig_t s;
    cfx_sbig_init(&s);

    cfx_sbig_from_i64(&s, -999);
    CFX_ASSERT(cfx_sbig_sign(&s) == -1);

    cfx_sbig_abs_eq(&s);
    CFX_ASSERT(cfx_sbig_sign(&s) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&s.mag, 999));

    /* abs of positive stays positive */
    cfx_sbig_abs_eq(&s);
    CFX_ASSERT(cfx_sbig_sign(&s) == 1);

    cfx_sbig_free(&s);
    PRINT_TEST(1);
}

static void test_cmp(void) {
    cfx_sbig_t a, b;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);

    /* positive vs positive */
    cfx_sbig_from_i64(&a, 100);
    cfx_sbig_from_i64(&b, 50);
    CFX_ASSERT(cfx_sbig_cmp(&a, &b) == 1);
    CFX_ASSERT(cfx_sbig_cmp(&b, &a) == -1);
    cfx_sbig_from_i64(&b, 100);
    CFX_ASSERT(cfx_sbig_cmp(&a, &b) == 0);

    /* negative vs negative */
    cfx_sbig_from_i64(&a, -100);
    cfx_sbig_from_i64(&b, -50);
    CFX_ASSERT(cfx_sbig_cmp(&a, &b) == -1);  /* -100 < -50 */
    CFX_ASSERT(cfx_sbig_cmp(&b, &a) == 1);

    /* positive vs negative */
    cfx_sbig_from_i64(&a, 1);
    cfx_sbig_from_i64(&b, -1000);
    CFX_ASSERT(cfx_sbig_cmp(&a, &b) == 1);  /* 1 > -1000 */

    /* zero comparisons */
    cfx_sbig_from_i64(&a, 0);
    cfx_sbig_from_i64(&b, 0);
    CFX_ASSERT(cfx_sbig_cmp(&a, &b) == 0);

    cfx_sbig_from_i64(&b, 1);
    CFX_ASSERT(cfx_sbig_cmp(&a, &b) == -1);  /* 0 < 1 */

    cfx_sbig_from_i64(&b, -1);
    CFX_ASSERT(cfx_sbig_cmp(&a, &b) == 1);   /* 0 > -1 */

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    PRINT_TEST(1);
}

static void test_add_same_sign(void) {
    cfx_sbig_t a, b, out;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);
    cfx_sbig_init(&out);

    /* positive + positive */
    cfx_sbig_from_i64(&a, 100);
    cfx_sbig_from_i64(&b, 50);
    cfx_sbig_add(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 150));

    /* negative + negative */
    cfx_sbig_from_i64(&a, -100);
    cfx_sbig_from_i64(&b, -50);
    cfx_sbig_add(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 150));

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    cfx_sbig_free(&out);
    PRINT_TEST(1);
}

static void test_add_diff_sign(void) {
    cfx_sbig_t a, b, out;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);
    cfx_sbig_init(&out);

    /* positive + negative, |a| > |b| */
    cfx_sbig_from_i64(&a, 100);
    cfx_sbig_from_i64(&b, -30);
    cfx_sbig_add(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 70));

    /* positive + negative, |a| < |b| */
    cfx_sbig_from_i64(&a, 30);
    cfx_sbig_from_i64(&b, -100);
    cfx_sbig_add(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 70));

    /* cancellation to zero */
    cfx_sbig_from_i64(&a, 42);
    cfx_sbig_from_i64(&b, -42);
    cfx_sbig_add(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_is_zero(&out));

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    cfx_sbig_free(&out);
    PRINT_TEST(1);
}

static void test_add_zero(void) {
    cfx_sbig_t a, b, out;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);
    cfx_sbig_init(&out);

    cfx_sbig_from_i64(&a, 0);
    cfx_sbig_from_i64(&b, 42);
    cfx_sbig_add(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 42));

    cfx_sbig_add(&out, &b, &a);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 42));

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    cfx_sbig_free(&out);
    PRINT_TEST(1);
}

static void test_sub(void) {
    cfx_sbig_t a, b, out;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);
    cfx_sbig_init(&out);

    /* 100 - 30 = 70 */
    cfx_sbig_from_i64(&a, 100);
    cfx_sbig_from_i64(&b, 30);
    cfx_sbig_sub(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 70));

    /* 30 - 100 = -70 */
    cfx_sbig_from_i64(&a, 30);
    cfx_sbig_from_i64(&b, 100);
    cfx_sbig_sub(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 70));

    /* -30 - (-100) = 70 */
    cfx_sbig_from_i64(&a, -30);
    cfx_sbig_from_i64(&b, -100);
    cfx_sbig_sub(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 70));

    /* 50 - 50 = 0 */
    cfx_sbig_from_i64(&a, 50);
    cfx_sbig_from_i64(&b, 50);
    cfx_sbig_sub(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_is_zero(&out));

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    cfx_sbig_free(&out);
    PRINT_TEST(1);
}

static void test_mul(void) {
    cfx_sbig_t a, b, out;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);
    cfx_sbig_init(&out);

    /* positive * positive = positive */
    cfx_sbig_from_i64(&a, 12);
    cfx_sbig_from_i64(&b, 11);
    cfx_sbig_mul(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 132));

    /* positive * negative = negative */
    cfx_sbig_from_i64(&a, 12);
    cfx_sbig_from_i64(&b, -11);
    cfx_sbig_mul(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 132));

    /* negative * negative = positive */
    cfx_sbig_from_i64(&a, -12);
    cfx_sbig_from_i64(&b, -11);
    cfx_sbig_mul(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&out) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&out.mag, 132));

    /* anything * zero = zero */
    cfx_sbig_from_i64(&a, 999);
    cfx_sbig_from_i64(&b, 0);
    cfx_sbig_mul(&out, &a, &b);
    CFX_ASSERT(cfx_sbig_is_zero(&out));

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    cfx_sbig_free(&out);
    PRINT_TEST(1);
}

static void test_divrem(void) {
    cfx_sbig_t a, b, q, r;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);
    cfx_sbig_init(&q);
    cfx_sbig_init(&r);

    /* 17 / 5 = 3 remainder 2 */
    cfx_sbig_from_i64(&a, 17);
    cfx_sbig_from_i64(&b, 5);
    cfx_sbig_divrem(&q, &r, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&q) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&q.mag, 3));
    CFX_ASSERT(cfx_sbig_sign(&r) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&r.mag, 2));

    /* -17 / 5 = -3 remainder -2 (C99 semantics) */
    cfx_sbig_from_i64(&a, -17);
    cfx_sbig_from_i64(&b, 5);
    cfx_sbig_divrem(&q, &r, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&q) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&q.mag, 3));
    CFX_ASSERT(cfx_sbig_sign(&r) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&r.mag, 2));

    /* 17 / -5 = -3 remainder 2 */
    cfx_sbig_from_i64(&a, 17);
    cfx_sbig_from_i64(&b, -5);
    cfx_sbig_divrem(&q, &r, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&q) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&q.mag, 3));
    CFX_ASSERT(cfx_sbig_sign(&r) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&r.mag, 2));

    /* -17 / -5 = 3 remainder -2 */
    cfx_sbig_from_i64(&a, -17);
    cfx_sbig_from_i64(&b, -5);
    cfx_sbig_divrem(&q, &r, &a, &b);
    CFX_ASSERT(cfx_sbig_sign(&q) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&q.mag, 3));
    CFX_ASSERT(cfx_sbig_sign(&r) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&r.mag, 2));

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    cfx_sbig_free(&q);
    cfx_sbig_free(&r);
    PRINT_TEST(1);
}

static void test_string_io(void) {
    cfx_sbig_t s;
    cfx_sbig_init(&s);

    /* positive */
    cfx_sbig_from_dec(&s, "12345");
    CFX_ASSERT(cfx_sbig_sign(&s) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&s.mag, 12345));

    char *str = cfx_sbig_dec_alloc(&s, NULL);
    CFX_ASSERT(strcmp(str, "12345") == 0);
    free(str);

    /* negative */
    cfx_sbig_from_dec(&s, "-98765");
    CFX_ASSERT(cfx_sbig_sign(&s) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&s.mag, 98765));

    str = cfx_sbig_dec_alloc(&s, NULL);
    CFX_ASSERT(strcmp(str, "-98765") == 0);
    free(str);

    /* zero */
    cfx_sbig_from_dec(&s, "0");
    CFX_ASSERT(cfx_sbig_is_zero(&s));
    str = cfx_sbig_dec_alloc(&s, NULL);
    CFX_ASSERT(strcmp(str, "0") == 0);
    free(str);

    /* hex with sign */
    cfx_sbig_from_hex(&s, "-ff");
    CFX_ASSERT(cfx_sbig_sign(&s) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&s.mag, 255));

    cfx_sbig_free(&s);
    PRINT_TEST(1);
}

static void test_in_place_ops(void) {
    cfx_sbig_t a, b;
    cfx_sbig_init(&a);
    cfx_sbig_init(&b);

    /* add_eq */
    cfx_sbig_from_i64(&a, 100);
    cfx_sbig_from_i64(&b, -30);
    cfx_sbig_add_eq(&a, &b);
    CFX_ASSERT(cfx_sbig_sign(&a) == 1);
    CFX_ASSERT(cfx_big_eq_u64(&a.mag, 70));

    /* sub_eq */
    cfx_sbig_from_i64(&a, 50);
    cfx_sbig_from_i64(&b, 80);
    cfx_sbig_sub_eq(&a, &b);
    CFX_ASSERT(cfx_sbig_sign(&a) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&a.mag, 30));

    /* mul_eq */
    cfx_sbig_from_i64(&a, -7);
    cfx_sbig_from_i64(&b, 6);
    cfx_sbig_mul_eq(&a, &b);
    CFX_ASSERT(cfx_sbig_sign(&a) == -1);
    CFX_ASSERT(cfx_big_eq_u64(&a.mag, 42));

    cfx_sbig_free(&a);
    cfx_sbig_free(&b);
    PRINT_TEST(1);
}

int main(void) {
    CFX_TEST(test_init_free);
    CFX_TEST(test_from_i64);
    CFX_TEST(test_neg);
    CFX_TEST(test_abs);
    CFX_TEST(test_cmp);
    CFX_TEST(test_add_same_sign);
    CFX_TEST(test_add_diff_sign);
    CFX_TEST(test_add_zero);
    CFX_TEST(test_sub);
    CFX_TEST(test_mul);
    CFX_TEST(test_divrem);
    CFX_TEST(test_string_io);
    CFX_TEST(test_in_place_ops);
    puts("OK");
    return 0;
}
