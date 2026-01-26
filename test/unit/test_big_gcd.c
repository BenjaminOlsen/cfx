/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_gcd.c - Tests for big integer GCD */

#include "test_common.h"

static void test_big_gcd_basic(void) {
    cfx_big_t a, b, g;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&g);

    /* gcd(12, 8) = 4 */
    cfx_big_from_limb(&a, 12);
    cfx_big_from_limb(&b, 8);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 4);

    /* gcd(17, 13) = 1 (coprime) */
    cfx_big_from_limb(&a, 17);
    cfx_big_from_limb(&b, 13);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 1);

    /* gcd(0, 5) = 5 */
    cfx_big_from_limb(&a, 0);
    cfx_big_from_limb(&b, 5);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 5);

    /* gcd(100, 0) = 100 */
    cfx_big_from_limb(&a, 100);
    cfx_big_from_limb(&b, 0);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 100);

    /* gcd(48, 18) = 6 */
    cfx_big_from_limb(&a, 48);
    cfx_big_from_limb(&b, 18);
    cfx_big_gcd(&g, &a, &b);
    CFX_ASSERT(g.n == 1 && g.limb[0] == 6);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&g);
}

static void test_big_gcd_large(void) {
    cfx_big_t a, b, g;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&g);

    /* a = 123456789 * 1000, b = 123456789 * 777 => gcd = 123456789 */
    cfx_big_from_str(&a, "123456789000");
    cfx_big_from_str(&b, "95925925053");
    cfx_big_gcd(&g, &a, &b);

    cfx_big_t expected;
    cfx_big_init(&expected);
    cfx_big_from_limb(&expected, 123456789);
    CFX_ASSERT(cfx_big_cmp(&g, &expected) == 0);

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&g);
    cfx_big_free(&expected);
}

static void test_modinv_basic(void) {
    cfx_big_t a, n, inv, t;
    cfx_big_init(&a);
    cfx_big_init(&n);
    cfx_big_init(&inv);
    cfx_big_init(&t);

    /* 3^(-1) mod 7 = 5  (since 3*5 = 15 = 2*7 + 1) */
    cfx_big_from_limb(&a, 3);
    cfx_big_from_limb(&n, 7);
    CFX_ASSERT(cfx_big_modinv(&inv, &a, &n) == 1);
    CFX_ASSERT(cfx_big_cmp_u64(&inv, 5) == 0);

    /* verify: a * inv mod n == 1 */
    cfx_big_mul(&t, &a, &inv);
    cfx_big_mod(&t, &t, &n);
    CFX_ASSERT(cfx_big_is_one(&t));

    /* 2^(-1) mod 11 = 6  (since 2*6 = 12 = 11 + 1) */
    cfx_big_from_limb(&a, 2);
    cfx_big_from_limb(&n, 11);
    CFX_ASSERT(cfx_big_modinv(&inv, &a, &n) == 1);
    cfx_big_mul(&t, &a, &inv);
    cfx_big_mod(&t, &t, &n);
    CFX_ASSERT(cfx_big_is_one(&t));

    cfx_big_free(&a);
    cfx_big_free(&n);
    cfx_big_free(&inv);
    cfx_big_free(&t);
    PRINT_TEST(1);
}

static void test_modinv_no_inverse(void) {
    cfx_big_t a, n, inv;
    cfx_big_init(&a);
    cfx_big_init(&n);
    cfx_big_init(&inv);

    /* 4^(-1) mod 8 doesn't exist (gcd = 4) */
    cfx_big_from_limb(&a, 4);
    cfx_big_from_limb(&n, 8);
    CFX_ASSERT(cfx_big_modinv(&inv, &a, &n) == 0);

    /* 6^(-1) mod 9 doesn't exist (gcd = 3) */
    cfx_big_from_limb(&a, 6);
    cfx_big_from_limb(&n, 9);
    CFX_ASSERT(cfx_big_modinv(&inv, &a, &n) == 0);

    cfx_big_free(&a);
    cfx_big_free(&n);
    cfx_big_free(&inv);
    PRINT_TEST(1);
}

static void test_modinv_large(void) {
    cfx_big_t a, n, inv, t;
    cfx_big_init(&a);
    cfx_big_init(&n);
    cfx_big_init(&inv);
    cfx_big_init(&t);

    /* larger numbers - use a prime modulus */
    cfx_big_from_str(&a, "123456789");
    cfx_big_from_str(&n, "1000000007");  /* prime */
    CFX_ASSERT(cfx_big_modinv(&inv, &a, &n) == 1);

    /* verify: a * inv mod n == 1 */
    cfx_big_mul(&t, &a, &inv);
    cfx_big_mod(&t, &t, &n);
    CFX_ASSERT(cfx_big_is_one(&t));

    cfx_big_free(&a);
    cfx_big_free(&n);
    cfx_big_free(&inv);
    cfx_big_free(&t);
    PRINT_TEST(1);
}

int main(void) {
    CFX_TEST(test_big_gcd_basic);
    CFX_TEST(test_big_gcd_large);
    CFX_TEST(test_modinv_basic);
    CFX_TEST(test_modinv_no_inverse);
    CFX_TEST(test_modinv_large);
    puts("OK");
    return 0;
}
