/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_pollard_rho.c - Tests for Pollard-Rho factorization */

#include "test_common.h"

static void test_big_pollard_rho_small_composite(void) {
    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* n = 15 = 3 * 5 */
    cfx_big_from_limb(&n, 15);
    cfx_big_pollard_rho(&factor, &n);
    /* Should find either 3 or 5 */
    CFX_ASSERT(factor.n == 1);
    CFX_ASSERT(factor.limb[0] == 3 || factor.limb[0] == 5);

    /* n = 21 = 3 * 7 */
    cfx_big_from_limb(&n, 21);
    cfx_big_pollard_rho(&factor, &n);
    CFX_ASSERT(factor.n == 1);
    CFX_ASSERT(factor.limb[0] == 3 || factor.limb[0] == 7);

    cfx_big_free(&n);
    cfx_big_free(&factor);
}

static void test_big_pollard_rho_semiprime(void) {
    cfx_big_t n, factor, quotient, check;
    cfx_big_init(&n);
    cfx_big_init(&factor);
    cfx_big_init(&quotient);
    cfx_big_init(&check);

    /* n = 1000003 * 1000033 = 1000036000099 (semiprime) */
    cfx_big_from_str(&n, "1000036000099");
    cfx_big_pollard_rho(&factor, &n);

    /* factor should be non-trivial (not 1 and not n) */
    CFX_ASSERT(!cfx_big_is_one(&factor));
    CFX_ASSERT(cfx_big_cmp(&factor, &n) != 0);

    /* Verify: n should be divisible by factor */
    cfx_big_copy(&quotient, &n);
    cfx_big_divrem_eq(&quotient, &factor, NULL);
    cfx_big_mul(&check, &quotient, &factor);
    CFX_ASSERT(cfx_big_cmp(&check, &n) == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);
    cfx_big_free(&quotient);
    cfx_big_free(&check);
}

static void test_big_pollard_rho_prime(void) {
    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* n = 104729 (prime) - Pollard-Rho should return n itself */
    cfx_big_from_limb(&n, 104729);
    cfx_big_pollard_rho(&factor, &n);
    /* For a prime, it should return the prime itself (can't factor) */
    CFX_ASSERT(cfx_big_cmp(&factor, &n) == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);
}

static void test_big_pollard_rho_even(void) {
    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* Even number should return 2 */
    cfx_big_from_limb(&n, 1234567890);
    cfx_big_pollard_rho(&factor, &n);
    CFX_ASSERT(factor.n == 1 && factor.limb[0] == 2);

    cfx_big_free(&n);
    cfx_big_free(&factor);
}

static void test_big_pollard_rho_large_semiprime(void) {
    cfx_big_t n, factor, quotient, check;
    cfx_big_init(&n);
    cfx_big_init(&factor);
    cfx_big_init(&quotient);
    cfx_big_init(&check);

    /* n = 10007 * 10009 = 100160063 */
    cfx_big_from_str(&n, "100160063");
    cfx_big_pollard_rho(&factor, &n);

    /* factor should be non-trivial */
    CFX_ASSERT(!cfx_big_is_one(&factor));
    CFX_ASSERT(cfx_big_cmp(&factor, &n) != 0);

    /* Verify divisibility */
    cfx_big_copy(&quotient, &n);
    cfx_big_divrem_eq(&quotient, &factor, NULL);
    cfx_big_mul(&check, &quotient, &factor);
    CFX_ASSERT(cfx_big_cmp(&check, &n) == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);
    cfx_big_free(&quotient);
    cfx_big_free(&check);
}

int main(void) {
    CFX_TEST(test_big_pollard_rho_small_composite);
    CFX_TEST(test_big_pollard_rho_semiprime);
    CFX_TEST(test_big_pollard_rho_prime);
    CFX_TEST(test_big_pollard_rho_even);
    CFX_TEST(test_big_pollard_rho_large_semiprime);
    puts("OK");
    return 0;
}
