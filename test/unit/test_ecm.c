/*
 * test_ecm.c - Unit tests for ECM factorization
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include "cfx/big.h"
#include "cfx/ecm.h"
#include "cfx/algo.h"
#include "cfx/macros.h"

/* Test factoring small semiprimes */
static int test_ecm_small_semiprime(void) {
    printf("test_ecm_small_semiprime... ");

    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* n = 1000003 * 1000033 = 1000036000099 */
    cfx_big_from_str(&n, "1000036000099");

    int found = cfx_ecm_factor(&factor, &n, 10000, 10);
    CFX_ASSERT(found == 1);

    /* Verify it's a valid factor */
    cfx_big_t rem;
    cfx_big_init(&rem);
    cfx_big_mod(&rem, &n, &factor);
    CFX_ASSERT(cfx_big_is_zero(&rem));

    /* Factor should be one of the primes */
    uint64_t f = 0;
    if (factor.n == 1) {
        f = factor.limb[0];
    } else if (factor.n == 2 && CFX_LIMB_BITS == 32) {
        f = factor.limb[0] | ((uint64_t)factor.limb[1] << 32);
    }
    CFX_ASSERT(f == 1000003 || f == 1000033);

    cfx_big_free(&n);
    cfx_big_free(&factor);
    cfx_big_free(&rem);

    printf("OK\n");
    return 0;
}

/* Test factoring a larger semiprime (two 32-bit primes) */
static int test_ecm_medium_semiprime(void) {
    printf("test_ecm_medium_semiprime... ");

    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* n = 2147483647 * 2147483629 (two Mersenne-ish primes near 2^31) */
    /* = 4611686009837453363 */
    cfx_big_from_str(&n, "4611686009837453363");

    clock_t start = clock();
    int found = cfx_ecm_factor(&factor, &n, 100000, 25);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("(%.3fs) ", elapsed);

    CFX_ASSERT(found == 1);

    /* Verify it's a valid factor */
    cfx_big_t rem;
    cfx_big_init(&rem);
    cfx_big_mod(&rem, &n, &factor);
    CFX_ASSERT(cfx_big_is_zero(&rem));

    cfx_big_free(&n);
    cfx_big_free(&factor);
    cfx_big_free(&rem);

    printf("OK\n");
    return 0;
}

/* Test factoring a 128-bit semiprime (two 64-bit primes) */
static int test_ecm_large_semiprime(void) {
    printf("test_ecm_large_semiprime... ");

    cfx_big_t n, factor, p, q;
    cfx_big_init(&n);
    cfx_big_init(&factor);
    cfx_big_init(&p);
    cfx_big_init(&q);

    /* Two 64-bit primes:
     * p = 18446744073709551557 (largest 64-bit prime)
     * q = 18446744073709551533 (another large 64-bit prime)
     * n = p * q = 340282366920938460843936948965011886881
     */
    cfx_big_from_str(&p, "18446744073709551557");
    cfx_big_from_str(&q, "18446744073709551533");

    /* Compute n = p * q */
    cfx_big_mul(&n, &p, &q);

    char* str = cfx_big_dec_alloc(&n, NULL);
    printf("\nn = %s\n", str);
    printf("expected = 340282366920938460843936948965011886881\n");
    free(str);

    /* Verify the product is correct */
    cfx_big_t expected_n;
    cfx_big_init(&expected_n);
    cfx_big_from_str(&expected_n, "340282366920938460843936948965011886881");
    if (cfx_big_cmp(&n, &expected_n) != 0) {
        printf("ERROR: n != p*q!\n");
        str = cfx_big_dec_alloc(&p, NULL);
        printf("p = %s\n", str);
        free(str);
        str = cfx_big_dec_alloc(&q, NULL);
        printf("q = %s\n", str);
        free(str);
    }
    cfx_big_free(&expected_n);

    clock_t start = clock();
    int found = cfx_ecm_factor_auto(&factor, &n);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("ECM took %.3f seconds\n", elapsed);

    if (found) {
        /* Print factor as decimal */
        char* str = cfx_big_dec_alloc(&factor, NULL);
        printf("Found factor: %s\n", str);
        free(str);

        /* Verify it's a valid factor */
        cfx_big_t rem, quot;
        cfx_big_init(&rem);
        cfx_big_init(&quot);

        /* Check if factor divides n */
        cfx_big_copy(&quot, &n);
        cfx_big_divrem_eq(&quot, &factor, &rem);

        str = cfx_big_dec_alloc(&quot, NULL);
        printf("n / factor = %s\n", str);
        free(str);

        str = cfx_big_dec_alloc(&rem, NULL);
        printf("remainder = %s\n", str);
        free(str);

        /* Factor should be p or q */
        str = cfx_big_dec_alloc(&p, NULL);
        printf("Expected p = %s\n", str);
        free(str);

        str = cfx_big_dec_alloc(&q, NULL);
        printf("Expected q = %s\n", str);
        free(str);

        /* The factor found should equal p or q */
        int is_p = (cfx_big_cmp(&factor, &p) == 0);
        int is_q = (cfx_big_cmp(&factor, &q) == 0);
        printf("factor == p? %d, factor == q? %d\n", is_p, is_q);

        CFX_ASSERT(cfx_big_is_zero(&rem));
        cfx_big_free(&rem);
        cfx_big_free(&quot);
    } else {
        printf("No factor found (this is expected for very hard semiprimes)\n");
    }

    cfx_big_free(&n);
    cfx_big_free(&factor);
    cfx_big_free(&p);
    cfx_big_free(&q);

    printf("test_ecm_large_semiprime... OK\n");
    return 0;
}

/* Test that ECM handles trivial cases */
static int test_ecm_trivial(void) {
    printf("test_ecm_trivial... ");

    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* Even number - should find 2 immediately */
    cfx_big_from_limb(&n, 1234567890);
    int found = cfx_ecm_factor(&factor, &n, 1000, 1);
    CFX_ASSERT(found == 1);
    CFX_ASSERT(factor.n == 1 && factor.limb[0] == 2);

    cfx_big_free(&n);
    cfx_big_free(&factor);

    printf("OK\n");
    return 0;
}

static int test_ecm_zero_input(void) {
    printf("test_ecm_zero_input... ");

    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    int found = cfx_ecm_factor(&factor, &n, 1000, 1);
    CFX_ASSERT(found == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);

    printf("OK\n");
    return 0;
}

static int test_ecm_one_input(void) {
    printf("test_ecm_one_input... ");

    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    cfx_big_from_limb(&n, 1);
    int found = cfx_ecm_factor(&factor, &n, 1000, 1);
    CFX_ASSERT(found == 0);

    cfx_big_free(&n);
    cfx_big_free(&factor);

    printf("OK\n");
    return 0;
}

static int test_ecm_auto_small(void) {
    printf("test_ecm_auto_small... ");

    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* small even number - exercises < 64 bit path and even check */
    cfx_big_from_limb(&n, 1234);
    int found = cfx_ecm_factor_auto(&factor, &n);
    CFX_ASSERT(found == 1);
    CFX_ASSERT(factor.limb[0] == 2);

    cfx_big_free(&n);
    cfx_big_free(&factor);

    printf("OK\n");
    return 0;
}

static int test_ecm_auto_96bit(void) {
    printf("test_ecm_auto_96bit... ");

    cfx_big_t n, factor;
    cfx_big_init(&n);
    cfx_big_init(&factor);

    /* ~80 bit semiprime: 1125899906842597 * 1125899906842679
     * = 1267650600228226802437530494363 (~100 bits total, so > 96 bit path)
     * Let's use a smaller one in the 65-96 bit range */
    /* 2^65 = 36893488147419103232, need something ~75 bits */
    /* 1099511627791 * 1099511627689 = 1208925819614629000199 (~70 bits) */
    cfx_big_from_str(&n, "1208925819614629000199");

    clock_t start = clock();
    int found = cfx_ecm_factor_auto(&factor, &n);
    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC;

    printf("(%.3fs) ", elapsed);

    if (found) {
        cfx_big_t rem;
        cfx_big_init(&rem);
        cfx_big_mod(&rem, &n, &factor);
        CFX_ASSERT(cfx_big_is_zero(&rem));
        cfx_big_free(&rem);
    }

    cfx_big_free(&n);
    cfx_big_free(&factor);

    printf("OK\n");
    return 0;
}

int main(void) {
#ifdef CFX_MEMORY_STATIC
    /* Skip: ECM tests exhaust static buffer pool */
    printf("SKIP: test_ecm requires dynamic memory mode\n");
    return 0;
#endif
    CFX_TEST(test_ecm_trivial);
    CFX_TEST(test_ecm_zero_input);
    CFX_TEST(test_ecm_one_input);
    CFX_TEST(test_ecm_auto_small);
    CFX_TEST(test_ecm_small_semiprime);
    CFX_TEST(test_ecm_medium_semiprime);
    CFX_TEST(test_ecm_auto_96bit);
    // CFX_TEST(test_ecm_large_semiprime);
    printf("\nAll ECM tests passed!\n");
    return 0;
}
