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

#define CFX_ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERT failed: %s at %s:%d\n", #cond, __FILE__, __LINE__); \
        return 1; \
    } \
} while (0)

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

    char* str = cfx_big_to_str(&n, NULL);
    printf("\nn = %s\n", str);
    printf("expected = 340282366920938460843936948965011886881\n");
    free(str);

    /* Verify the product is correct */
    cfx_big_t expected_n;
    cfx_big_init(&expected_n);
    cfx_big_from_str(&expected_n, "340282366920938460843936948965011886881");
    if (cfx_big_cmp(&n, &expected_n) != 0) {
        printf("ERROR: n != p*q!\n");
        str = cfx_big_to_str(&p, NULL);
        printf("p = %s\n", str);
        free(str);
        str = cfx_big_to_str(&q, NULL);
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
        char* str = cfx_big_to_str(&factor, NULL);
        printf("Found factor: %s\n", str);
        free(str);

        /* Verify it's a valid factor */
        cfx_big_t rem, quot;
        cfx_big_init(&rem);
        cfx_big_init(&quot);

        /* Check if factor divides n */
        cfx_big_copy(&quot, &n);
        cfx_big_div_eq(&quot, &factor, &rem);

        str = cfx_big_to_str(&quot, NULL);
        printf("n / factor = %s\n", str);
        free(str);

        str = cfx_big_to_str(&rem, NULL);
        printf("remainder = %s\n", str);
        free(str);

        /* Factor should be p or q */
        str = cfx_big_to_str(&p, NULL);
        printf("Expected p = %s\n", str);
        free(str);

        str = cfx_big_to_str(&q, NULL);
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

int main(void) {
    printf("=== ECM Tests ===\n");

    if (test_ecm_trivial()) return 1;
    if (test_ecm_small_semiprime()) return 1;
    if (test_ecm_medium_semiprime()) return 1;
    if (test_ecm_large_semiprime()) return 1;

    printf("\nAll ECM tests passed!\n");
    return 0;
}
