/*
 * Debug test for division issue
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>

#include "cfx/big.h"
#include "cfx/arith.h"

static void print_limbs(const char* name, const cfx_big_t* b) {
    printf("%s: n=%zu, limbs = {", name, b->n);
    for (size_t i = 0; i < b->n; i++) {
        printf("0x%08" PRIx32 "%s", (uint32_t)b->limb[i], i+1 < b->n ? ", " : "");
    }
    printf("}\n");
}

int main(void) {
    printf("=== Division Debug Test ===\n");
    printf("CFX_LIMB_BITS = %d\n\n", CFX_LIMB_BITS);

    cfx_big_t n, p, q, rem;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&q);
    cfx_big_init(&rem);

    /* Simple test first: 15 / 3 = 5 */
    printf("=== Simple test: 15 / 3 ===\n");
    cfx_big_from_limb(&n, 15);
    cfx_big_from_limb(&p, 3);
    cfx_big_divrem(&q, &rem, &n, &p);
    printf("15 / 3 = %u, rem = %u (expected 5, 0)\n\n",
           (unsigned)q.limb[0], cfx_big_is_zero(&rem) ? 0 : (unsigned)rem.limb[0]);

    /* Test with 2-limb numbers: 0xFFFFFFFF * 2 / 2 */
    printf("=== Test: (2^32-1)*2 / 2 ===\n");
    cfx_big_from_str(&n, "8589934590");  /* (2^32-1)*2 */
    cfx_big_from_limb(&p, 2);
    cfx_big_divrem(&q, &rem, &n, &p);
    char* str = cfx_big_to_str(&q, NULL);
    printf("8589934590 / 2 = %s (expected 4294967295)\n", str);
    free(str);
    printf("rem = %s (expected 0)\n\n", cfx_big_is_zero(&rem) ? "0" : "nonzero");

    /* Test with 2-limb divisor */
    printf("=== Test: 2-limb divisor ===\n");
    cfx_big_from_str(&n, "18446744073709551616");  /* 2^64 */
    cfx_big_from_str(&p, "4294967296");            /* 2^32 */
    cfx_big_divrem(&q, &rem, &n, &p);
    str = cfx_big_to_str(&q, NULL);
    printf("2^64 / 2^32 = %s (expected 4294967296)\n", str);
    free(str);
    printf("rem = %s (expected 0)\n\n", cfx_big_is_zero(&rem) ? "0" : "nonzero");

    /* Test with high limb values - potential edge case */
    printf("=== Test: High limb edge case ===\n");
    /* n = 0xFFFFFFFF * 0xFFFFFFFF = 0xFFFFFFFE00000001, divide by 0xFFFFFFFF */
    cfx_big_from_str(&n, "18446744065119617025");  /* 0xFFFFFFFE00000001 */
    cfx_big_from_str(&p, "4294967295");            /* 0xFFFFFFFF */
    cfx_big_divrem(&q, &rem, &n, &p);
    str = cfx_big_to_str(&q, NULL);
    printf("0xFFFFFFFE00000001 / 0xFFFFFFFF = %s (expected 4294967295)\n", str);
    free(str);
    str = cfx_big_to_str(&rem, NULL);
    printf("rem = %s (expected 2)\n\n", str);
    free(str);

    /* Test: 2-limb by 2-limb division with max values */
    printf("=== Test: 2-limb by 2-limb max values ===\n");
    cfx_big_from_str(&n, "340282366920938463463374607431768211455");  /* 2^128 - 1 */
    cfx_big_from_str(&p, "18446744073709551615");                    /* 2^64 - 1 */
    cfx_big_divrem(&q, &rem, &n, &p);
    str = cfx_big_to_str(&q, NULL);
    printf("(2^128-1) / (2^64-1) = %s (expected 18446744073709551617)\n", str);
    free(str);
    str = cfx_big_to_str(&rem, NULL);
    printf("rem = %s (expected 2)\n\n", str);
    free(str);

    /* n = p * q where p and q are 64-bit primes */
    printf("=== Main test: large semiprime ===\n");
    cfx_big_from_str(&n, "340282366920938460843936948965011886881");
    cfx_big_from_str(&p, "18446744073709551557");

    print_limbs("n", &n);
    print_limbs("p", &p);

    printf("\nn = ");
    str = cfx_big_to_str(&n, NULL);
    printf("%s\n", str);
    free(str);

    printf("p = ");
    str = cfx_big_to_str(&p, NULL);
    printf("%s\n", str);
    free(str);

    /* First try using cfx_big_divrem directly to rule out cfx_big_div_eq */
    printf("\nUsing cfx_big_divrem...\n");
    cfx_big_divrem(&q, &rem, &n, &p);

    print_limbs("q", &q);
    print_limbs("rem", &rem);

    printf("q = n / p = ");
    str = cfx_big_to_str(&q, NULL);
    printf("%s\n", str);
    free(str);

    printf("rem = n %% p = ");
    str = cfx_big_to_str(&rem, NULL);
    printf("%s\n", str);
    free(str);

    printf("\nExpected q = 18446744073709551533\n");
    printf("Expected rem = 0\n");

    /* Verify */
    cfx_big_t expected_q;
    cfx_big_init(&expected_q);
    cfx_big_from_str(&expected_q, "18446744073709551533");

    int q_ok = (cfx_big_cmp(&q, &expected_q) == 0);
    int rem_ok = cfx_big_is_zero(&rem);

    printf("\nq correct: %s\n", q_ok ? "YES" : "NO");
    printf("rem correct: %s\n", rem_ok ? "YES" : "NO");

    /* Verify: q * p + rem should equal n */
    printf("\nVerification: q * p + rem should equal n\n");
    cfx_big_t check, n_orig;
    cfx_big_init(&check);
    cfx_big_init(&n_orig);
    cfx_big_from_str(&n_orig, "340282366920938460843936948965011886881");
    cfx_big_mul(&check, &q, &p);
    print_limbs("q * p", &check);
    cfx_big_add(&check, &rem);
    str = cfx_big_to_str(&check, NULL);
    printf("q * p + rem = %s\n", str);
    free(str);
    str = cfx_big_to_str(&n_orig, NULL);
    printf("n           = %s\n", str);
    free(str);
    int verify_ok = (cfx_big_cmp(&check, &n_orig) == 0);
    printf("Verification: %s\n", verify_ok ? "PASS" : "FAIL");
    cfx_big_free(&check);
    cfx_big_free(&n_orig);

    cfx_big_free(&n);
    cfx_big_free(&p);
    cfx_big_free(&q);
    cfx_big_free(&rem);
    cfx_big_free(&expected_q);

    return (q_ok && rem_ok) ? 0 : 1;
}
