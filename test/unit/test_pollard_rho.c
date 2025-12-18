/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/algo.h"
#include "cfx/arith.h"
#include "cfx/macros.h"

#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>


static int is_valid_factor(uint64_t n, uint64_t d) {
    return (d > 1) && (d < n) && (n % d == 0);
}

static void expect_factor(uint64_t n) {
    /* Rho is only defined/useful for composites; guard primes */
    if (cfx_is_prime_u64(n)) {
        /* For primes, cfx_pollard_rho_brent may return n or 0, but we don't require anything. */
        /* Just ensure it *doesn't* falsely report a composite factor. */
        uint64_t d = cfx_pollard_rho_brent(n);
        if (is_valid_factor(n, d)) {
            fprintf(stderr,
                "cfx_pollard_rho_brent returned a nontrivial factor for prime %" PRIu64 "\n",
                n);
            CFX_FAIL();
        }
        return;
    }

    srand(123456u);

    /* Try a few times in case the internal random choices hit a bad cycle */
    for (int attempts = 0; attempts < 5; ++attempts) {
        uint64_t d = cfx_pollard_rho_brent(n);
        printf("cfx_pollard_rho_brent(%" PRIu64 ") = %" PRIu64 "\n", n, d);
        if (is_valid_factor(n, d)) {
            /* ? verify cofactor primality or at least correctness */
            uint64_t m = n / d;
            CFX_ASSERT(n == d * m);
            /* Either side may still be composite; that's fine (we only test rho's split) */
            return;
        }
        srand(123456u + (unsigned)attempts + 1);
    }
    fprintf(stderr, "cfx_pollard_rho_brent failed to find a factor for %" PRIu64 "\n",
            n);
    CFX_ASSERT(0);
}

static void test_small_semiprimes(void) {
    expect_factor(91ULL);
    expect_factor(221ULL);
    expect_factor(10403ULL);
    expect_factor(11021ULL);
}

static void test_carmichael_and_square_semiprimes(void) {
    expect_factor(561ULL);
    expect_factor(1009ULL * 1009ULL);
}

static void test_repeated_prime_powers(void) {
    uint64_t n1 = 1;
    for (int i = 0; i < 10; ++i) n1 *= 3ULL;
    expect_factor(n1);

    uint64_t n2 = UINT64_C(1) << 40;
    expect_factor(n2);
}

static void test_large_64bit_semiprime(void) {
    const uint64_t p = 4294967291ULL; /* 2^32 - 5, prime */
    const uint64_t q = 4294967279ULL; /* 2^32 - 17, prime */

    uint64_t n = p * q; /* still < 2^64 */

    expect_factor(n);
}

static void test_primes_do_not_yield_factors(void) {
    /* These should NOT yield nontrivial factors */
    CFX_ASSERT(cfx_is_prime_u64(29));
    CFX_ASSERT(cfx_is_prime_u64(97));
    CFX_ASSERT(cfx_is_prime_u64(257));
    CFX_ASSERT(cfx_is_prime_u64(65537));

    /* cfx_pollard_rho_brent may return 0 or n — both are acceptable "no factor" signals. */
    /* We only CFX_ASSERT it does NOT return a valid factor. */
    uint64_t primes[] = {29, 97, 257, 65537};
    for (size_t i = 0; i < sizeof(primes)/sizeof(primes[0]); ++i) {
        srand(42u);
        uint64_t d = cfx_pollard_rho_brent(primes[i]);
        CFX_ASSERT(!is_valid_factor(primes[i], d));
    }
}

int main(void) {
    CFX_TEST(test_small_semiprimes);
    CFX_TEST(test_carmichael_and_square_semiprimes);
    CFX_TEST(test_repeated_prime_powers);
    CFX_TEST(test_large_64bit_semiprime);
    CFX_TEST(test_primes_do_not_yield_factors);

    puts("cfx_pollard_rho_brent tests: OK");
    return 0;
}
