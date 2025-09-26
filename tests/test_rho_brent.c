#include "cfx/algo.h"
#include "cfx/num.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>


static int is_valid_factor(cfx_limb_t n, cfx_limb_t d) {
    return (d > 1) && (d < n) && (n % d == 0);
}

static void expect_factor(cfx_limb_t n) {
    // Rho is only defined/useful for composites; guard primes
    if (cfx_is_prime_u64(n)) {
        // For primes, cfx_rho_brent may return n or 0, but we don't require anything.
        // Just ensure it *doesn't* falsely report a composite factor.
        cfx_limb_t d = cfx_rho_brent(n);
        if (is_valid_factor(n, d)) {
            fprintf(stderr, "cfx_rho_brent returned a nontrivial factor for prime "CFX_PRIuLIMB"\n",
                    n);
            assert(0);
        }
        return;
    }

    // Rho is randomized; fix seed for reproducibility
    srand(123456u);

    // Try a few times in case the internal random choices hit a bad cycle
    for (int attempts = 0; attempts < 5; ++attempts) {
        cfx_limb_t d = cfx_rho_brent(n);
        printf("cfx_rho_brent("CFX_PRIuLIMB") = "CFX_PRIuLIMB"\n", n, d);
        if (is_valid_factor(n, d)) {
            // ? verify cofactor primality or at least correctness
            cfx_limb_t m = n / d;
            assert(n == d * m);
            // Either side may still be composite; that's fine (we only test rho’s split)
            return;
        }
        srand(123456u + (unsigned)attempts + 1);
    }
    fprintf(stderr, "cfx_rho_brent failed to find a factor for "CFX_PRIuLIMB"\n",
            n);
    assert(0);
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
    cfx_limb_t n1 = 1;
    for (int i = 0; i < 10; ++i) n1 *= 3ULL;
    expect_factor(n1);

    cfx_limb_t n2 = 1ULL << 40;
    expect_factor(n2);
}

static void test_large_64bit_semiprime(void) {
    const cfx_limb_t p = 4294967291ULL; // 2^32 - 5, prime
    const cfx_limb_t q = 4294967279ULL; // 2^32 - 17, prime

    cfx_acc_t prod = ( cfx_acc_t)p * ( cfx_acc_t)q;
    cfx_limb_t n = (cfx_limb_t)prod; // still < 2^64

    expect_factor(n);
}

static void test_primes_do_not_yield_factors(void) {
    // These should NOT yield nontrivial factors
    assert(cfx_is_prime_u64(29));
    assert(cfx_is_prime_u64(97));
    assert(cfx_is_prime_u64(257));
    assert(cfx_is_prime_u64(65537));

    // cfx_rho_brent may return 0 or n — both are acceptable “no factor” signals.
    // We only assert it does NOT return a valid factor.
    cfx_limb_t primes[] = {29, 97, 257, 65537};
    for (size_t i = 0; i < sizeof(primes)/sizeof(primes[0]); ++i) {
        srand(42u);
        cfx_limb_t d = cfx_rho_brent(primes[i]);
        assert(!is_valid_factor(primes[i], d));
    }
}

int main(void) {
    test_small_semiprimes();
    test_carmichael_and_square_semiprimes();
    test_repeated_prime_powers();
    test_large_64bit_semiprime();
    test_primes_do_not_yield_factors();

    puts("cfx_rho_brent tests: OK");
    return 0;
}
