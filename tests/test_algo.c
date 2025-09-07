#include "cfx/algo.h"
#include "cfx/primes.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <assert.h>

extern const size_t cfx_primes_len;
extern const uint32_t cfx_primes[];

static void test_prime_sieve(void) {
    uint32_t maxp = cfx_primes[cfx_primes_len-1];
    cfx_vec_t primes = cfx_sieve_primes(maxp);
    // CFX_PRINT_DBG("len: %zu, maxp: %u\n", cfx_primes_len, maxp);
    assert(primes.size == cfx_primes_len);
    for (size_t i = 0; i < primes.size; ++i) {
        assert(primes.data[i] == cfx_primes[i]);
    }
    cfx_vec_free(&primes);
    CFX_PRINT_DBG("ok\n");
}

static void test_primality_test(void) {

    int fpos = 0;
    int fneg = 0;
    CFX_PRINT_DBG("cfx_primes_len: %zu\n", cfx_primes_len);
    for (size_t i = 0; i < cfx_primes_len - 1; ++i) {
        uint64_t p1 = cfx_primes[i];
        uint64_t p2 = cfx_primes[i+1];
        int isprime = cfx_is_prime_u64(p1);
        fpos += (1 - isprime);
        // CFX_PRINT_DBG("%llu is prime: %d\n", p1, isprime);
        assert(isprime == 1);

        for (uint64_t n = p1 + 1; n < p2; ++n) {
            int iscomp = !cfx_is_prime_u64(n);
            fneg += (1 - iscomp);
            // CFX_PRINT_DBG("%llu is prime: %d\n", n, !iscomp);
            assert(iscomp == 1);
        }
    }
    CFX_PRINT_DBG("false pos: %d, false neg: %d\n", fpos, fneg);

    assert(cfx_is_prime_u64(UINT64_MAX) == 0);
    /* largest prime < 2^64 */
    assert(cfx_is_prime_u64(18446744073709551557ULL) == 1);
}

static void test_mulmod_basic(void) {
    assert(cfx_mulmod_u64(0, 0, 7) == 0);
    assert(cfx_mulmod_u64(1, 0, 7) == 0);
    assert(cfx_mulmod_u64(3, 4, 7) == 5);     // 12 % 7 = 5
    assert(cfx_mulmod_u64(6, 6, 7) == 1);     // 36 % 7 = 1

    // m = 1 => always 0
    assert(cfx_mulmod_u64(UINT64_C(123), UINT64_C(456), UINT64_C(1)) == 0);

    // large operands, prime modulus near 2^64
    const uint64_t p = UINT64_C(18446744073709551557); // known 64-bit prime
    uint64_t a = UINT64_C(18446744073709551615) - 12345;
    uint64_t b = UINT64_C(18446744073709551615) - 67890;
    assert(cfx_mulmod_u64(a,b,p) == 833451784);
}

static void test_powmod_edges(void) {
    assert(cfx_powmod_u64(0, 0, 1) == 0);
    assert(cfx_powmod_u64(0, 0, 2) == 1);
    assert(cfx_powmod_u64(5, 0, 13) == 1);

    assert(cfx_powmod_u64(0, 1, 13) == 0);
    assert(cfx_powmod_u64(0, 5, 13) == 0);

    assert(cfx_powmod_u64(2, 10, 1000) == 24); // 1024 % 1000
    assert(cfx_powmod_u64(3, 12345, 67891) == 31627);

    const uint64_t p = UINT64_C(18446744073709551557);
    uint64_t a = UINT64_C(18446744073709551615) - 424242;
    uint64_t e = UINT64_C(1234567890123456789);
    assert(cfx_powmod_u64(a,e,p) == UINT64_C(14322199550612880112));
}

static void test_fermat_primes(void) {
    const uint64_t primes[] = {
        2, 3, 5, 17, 257, 65537, UINT64_C(18446744073709551557)
    };
    for (size_t i = 0; i < sizeof(primes)/sizeof(primes[0]); ++i) {
        uint64_t p = primes[i];
        if (p == 2) {
            // pick a = 1
            assert(cfx_powmod_u64(1, p-1, p) == 1 % p);
        } else {
            // a not divisible by p
            uint64_t a = 2;
            if (p % 2 == 0) a = 3;
            assert(cfx_powmod_u64(a, p-1, p) == 1);
        }
    }
}

int main(void) {
    test_prime_sieve();
    test_primality_test();
    test_mulmod_basic();
    test_powmod_edges();
    test_fermat_primes();
    return 0;
}

