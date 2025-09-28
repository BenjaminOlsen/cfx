
#include "cfx/algo.h"
#include "cfx/primes.h"
#include "cfx/macros.h"

#include <stdio.h>

extern const size_t cfx_primes_len;
extern const uint32_t cfx_primes[];

static void test_prime_sieve(void) {
    uint32_t maxp = cfx_primes[cfx_primes_len-1];
    cfx_vec_t primes = cfx_sieve_primes(maxp);
    // CFX_PRINT_DBG("len: %zu, maxp: %u\n", cfx_primes_len, maxp);
    assert(primes.size == cfx_primes_len);
    for (size_t i = 0; i < primes.size; ++i) {
        CFX_ASSERT(primes.data[i] == cfx_primes[i]);
    }
    cfx_vec_free(&primes);
    CFX_PRINT_DBG("ok\n");
}

static void test_primality_test(void) {

    int fpos = 0;
    int fneg = 0;
    CFX_PRINT_DBG("cfx_primes_len: %zu\n", cfx_primes_len);
    for (size_t i = 0; i < cfx_primes_len - 1; ++i) {
        cfx_limb_t p1 = cfx_primes[i];
        cfx_limb_t p2 = cfx_primes[i+1];
        int isprime = cfx_is_prime_u64(p1);
        fpos += (1 - isprime);
        // CFX_PRINT_DBG(""CFX_PRIuLIMB" is prime: %d\n", p1, isprime);
        assert(isprime == 1);

        for (cfx_limb_t n = p1 + 1; n < p2; ++n) {
            int iscomp = !cfx_is_prime_u64(n);
            fneg += (1 - iscomp);
            // CFX_PRINT_DBG(""CFX_PRIuLIMB" is prime: %d\n", n, !iscomp);
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
    const cfx_limb_t p = UINT64_C(18446744073709551557); // known 64-bit prime
    cfx_limb_t a = UINT64_C(18446744073709551615) - 12345;
    cfx_limb_t b = UINT64_C(18446744073709551615) - 67890;
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

    const cfx_limb_t p = UINT64_C(18446744073709551557);
    cfx_limb_t a = UINT64_C(18446744073709551615) - 424242;
    cfx_limb_t e = UINT64_C(1234567890123456789);
    assert(cfx_powmod_u64(a,e,p) == UINT64_C(14322199550612880112));
}

static void test_fermat_primes(void) {
    const cfx_limb_t primes[] = {
        2, 3, 5, 17, 257, 65537, UINT64_C(18446744073709551557)
    };
    for (size_t i = 0; i < sizeof(primes)/sizeof(primes[0]); ++i) {
        cfx_limb_t p = primes[i];
        if (p == 2) {
            // pick a = 1
            assert(cfx_powmod_u64(1, p-1, p) == 1 % p);
        } else {
            // a not divisible by p
            cfx_limb_t a = 2;
            if (p % 2 == 0) a = 3;
            assert(cfx_powmod_u64(a, p-1, p) == 1);
        }
    }
}


/* --- simple reference multiply with carry propagation --- */
static void mul_ref(const cfx_limb_t* A, size_t na,
                    const cfx_limb_t* B, size_t nb,
                    cfx_limb_t* R) {
    size_t nout = na + nb;
    memset(R, 0, nout * sizeof(cfx_limb_t));

    for (size_t i = 0; i < na; ++i) {
        cfx_acc_t carry = 0;
        for (size_t j = 0; j < nb; ++j) {
            cfx_acc_t cur = (cfx_acc_t)A[i] * B[j]
                            + (cfx_acc_t)R[i + j]
                            + carry;
            R[i + j] = (cfx_limb_t)cur;
            carry    = cur >> CFX_LIMB_BITS;
        }
        size_t k = i + nb;
        while (carry) {
            cfx_acc_t t = (cfx_acc_t)R[k] + carry;
            R[k]   = (cfx_limb_t)t;
            carry  = t >> CFX_LIMB_BITS;
            ++k;
        }
    }
    // CFX_PRINT_DBG("A:\t"); 
    // PRINT_ARR(A, na);
    // CFX_PRINT_DBG("B:\t"); 
    // PRINT_ARR(B, nb);
    // CFX_PRINT_DBG("R:\t");
    // PRINT_ARR(R, nout);
}

/* --- helpers --- */
static void assert_eq_arr(const cfx_limb_t* X, const cfx_limb_t* Y, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        assert(X[i] == Y[i]);
    }
}

static void test_zero(void) {
    cfx_limb_t A[] = {0};
    cfx_limb_t B[] = {0, 0, 0};
    size_t na = 1, nb = 3, nout = na + nb;
    cfx_limb_t R1[4], R2[4];

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_one_by_small(void) {
    cfx_limb_t A[] = {1};
    cfx_limb_t B[] = {0x0123456789abcdefULL, 0xfedcba9876543210ULL};
    size_t na = 1, nb = 2, nout = na + nb;
    cfx_limb_t R1[3], R2[3];

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_single_limb_carry(void) {
    cfx_limb_t A[] = {0xffffffffffffffffULL};
    cfx_limb_t B[] = {0xffffffffffffffffULL};
    size_t na = 1, nb = 1, nout = na + nb;
    cfx_limb_t R1[2], R2[2];

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_mixed_lengths(void) {
    cfx_limb_t A[] = {0x0000000000000001ULL, 0x0000000000000000ULL, 0x0000000000000001ULL}; // 2^128 + 1
    cfx_limb_t B[] = {0x0000000000000000ULL, 0x0000000000000002ULL};                         // 2^65
    size_t na = 3, nb = 2, nout = na + nb;
    cfx_limb_t R1[5], R2[5];

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_all_ones_multi(void) {
    /* (2^192 - 1) * (2^128 - 1) */
    cfx_limb_t A[] = {0xffffffffffffffffULL, 0xffffffffffffffffULL, 0xffffffffffffffffULL};
    cfx_limb_t B[] = {0xffffffffffffffffULL, 0xffffffffffffffffULL};
    size_t na = 3, nb = 2, nout = na + nb;
    cfx_limb_t R1[5], R2[5];

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_power_of_two_alignment(void) {
    /* A = 2^64 + 1, B = 2^128 + 2^64 + 1 */
    cfx_limb_t A[] = {1ULL, 1ULL};
    cfx_limb_t B[] = {1ULL, 1ULL, 1ULL};
    size_t na = 2, nb = 3, nout = na + nb;
    cfx_limb_t R1[5], R2[5];

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

/* a tiny deterministic fuzz */
static cfx_limb_t xorshift64(cfx_limb_t* s) {
    cfx_limb_t x = *s;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    *s = x;
    return x;
}
static void test_fuzz_small(void) {
    cfx_limb_t seed = 0x123456789abcdef0ULL;
    for (int t = 0; t < 100; ++t) {
        size_t na = 1 + (xorshift64(&seed) % 4);
        size_t nb = 1 + (xorshift64(&seed) % 4);
        cfx_limb_t A[4], B[4], R1[8], R2[8];

        for (size_t i = 0; i < na; ++i) A[i] = xorshift64(&seed);
        for (size_t j = 0; j < nb; ++j) B[j] = xorshift64(&seed);

        cfx_mul_csa_portable(A, na, B, nb, R1);
        mul_ref(A, na, B, nb, R2);
        assert_eq_arr(R1, R2, na + nb);
    }
    PRINT_TEST(1);
}

int main(void) {
    test_prime_sieve();
    test_primality_test();
    test_mulmod_basic();
    test_powmod_edges();
    test_fermat_primes();
    // csa tests:
    test_zero();
    test_one_by_small();
    test_single_limb_carry();
    test_mixed_lengths();
    test_all_ones_multi();
    test_power_of_two_alignment();
    test_fuzz_small();
    puts("OK: cfx_mul_csa_portable tests passed.\n");
    return 0;
}

