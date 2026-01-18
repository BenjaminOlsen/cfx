/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/algo.h"
#include "cfx/primes.h"
#include "cfx/macros.h"
#include "cfx/mont.h"

#include <stdio.h>


static void test_prime_sieve(void) {
    uint32_t maxp = cfx_primes[cfx_primes_len-1];
    cfx_vec_t primes = cfx_sieve_primes(maxp);
    /* CFX_PRINT_DBG("len: %zu, maxp: %u\n", cfx_primes_len, maxp); */
    CFX_ASSERT(primes.size == cfx_primes_len);
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
        /* CFX_PRINT_DBG(""CFX_PRIuLIMB" is prime: %d\n", p1, isprime); */
        CFX_ASSERT(isprime == 1);

        for (cfx_limb_t n = p1 + 1; n < p2; ++n) {
            int iscomp = !cfx_is_prime_u64(n);
            fneg += (1 - iscomp);
            /* CFX_PRINT_DBG(""CFX_PRIuLIMB" is prime: %d\n", n, !iscomp); */
            CFX_ASSERT(iscomp == 1);
        }
    }
    CFX_PRINT_DBG("false pos: %d, false neg: %d\n", fpos, fneg);

    CFX_ASSERT(cfx_is_prime_u64(UINT64_MAX) == 0);
    /* largest prime < 2^64 */
    CFX_ASSERT(cfx_is_prime_u64(18446744073709551557ULL) == 1);
}

static void test_mulmod_basic(void) {
    CFX_ASSERT(cfx_mulmod_u64(0, 0, 7) == 0);
    CFX_ASSERT(cfx_mulmod_u64(1, 0, 7) == 0);
    CFX_ASSERT(cfx_mulmod_u64(3, 4, 7) == 5);     /* 12 % 7 = 5 */
    CFX_ASSERT(cfx_mulmod_u64(6, 6, 7) == 1);     /* 36 % 7 = 1 */

    /* m = 1 => always 0 */
    CFX_ASSERT(cfx_mulmod_u64(UINT64_C(123), UINT64_C(456), UINT64_C(1)) == 0);

    /* large operands, prime modulus near 2^64 */
    const uint64_t p = UINT64_C(18446744073709551557); /* known 64-bit prime */
    uint64_t a = UINT64_C(18446744073709551615) - 12345;
    uint64_t b = UINT64_C(18446744073709551615) - 67890;
    CFX_ASSERT(cfx_mulmod_u64(a,b,p) == 833451784);
}

static void test_powmod_edges(void) {
    CFX_ASSERT(cfx_powmod_u64(0, 0, 1) == 0);
    CFX_ASSERT(cfx_powmod_u64(0, 0, 2) == 1);
    CFX_ASSERT(cfx_powmod_u64(5, 0, 13) == 1);

    CFX_ASSERT(cfx_powmod_u64(0, 1, 13) == 0);
    CFX_ASSERT(cfx_powmod_u64(0, 5, 13) == 0);

    CFX_ASSERT(cfx_powmod_u64(2, 10, 1000) == 24); /* 1024 % 1000 */
    CFX_ASSERT(cfx_powmod_u64(3, 12345, 67891) == 31627);

    const uint64_t p = UINT64_C(18446744073709551557);
    uint64_t a = UINT64_C(18446744073709551615) - 424242;
    uint64_t e = UINT64_C(1234567890123456789);
    CFX_ASSERT(cfx_powmod_u64(a,e,p) == UINT64_C(14322199550612880112));
}

/* Tests for Montgomery 64-bit operations - covers fallback path on 32-bit platforms */
static void test_modinv64(void) {
    /* Small cases */
    CFX_ASSERT(cfx_modinv64(3, 7) == 5);     /* 3*5 = 15 ≡ 1 (mod 7) */
    CFX_ASSERT(cfx_modinv64(2, 5) == 3);     /* 2*3 = 6 ≡ 1 (mod 5) */
    CFX_ASSERT(cfx_modinv64(1, 7) == 1);     /* 1*1 = 1 */

    /* Edge case: no inverse exists */
    CFX_ASSERT(cfx_modinv64(2, 4) == 0);     /* gcd(2,4) = 2 != 1 */
    CFX_ASSERT(cfx_modinv64(0, 7) == 0);     /* 0 has no inverse */

    /* Verify via multiplication for a range of small primes */
    const uint64_t primes[] = {7, 11, 13, 17, 19, 23, 97, 101, 350381};
    for (size_t i = 0; i < sizeof(primes)/sizeof(primes[0]); i++) {
        uint64_t p = primes[i];
        for (uint64_t a = 1; a < 10 && a < p; a++) {
            uint64_t inv = cfx_modinv64(a, p);
            CFX_ASSERT(cfx_mulmod64(a, inv, p) == 1);
        }
    }

    /* Large prime near 2^64: p = 2^64 - 59, so r = 2^64 mod p = 59.
     * This is the critical case for Montgomery R^(-1) computation. */
    const uint64_t big_p = UINT64_C(18446744073709551557);
    uint64_t inv = cfx_modinv64(59, big_p);
    CFX_ASSERT(cfx_mulmod64(59, inv, big_p) == 1);

    /* Odd values and powers of 2 work correctly with the binary GCD */
    uint64_t test_vals[] = {2, 3, 4, 5, 7, 8, 9, 11, 13, 15, 16, 21, 25, 27, 32};
    for (size_t i = 0; i < sizeof(test_vals)/sizeof(test_vals[0]); i++) {
        inv = cfx_modinv64(test_vals[i], big_p);
        CFX_ASSERT(cfx_mulmod64(test_vals[i], inv, big_p) == 1);
    }
}

static void test_mont64_roundtrip(void) {
    /* Test Montgomery to/from conversion for various moduli */
    const uint64_t moduli[] = {
        7, 11, 17, 101, 350381,  /* small primes */
        UINT64_C(18446744073709551557)  /* largest 64-bit prime */
    };

    for (size_t i = 0; i < sizeof(moduli)/sizeof(moduli[0]); i++) {
        uint64_t n = moduli[i];
        cfx_mont64_t mont;
        cfx_mont64_init(&mont, n);

        /* Verify r * r_inv ≡ 1 (mod n) */
        CFX_ASSERT(cfx_mulmod64(mont.r, mont.r_inv, n) == 1);

        /* Test round-trip for several values */
        for (uint64_t a = 1; a < 10 && a < n; a++) {
            uint64_t aR = cfx_mont64_to(&mont, a);
            uint64_t back = cfx_mont64_from(&mont, aR);
            CFX_ASSERT(back == a);
        }
    }

    /* Edge case: values near n-1 for the large prime */
    const uint64_t big_p = UINT64_C(18446744073709551557);
    cfx_mont64_t mont;
    cfx_mont64_init(&mont, big_p);

    uint64_t test_vals[] = {1, 2, big_p - 1, big_p - 2, UINT64_C(12345678901234567)};
    for (size_t i = 0; i < sizeof(test_vals)/sizeof(test_vals[0]); i++) {
        uint64_t a = test_vals[i] % big_p;
        if (a == 0) continue;
        uint64_t aR = cfx_mont64_to(&mont, a);
        uint64_t back = cfx_mont64_from(&mont, aR);
        CFX_ASSERT(back == a);
    }
}

static void test_fermat_primes(void) {
    const uint64_t primes[] = {
        2, 3, 5, 17, 257, 65537, UINT64_C(18446744073709551557)
    };
    for (size_t i = 0; i < sizeof(primes)/sizeof(primes[0]); ++i) {
        uint64_t p = primes[i];
        if (p == 2) {
            /* pick a = 1 */
            CFX_ASSERT(cfx_powmod_u64(1, p-1, p) == 1 % p);
        } else {
            /* a not divisible by p */
            uint64_t a = 2;
            if (p % 2 == 0) a = 3;
            CFX_ASSERT(cfx_powmod_u64(a, p-1, p) == 1);
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
    /* CFX_PRINT_DBG("A:\t");  */
    /* PRINT_ARR(A, na); */
    /* CFX_PRINT_DBG("B:\t");  */
    /* PRINT_ARR(B, nb); */
    /* CFX_PRINT_DBG("R:\t"); */
    /* PRINT_ARR(R, nout); */
}

/* --- helpers --- */
static void assert_eq_arr(const cfx_limb_t* X, const cfx_limb_t* Y, size_t n) {
    for (size_t i = 0; i < n; ++i) {
        CFX_ASSERT(X[i] == Y[i]);
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
#if CFX_LIMB_BITS == 64
    cfx_limb_t B[] = {0x0123456789abcdefULL, 0xfedcba9876543210ULL};
    size_t na = 1, nb = 2, nout = na + nb;
    cfx_limb_t R1[3], R2[3];
#else
    cfx_limb_t B[] = {0x89abcdef, 0x01234567, 0x76543210, 0xfedcba98};
    size_t na = 1, nb = 4, nout = na + nb;
    cfx_limb_t R1[5], R2[5];
#endif

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_single_limb_carry(void) {
    cfx_limb_t A[] = {CFX_LIMB_MAX};
    cfx_limb_t B[] = {CFX_LIMB_MAX};
    size_t na = 1, nb = 1, nout = na + nb;
    cfx_limb_t R1[2], R2[2];

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_mixed_lengths(void) {
#if CFX_LIMB_BITS == 64
    cfx_limb_t A[] = {0x0000000000000001ULL, 0x0000000000000000ULL, 0x0000000000000001ULL}; /* 2^128 + 1 */
    cfx_limb_t B[] = {0x0000000000000000ULL, 0x0000000000000002ULL};                         /* 2^65 */
    size_t na = 3, nb = 2, nout = na + nb;
    cfx_limb_t R1[5], R2[5];
#else
    cfx_limb_t A[] = {1, 0, 0, 0, 1, 0}; /* 2^128 + 1 in 32-bit limbs */
    cfx_limb_t B[] = {0, 0, 2, 0};       /* 2^65 in 32-bit limbs */
    size_t na = 6, nb = 4, nout = na + nb;
    cfx_limb_t R1[10], R2[10];
#endif

    cfx_mul_csa_portable(A, na, B, nb, R1);
    mul_ref(A, na, B, nb, R2);
    assert_eq_arr(R1, R2, nout);
    PRINT_TEST(1);
}

static void test_all_ones_multi(void) {
#if CFX_LIMB_BITS == 64
    /* (2^192 - 1) * (2^128 - 1) */
    cfx_limb_t A[] = {CFX_LIMB_MAX, CFX_LIMB_MAX, CFX_LIMB_MAX};
    cfx_limb_t B[] = {CFX_LIMB_MAX, CFX_LIMB_MAX};
    size_t na = 3, nb = 2, nout = na + nb;
    cfx_limb_t R1[5], R2[5];
#else
    cfx_limb_t A[] = {CFX_LIMB_MAX, CFX_LIMB_MAX, CFX_LIMB_MAX, CFX_LIMB_MAX, CFX_LIMB_MAX, CFX_LIMB_MAX};
    cfx_limb_t B[] = {CFX_LIMB_MAX, CFX_LIMB_MAX, CFX_LIMB_MAX, CFX_LIMB_MAX};
    size_t na = 6, nb = 4, nout = na + nb;
    cfx_limb_t R1[10], R2[10];
#endif

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
static cfx_limb_t xorshift(cfx_limb_t* s) {
    cfx_limb_t x = *s;
#if CFX_LIMB_BITS == 64
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
#else
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
#endif
    *s = x;
    return x;
}
static void test_fuzz_small(void) {
    cfx_limb_t seed = (cfx_limb_t)0x12345678;
    for (int t = 0; t < 100; ++t) {
        size_t na = 1 + (xorshift(&seed) % 4);
        size_t nb = 1 + (xorshift(&seed) % 4);
        cfx_limb_t A[4], B[4], R1[8], R2[8];

        for (size_t i = 0; i < na; ++i) A[i] = xorshift(&seed);
        for (size_t j = 0; j < nb; ++j) B[j] = xorshift(&seed);

        cfx_mul_csa_portable(A, na, B, nb, R1);
        mul_ref(A, na, B, nb, R2);
        assert_eq_arr(R1, R2, na + nb);
    }
    PRINT_TEST(1);
}

static void test_cfz(void) {
    cfx_limb_t l0 = 0;
    CFX_ASSERT(cfx_clz(l0) == CFX_LIMB_BITS);
    cfx_limb_t l1 = 1;
    CFX_ASSERT(cfx_clz(l1) == CFX_LIMB_BITS-1);
    for (size_t k = 2; k < CFX_LIMB_BITS-1; ++k) {
        l1 <<= 1;
        CFX_ASSERT(cfx_clz(l1) == CFX_LIMB_BITS-k);
    }

    PRINT_TEST(1);
}

static void test_sqrt(void) {
    CFX_ASSERT (cfx_isqrt_nr(100) == 10);
    CFX_ASSERT (cfx_isqrt_nr(25) == 5);
    CFX_ASSERT (cfx_isqrt_nr(0) == 0);

#if CFX_LIMB_BITS == 64
    /* These tests require 64-bit limbs */
    CFX_ASSERT (cfx_isqrt_nr(0xFFFFFFFFFFFFFFULL) == 268435455);
    CFX_ASSERT (cfx_isqrt_nr(0xFFFFFFF0000000ULL) == 268435455);
    CFX_ASSERT (cfx_isqrt_nr(0xabcdefabadaULL) == 3436031);
#else
    /* 32-bit limb tests */
    CFX_ASSERT (cfx_isqrt_nr(0xFFFFFFFF) == 65535);
    CFX_ASSERT (cfx_isqrt_nr(0x10000) == 256);
#endif

    PRINT_TEST(1);
}


static void test_factor_u64(void) {
    cfx_vec_t primes, exps;

    /* should return -1 */
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, 0) == -1);

    /* should return 0 with empty vectors */
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, 1) == 0);
    CFX_ASSERT(primes.size == 0);
    CFX_ASSERT(exps.size == 0);
    cfx_vec_free(&primes);
    cfx_vec_free(&exps);

    /* test prime: 17 */
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, 17) == 0);
    CFX_ASSERT(primes.size == 1);
    CFX_ASSERT(primes.data[0] == 17);
    CFX_ASSERT(exps.data[0] == 1); 
    cfx_vec_free(&primes);
    cfx_vec_free(&exps);

    /* test power of 2: 64 = 2^6 */
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, 64) == 0);
    CFX_ASSERT(primes.size == 1);
    CFX_ASSERT(primes.data[0] == 2);
    CFX_ASSERT(exps.data[0] == 6);
    cfx_vec_free(&primes);
    cfx_vec_free(&exps);

    /* test composite: 60 = 2^2 * 3 * 5 */
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, 60) == 0);
    CFX_ASSERT(primes.size == 3);
    CFX_ASSERT(primes.data[0] == 2 && exps.data[0] == 2);
    CFX_ASSERT(primes.data[1] == 3 && exps.data[1] == 1);
    CFX_ASSERT(primes.data[2] == 5 && exps.data[2] == 1);
    cfx_vec_free(&primes);
    cfx_vec_free(&exps);

    /* test semiprime: 15 = 3 * 5 */
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, 15) == 0);
    CFX_ASSERT(primes.size == 2);
    CFX_ASSERT(primes.data[0] == 3 && exps.data[0] == 1);
    CFX_ASSERT(primes.data[1] == 5 && exps.data[1] == 1);
    cfx_vec_free(&primes);
    cfx_vec_free(&exps);

    /* test larger: 1000003 * 1000033 = 1000036000099 */
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, 1000036000099ULL) == 0);
    CFX_ASSERT(primes.size == 2);
    CFX_ASSERT(primes.data[0] == 1000003 && exps.data[0] == 1);
    CFX_ASSERT(primes.data[1] == 1000033 && exps.data[1] == 1);
    cfx_vec_free(&primes);
    cfx_vec_free(&exps);

    /* round-trip: product of factors equals original */
    uint64_t n = 123456789ULL;
    CFX_ASSERT(cfx_factor_u64(&primes, &exps, n) == 0);
    uint64_t product = 1;
    for (size_t i = 0; i < primes.size; i++) {
        for (uint64_t j = 0; j < exps.data[i]; j++) {
            product *= primes.data[i];
        }
    }
    CFX_ASSERT(product == n);
    cfx_vec_free(&primes);
    cfx_vec_free(&exps);
}

static void test_mul_csa_portable_fast(void) {
    cfx_mul_scratch_t scratch = {0};  /* must be zero-initialized before first alloc */
    cfx_limb_t A[4], B[4], R[8], R_ref[8];

    /* simple 1x1 multiplication */
    A[0] = 7;
    B[0] = 6;
    cfx_mul_scratch_alloc(&scratch, 2);
    cfx_mul_scratch_zero(&scratch, 2);
    cfx_mul_csa_portable_fast(A, 1, B, 1, R, &scratch);
    CFX_ASSERT(R[0] == 42);
    cfx_mul_scratch_free(&scratch);

    /* 2x2 multiplication with known result */
    /* (B + 1) * (B + 1) = B^2 + 2B + 1, where B = 2^LIMB_BITS */
    /* Result: R[0]=1, R[1]=2, R[2]=1, R[3]=0 */
    A[0] = 1; A[1] = 1;
    B[0] = 1; B[1] = 1;
    cfx_mul_scratch_alloc(&scratch, 4);
    cfx_mul_scratch_zero(&scratch, 4);
    cfx_mul_csa_portable_fast(A, 2, B, 2, R, &scratch);
    CFX_ASSERT(R[0] == 1);
    CFX_ASSERT(R[1] == 2);
    CFX_ASSERT(R[2] == 1);
    CFX_ASSERT(R[3] == 0);
    cfx_mul_scratch_free(&scratch);

    /* compare with cfx_mul_csa_portable (reference) using limb-safe values */
    A[0] = 0xDEADBEEF; A[1] = 0xCAFEBABE;
    B[0] = 0x12345678; B[1] = 0x9ABCDEF0;

    cfx_mul_csa_portable(A, 2, B, 2, R_ref);

    cfx_mul_scratch_alloc(&scratch, 4);
    cfx_mul_scratch_zero(&scratch, 4);
    cfx_mul_csa_portable_fast(A, 2, B, 2, R, &scratch);

    for (size_t i = 0; i < 4; i++) {
        CFX_ASSERT(R[i] == R_ref[i]);
    }
    cfx_mul_scratch_free(&scratch);

    /* asymmetric sizes with max limb values */
    A[0] = CFX_LIMB_MAX;
    B[0] = CFX_LIMB_MAX; B[1] = CFX_LIMB_MAX;

    cfx_mul_csa_portable(A, 1, B, 2, R_ref);

    cfx_mul_scratch_alloc(&scratch, 3);
    cfx_mul_scratch_zero(&scratch, 3);
    cfx_mul_csa_portable_fast(A, 1, B, 2, R, &scratch);

    for (size_t i = 0; i < 3; i++) {
        CFX_ASSERT(R[i] == R_ref[i]);
    }
    cfx_mul_scratch_free(&scratch);
}

int main(void) {
    CFX_TEST(test_prime_sieve);
    CFX_TEST(test_primality_test);
    CFX_TEST(test_mulmod_basic);
    CFX_TEST(test_powmod_edges);
    CFX_TEST(test_modinv64);
    CFX_TEST(test_mont64_roundtrip);
    CFX_TEST(test_fermat_primes);
    CFX_TEST(test_zero);
    CFX_TEST(test_one_by_small);
    CFX_TEST(test_single_limb_carry);
    CFX_TEST(test_mixed_lengths);
    CFX_TEST(test_all_ones_multi);
    CFX_TEST(test_power_of_two_alignment);
    CFX_TEST(test_fuzz_small);
    CFX_TEST(test_cfz);
    CFX_TEST(test_sqrt);
    CFX_TEST(test_factor_u64);
    CFX_TEST(test_mul_csa_portable_fast);
    puts("OK");
    return 0;
}

