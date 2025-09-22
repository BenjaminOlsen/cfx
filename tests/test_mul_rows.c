// test_mul_rows.c

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "cfx/big.h" 
#include "cfx/macros.h"


// ---- Tiny helpers -----------------------------------------------------------

static void big_set_limbs(cfx_big_t* x, const uint64_t* limbs, size_t n)
{
    cfx_big_reserve(x, n);
    if (n) memcpy(x->limb, limbs, n * sizeof(uint64_t));
    x->n = n;
    // Normalize: trim leading zeros in case caller passed any
    while (x->n && x->limb[x->n - 1] == 0) x->n--;
}

static int big_equal(const cfx_big_t* a, const cfx_big_t* b)
{
    if (a->n != b->n) return 0;
    for (size_t i = 0; i < a->n; ++i)
        if (a->limb[i] != b->limb[i]) return 0;
    return 1;
}

static void big_print(const char* tag, const cfx_big_t* x)
{
    printf("%s n=%zu [", tag, x->n);
    for (size_t i = x->n; i-- > 0;) {
        printf("%016" PRIx64 "", (unsigned long long)x->limb[i]);
        if (i) printf("_");
    }
    printf("]\n");
}

// Reference schoolbook multiply: out = a * b  (no threading, exact)
static void big_mul_ref(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b)
{
    cfx_big_reserve(out, a->n + b->n + 1);
    memset(out->limb, 0, (a->n + b->n + 1) * sizeof(uint64_t));
    out->n = a->n + b->n + 1;

#if defined(__SIZEOF_INT128__)
    for (size_t j = 0; j < b->n; ++j) {
        __uint128_t carry = 0;
        for (size_t i = 0; i < a->n; ++i) {
            __uint128_t t = (__uint128_t)a->limb[i] * b->limb[j]
                          + out->limb[i + j] + carry;
            out->limb[i + j] = (uint64_t)t;
            carry = t >> 64;
        }
        size_t k = a->n + j;
        while (carry) {
            __uint128_t t = (__uint128_t)out->limb[k] + carry;
            out->limb[k] = (uint64_t)t;
            carry = t >> 64;
            ++k;
        }
    }
#else
#  error "This reference multiply needs __int128. For MSVC, port to _umul128/_addcarry_u64."
#endif

    // normalize
    while (out->n && out->limb[out->n - 1] == 0) out->n--;
}

// Deterministic PRNG (SplitMix64)
static uint64_t splitmix64(uint64_t* s) {
    uint64_t z = (*s += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}

static void big_rand(cfx_big_t* x, size_t n, uint64_t* seed)
{
    cfx_big_reserve(x, n);
    for (size_t i = 0; i < n; ++i) x->limb[i] = splitmix64(seed);
    x->n = n;
    // Ensure top limb is nonzero (if n>0) to avoid degenerate trims
    if (n) { if (x->limb[n - 1] == 0) x->limb[n - 1] = 1; }
}

// ---- Tests ------------------------------------------------------------------

static void test_mul_zero_zero(void)
{
    cfx_big_t a, b;
    cfx_big_init(&a); cfx_big_init(&b);
    cfx_big_from_u64(&a, 0);
    cfx_big_from_u64(&b, 0);

    cfx_big_mul_rows_pthreads(&a, &b, 1);
    assert(a.n == 0 || (a.n == 1 && a.limb[0] == 0));
    cfx_big_free(&a); cfx_big_free(&b);
}

static void test_mul_zero_x(void)
{
    cfx_big_t a, b;
    cfx_big_init(&a); cfx_big_init(&b);
    cfx_big_from_u64(&a, 0);
    uint64_t one = 123456789ULL;
    big_set_limbs(&b, &one, 1);

    cfx_big_mul_rows_pthreads(&a, &b, 4);
    assert(a.n == 0 || (a.n == 1 && a.limb[0] == 0));
    cfx_big_free(&a); cfx_big_free(&b);
}

static void test_mul_x_zero(void)
{
    cfx_big_t a, b;
    cfx_big_init(&a); cfx_big_init(&b);
    uint64_t one = 987654321ULL;
    big_set_limbs(&a, &one, 1);
    cfx_big_from_u64(&b, 0);

    cfx_big_mul_rows_pthreads(&a, &b, 2);
    assert(a.n == 0 || (a.n == 1 && a.limb[0] == 0));
    cfx_big_free(&a); cfx_big_free(&b);
}

static void test_mul_by_one(void)
{
    cfx_big_t a, m, ref;
    cfx_big_init(&a); cfx_big_init(&m); cfx_big_init(&ref);

    uint64_t limbs[] = {0x0123456789abcdefULL, 0xfedcba9876543210ULL};
    big_set_limbs(&a, limbs, 2);
    cfx_big_from_u64(&m, 1);

    // reference: a * 1 = a
    big_set_limbs(&ref, limbs, 2);

    cfx_big_mul_rows_pthreads(&a, &m, 8);
    assert(big_equal(&a, &ref));

    cfx_big_free(&a); cfx_big_free(&m); cfx_big_free(&ref);
}

static void test_cross_limb_carry_small(void)
{
    // (2^64 - 1)^2 = high: 0xfffffffffffffffe, low: 0x0000000000000001
    cfx_big_t a, m, ref;
    cfx_big_init(&a); cfx_big_init(&m); cfx_big_init(&ref);

    uint64_t max = ~0ULL;
    big_set_limbs(&a, &max, 1);
    big_set_limbs(&m, &max, 1);

    uint64_t expect[] = { 0x0000000000000001ULL, 0xfffffffffffffffeULL };
    big_set_limbs(&ref, expect, 2);

    cfx_big_mul_rows_pthreads(&a, &m, 1);
    assert(big_equal(&a, &ref));

    cfx_big_free(&a); cfx_big_free(&m); cfx_big_free(&ref);
}

static void test_small_vector_known(void)
{
    // [1,1] * [1,1] = [1,2,1]  (i.e., (1 + B)^2)
    cfx_big_t a, m, ref;
    cfx_big_init(&a); cfx_big_init(&m); cfx_big_init(&ref);

    uint64_t v[] = {1ULL, 1ULL};
    big_set_limbs(&a, v, 2);
    big_set_limbs(&m, v, 2);

    uint64_t expect[] = {1ULL, 2ULL, 1ULL};
    big_set_limbs(&ref, expect, 3);

    cfx_big_mul_rows_pthreads(&a, &m, 3);
    assert(big_equal(&a, &ref));

    cfx_big_free(&a); cfx_big_free(&m); cfx_big_free(&ref);
}

static void test_random_compare_ref(size_t na, size_t nb, int threads, uint64_t seed_init, char* msg)
{
    cfx_big_t a, b, ref, tmpa;
    cfx_big_init(&a); cfx_big_init(&b); cfx_big_init(&ref); cfx_big_init(&tmpa);

    uint64_t seed = seed_init;
    big_rand(&a, na, &seed);
    big_rand(&b, nb, &seed);

    // ref = a * b  (commutative)
    big_mul_ref(&ref, &a, &b);

    // tmpa = a; tmpa *= b via threaded mul
    // note: API under test overwrites first arg
    big_set_limbs(&tmpa, a.limb, a.n);
    cfx_big_mul_rows_pthreads(&tmpa, &b, threads);

    if (!big_equal(&tmpa, &ref)) {
        big_print("A", &a);
        big_print("B", &b);
        big_print("REF", &ref);
        big_print("GOT", &tmpa);
    }
    assert(big_equal(&tmpa, &ref));
    if (msg) printf("%s\n", msg);
    cfx_big_free(&a); cfx_big_free(&b); cfx_big_free(&ref); cfx_big_free(&tmpa);
}

static void test_thread_counts_agree(void)
{
    cfx_big_t a, b, t0, t1, t8, t32;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&t0);
    cfx_big_init(&t1);
    cfx_big_init(&t8);
    cfx_big_init(&t32);

    uint64_t seed = 12345;
    big_rand(&a, 764, &seed);
    big_rand(&b, 857, &seed);

    /* normal mul sanity check */
    cfx_big_copy(&t0, &a);
    cfx_big_mul(&t0, &b);

    big_set_limbs(&t1, a.limb, a.n);
    cfx_big_mul_rows_pthreads(&t1, &b, 1);

    big_set_limbs(&t8, a.limb, a.n);
    cfx_big_mul_rows_pthreads(&t8, &b, 8);

    big_set_limbs(&t32, a.limb, a.n);
    cfx_big_mul_rows_pthreads(&t32, &b, 32);

    assert(big_equal(&t0, &t1));
    assert(big_equal(&t1, &t8));
    assert(big_equal(&t32, &t8));

    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_free(&t0);
    cfx_big_free(&t1);
    cfx_big_free(&t8);
    cfx_big_free(&t32);
}

// ---- Main -------------------------------------------------------------------

#define TEST_RAND(na, nb, threads, seed_init) \
    test_random_compare_ref(na, nb, threads, seed_init, \
        "test_random_compare_ref(" STR(na) ", " STR(nb) ", " STR(threads) ", " STR(seed_init) ") - OK")

int main(void)
{
    CFX_TEST(test_mul_zero_zero);
    CFX_TEST(test_mul_zero_x);
    CFX_TEST(test_mul_x_zero);
    CFX_TEST(test_mul_by_one);
    CFX_TEST(test_cross_limb_carry_small);
    CFX_TEST(test_small_vector_known);
    CFX_TEST(test_thread_counts_agree);

    // Randoms vs reference (sizes & threads)
    TEST_RAND(1, 1, 1, 0xC0FFEEF);
    TEST_RAND(3, 2, 2, 0xDEADBAEF);
    TEST_RAND(8, 9, 4, 0x12345678);
    TEST_RAND(32, 17, 8, 0xB1DC0DE);
    TEST_RAND(64, 64, 4, 0xFEEDFACE);
    TEST_RAND(128, 64, 4, 0xFEADFACA);
    TEST_RAND(4000, 6410, 4, 0xAC0FFEE);
    TEST_RAND(6430, 6420, 8, 0xFAEDFFCE);
    TEST_RAND(14302, 14203, 8, 0xFE000CE);
    TEST_RAND(14302, 14203, 18, 0xFABADA);
    printf("OK: all tests passed.\n");
    return 0;
}
