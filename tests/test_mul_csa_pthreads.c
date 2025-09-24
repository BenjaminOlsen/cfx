// test_mul_csa_pthreads.c

#include "cfx/big.h"
#include "cfx/fmt.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/* --- helpers --- */
static void ensure_cap(cfx_big_t* b, size_t need) {
    // PRINT_BIG("ensure_cap BEFORE", b);
    cfx_big_reserve(b, need);
    // PRINT_BIG("ensure_cap AFTER", b);
}

static void fill_zero(cfx_big_t* b, size_t nlimbs) {
    ensure_cap(b, nlimbs);
    memset(b->limb, 0, nlimbs * sizeof(cfx_u64_t));
    b->n = nlimbs ? 1 : 0; // canonical zero usually n=0, but harmless if 1 with limb[0]==0
    if (nlimbs) b->limb[0] = 0;
    b->n = 0; // prefer canonical zero
}

static void fill_val(cfx_big_t* b, size_t nlimbs, cfx_u64_t v) {
    ensure_cap(b, nlimbs);
    memset(b->limb, 0, nlimbs * sizeof(cfx_u64_t));
    if (nlimbs == 0) { b->n = 0; return; }
    b->limb[0] = v;
    b->n = (v == 0) ? 0 : 1;
}

static void fill_ones(cfx_big_t* b, size_t nlimbs) {
    ensure_cap(b, nlimbs);
    for (size_t i = 0; i < nlimbs; ++i) b->limb[i] = UINT64_MAX;
    b->n = nlimbs;
}

static cfx_u64_t xorshift64(cfx_u64_t* s) {
    cfx_u64_t x = *s;
    x ^= x << 13; x ^= x >> 7; x ^= x << 17;
    *s = x;
    return x;
}

static void fill_rand(cfx_big_t* b, size_t nlimbs, cfx_u64_t seed) {
    ensure_cap(b, nlimbs);
    cfx_u64_t s = seed ? seed : 0x123456789abcdef0ULL;
    for (size_t i = 0; i < nlimbs; ++i) b->limb[i] = xorshift64(&s);
    if (nlimbs) b->limb[nlimbs-1] |= (1ULL << 63); // avoid leading zeros
    b->n = nlimbs;
}

static int big_equal(const cfx_big_t* a, const cfx_big_t* b, size_t upto /*nlimbs*/) {
    size_t n = upto;
    for (size_t i = 0; i < n; ++i) {
        cfx_u64_t av = (i < a->n) ? a->limb[i] : 0;
        cfx_u64_t bv = (i < b->n) ? b->limb[i] : 0;
        if (av != bv) return 0;
    }
    return 1;
}

static void print_limbs(const char* tag, const cfx_big_t* x, size_t upto) {
    fprintf(stderr, "%s: [", tag);
    for (size_t i = 0; i < upto; ++i) {
        cfx_u64_t v = (i < x->n) ? x->limb[i] : 0;
        fprintf(stderr, "%s"U64F"", (i ? ", " : ""), (unsigned long long)v);
    }
    fprintf(stderr, "]\n");
}


static void run_case(const char* name,
                     const cfx_big_t* b0_in, const cfx_big_t* m_in,
                     int threads)
{
    printf("\n------------------- %s -------------------\n", name);
    // Make working copies
    cfx_big_t b_ref;
    cfx_big_init(&b_ref);
    cfx_big_t b_pt;
    cfx_big_init(&b_pt);
    cfx_big_t m_copy;
    cfx_big_init(&m_copy);

    // Capacity for outputs: up to nout = nb + nm
    const size_t nout = b0_in->n + m_in->n;
    printf("nout = %zu\n", nout);
    printf("run_case: %s, b0_in->cap %zu, b0_in->n: %zu, m_in->cap %zu, m_in->n: %zu\n",
        name, b0_in->cap, b0_in->n, m_in->cap, m_in->n);
    PRINT_BIG("b0_in", b0_in);
    cfx_big_copy(&b_ref, b0_in);
    PRINT_BIG("b_ref", &b_ref);
    cfx_big_copy(&b_pt,  b0_in);
    PRINT_BIG("b_pt", &b_pt);
    cfx_big_copy(&m_copy, m_in); // to test immutability
    PRINT_BIG("m_copy", &m_copy);

    ensure_cap(&b_ref, nout);
    PRINT_BIG("b_ref after ensure", &b_ref);
    ensure_cap(&b_pt,  nout);
    PRINT_BIG("b_pt after ensure", &b_pt);

    // reference multiply (single-thread schoolbook)
    cfx_big_mul(&b_ref, m_in);

    // pthreads CSA multiply
    cfx_big_mul_csa_pthreads(&b_pt, m_in, threads);

    // Check equality (compare first nout limbs)
    if (!big_equal(&b_ref, &b_pt, nout)) {
        fprintf(stderr, "[FAIL] %s (threads=%d)\n", name, threads);
        print_limbs("b0  ", b0_in, b0_in->n);
        print_limbs("m   ", m_in,  m_in->n);
        print_limbs("ref ", &b_ref, nout);
        print_limbs("pth ", &b_pt,  nout);
        assert(0 && "mismatch between reference and pthreads CSA");
    }

    // Ensure m was not modified
    if (!big_equal(&m_copy, m_in, m_in->n)) {
        fprintf(stderr, "[FAIL] %s: multiplicand was modified!\n", name);
        assert(0 && "multiplicand modified");
    }

    cfx_big_free(&b_ref);
    cfx_big_free(&b_pt);
    cfx_big_free(&m_copy);
    printf(">>>>>>>>>>> [OK] %s\n", name);
}

/* --- Individual tests (happy paths) --- */

static void test_zero_cases(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);

    // 0 * X
    fill_zero(&b, 0);
    fill_rand(&m, 3, 0xA1);
    run_case("zero_times_rand", &b, &m, 1);
    run_case("zero_times_rand_t4", &b, &m, 4);

    // X * 0
    fill_rand(&b, 4, 0xB2);
    fill_zero(&m, 0);
    run_case("rand_times_zero", &b, &m, 1);
    run_case("rand_times_zero_t3", &b, &m, 3);

    cfx_big_free(&b);
    cfx_big_free(&m);
}

static void test_one_cases(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);

    // 1 * X
    fill_val(&b, 1, 1);
    fill_rand(&m, 5, 0xC3);
    run_case("one_times_rand", &b, &m, 1);
    run_case("one_times_rand_t6", &b, &m, 6);

    // X * 1
    fill_rand(&b, 6, 0xD4);
    fill_val(&m, 1, 1);
    run_case("rand_times_one", &b, &m, 1);
    run_case("rand_times_one_t2", &b, &m, 2);

    cfx_big_free(&b);
    cfx_big_free(&m);
}

static void test_all_ones_small(void) {
    cfx_big_t b, m; 
    cfx_big_init(&b);
    cfx_big_init(&m);

    fill_ones(&b, 3); // (2^192 - 1)
    fill_ones(&m, 2); // (2^128 - 1)
    run_case("all_ones_3x2_t1", &b, &m, 1);
    run_case("all_ones_3x2_t4", &b, &m, 4);

    cfx_big_free(&b);
    cfx_big_free(&m);
}

static void test_powers_of_two_alignment(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);

    // b = 2^128 + 1, m = 2^65
    ensure_cap(&b, 3); memset(b.limb, 0, 3*sizeof(cfx_u64_t)); b.limb[0]=1; b.limb[2]=1; b.n=3;
    ensure_cap(&m, 2); memset(m.limb, 0, 2*sizeof(cfx_u64_t)); m.limb[1]=2; m.n=2;

    run_case("power2_align_t1", &b, &m, 1);
    run_case("power2_align_t8", &b, &m, 8);

    cfx_big_free(&b);
    cfx_big_free(&m);
}

static void test_asymmetric_sizes(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);

    fill_rand(&b, 7, 0xEE);
    fill_rand(&m, 2, 0xAD);
    run_case("asym_7x2_t1", &b, &m, 1);
    run_case("asym_7x2_t3", &b, &m, 3);

    fill_rand(&b, 2, 0x12);
    fill_rand(&m, 9, 0x34);
    run_case("asym_2x9_t1", &b, &m, 1);
    run_case("asym_2x9_t4", &b, &m, 4);

    cfx_big_free(&b);
    cfx_big_free(&m);
}

static void test_small_fuzz(void) {
    const int threads_vec[] = {1,2,3,4};
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);
    cfx_u64_t seed = 0xCAFEBABE12345678ULL;

    for (int t = 0; t < 20; ++t) {
        size_t nb = 1 + (xorshift64(&seed) % 6);   // 1..6 limbs
        size_t nm = 1 + (xorshift64(&seed) % 6);
        fill_rand(&b, nb, xorshift64(&seed));
        fill_rand(&m, nm, xorshift64(&seed));
        for (size_t i = 0; i < sizeof(threads_vec)/sizeof(threads_vec[0]); ++i) {
            int th = threads_vec[i];
            run_case("fuzz_small", &b, &m, th);
        }
    }

    cfx_big_free(&b);
    cfx_big_free(&m);
}

/* a medium case to shake out nout handling */
static void test_medium_case_128x128(void) {
    cfx_big_t b, m;
    cfx_big_init(&b);
    cfx_big_init(&m);
    fill_rand(&b, 128, 0x1111);
    fill_rand(&m, 128, 0x2222);
    run_case("medium_128x128_t1",  &b, &m, 1);
    run_case("medium_128x128_t4",  &b, &m, 4);
    run_case("medium_128x128_t12", &b, &m, 12);
    cfx_big_free(&b); 
    cfx_big_free(&m);
}

/* --- main --- */
int main(void) {
    test_zero_cases();
    test_one_cases();
    test_all_ones_small();
    test_powers_of_two_alignment();
    test_asymmetric_sizes();
    test_small_fuzz();
    test_medium_case_128x128();
    puts("OK: cfx_big_mul_csa_pthreads happy-path tests passed.");
    return 0;
}
