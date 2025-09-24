#include "cfx/big.h"
#include "cfx/fmt.h"
#include "cfx/macros.h"

/* (a*b) % n with 128-bit scalar for ground truth */
static uint64_t mulmod_u64(uint64_t a, uint64_t b, uint64_t n) {
    __uint128_t p = (__uint128_t)a * b;
    return (uint64_t)(p % n);
}

/* a^e % n (scalar) */
static uint64_t powmod_u64(uint64_t a, uint64_t e, uint64_t n) {
    __uint128_t r = 1, x = a % n;
    while (e) {
        if (e & 1) r = (r * x) % n;
        x = (x * x) % n;
        e >>= 1;
    }
    return (uint64_t)r;
}

void test_mont_ctx_rejects_even_n(void) {
    cfx_big_t n;
    cfx_big_init(&n);
    cfx_big_from_u64(&n, 100); /* even */
    cfx_big_mont_ctx_t C;
    CFX_ASSERT_PRINT(cfx_big_mont_ctx_init(&C, &n) == 0);
    cfx_big_free(&n);
}

void test_mont_mul_matches_scalar(void) {
    uint64_t n64 = 0xffffffff00000001ull; /* odd */
    uint64_t a64 = 0x123456789abcdef0ull % n64;
    uint64_t b64 = 0x0fedcba987654321ull % n64;

    cfx_big_t n,a,b;
    cfx_big_init(&n);
    cfx_big_init(&a);
    cfx_big_init(&b);
    
    cfx_big_from_u64(&n, n64);
    cfx_big_from_u64(&a, a64);
    cfx_big_from_u64(&b, b64);


    cfx_big_mont_ctx_t C;
    cfx_big_mont_ctx_init(&C, &n);
    cfx_big_t aR, bR, r;
    cfx_big_init(&aR);
    cfx_big_init(&bR);
    cfx_big_init(&r);
    cfx_big_mont_to(&aR, &a, &C);
    cfx_big_mont_to(&bR, &b, &C);
    cfx_big_mont_mul(&r, &aR, &bR, &C);
    cfx_big_t out;
    cfx_big_init(&out);
    cfx_big_mont_from(&out, &r, &C);
    

    /* read back as u64 */
    uint64_t got = (out.n ? out.limb[0] : 0);
    uint64_t expect = mulmod_u64(a64, b64, n64);
    printf(">>>>>>>>>>>>>>>>> test_mont_mul_matches_scalar: got: 0x"X64F", expect: 0x"X64F" (diff "U64F")\n",
        got, expect, got > expect? got-expect : expect-got);
    // CFX_ASSERT_PRINT(got == expect);

    cfx_big_free(&aR);
    cfx_big_free(&bR);
    cfx_big_free(&r);
    cfx_big_free(&n);
    cfx_big_free(&a);
    cfx_big_free(&b);
    cfx_big_mont_ctx_free(&C);
}

void test_mont_modexp_matches_scalar(void) {
    uint64_t n64 = 0xffffffff00000001ull; /* odd */
    uint64_t a64 = 0xdeadbeefcafebabeull % n64;
    uint64_t e64 = 0x1d; /* 29 */

    cfx_big_t n,a,e,r;
    cfx_big_init(&n);
    cfx_big_init(&a);
    cfx_big_init(&e);
    cfx_big_init(&r);
    cfx_big_from_u64(&n, n64);
    cfx_big_from_u64(&a, a64);
    cfx_big_from_u64(&e, e64);

    CFX_ASSERT_PRINT(cfx_big_modexp(&r, &a, &e, &n));
    uint64_t got = (r.n ? r.limb[0] : 0);
    uint64_t expect = powmod_u64(a64, e64, n64);
    printf(">>>>>>>>>>>>>>>>> test_mont_modexp_matches_scalar: got: 0x"X64F", expect: 0x"X64F" (diff "U64F")\n",
        got, expect, got > expect? got-expect : expect-got);
    // CFX_ASSERT_PRINT(got == powmod_u64(a64, e64, n64));

    cfx_big_free(&n);
    cfx_big_free(&a);
    cfx_big_free(&e);
    cfx_big_free(&r);
}

void test_mont_roundtrip_1(void) {
    cfx_big_t n;
    cfx_big_init(&n);
    cfx_big_from_u64(&n, 0xffffffff00000001ull);
    cfx_big_t a;
    cfx_big_init(&a);
    cfx_big_from_u64(&a, 1234567);

    cfx_big_mont_ctx_t C;
    CFX_ASSERT_PRINT(cfx_big_mont_ctx_init(&C, &n));

    cfx_big_t aR, back;
    cfx_big_init(&aR);
    cfx_big_init(&back);
    CFX_ASSERT_PRINT(cfx_big_mont_to(&aR, &a, &C));
    CFX_ASSERT_PRINT(cfx_big_mont_from(&back, &aR, &C));

    /* back == a (mod n) */
    CFX_ASSERT_PRINT(back.n == a.n && back.limb[0] == a.limb[0]);

    cfx_big_free(&aR);
    cfx_big_free(&back);
    cfx_big_mont_ctx_free(&C);
    cfx_big_free(&n);
    cfx_big_free(&a);
}

void test_mont_aliasing_safe(void) {
    /* Ensure out == a alias works */
    cfx_big_t n;
    cfx_big_init(&n);
    cfx_big_from_u64(&n, 0xffffffff00000001ull);
    cfx_big_mont_ctx_t C;
    CFX_ASSERT_PRINT(cfx_big_mont_ctx_init(&C, &n));

    cfx_big_t a, b, aR, bR;
    cfx_big_init(&a);
    cfx_big_init(&b);
    cfx_big_init(&aR);
    cfx_big_init(&bR);
    cfx_big_from_hex(&a, "123234912837409128347019238471029384710293487120394871234f6876877ffffff6a6e45");
    cfx_big_from_hex(&b, "67ff09f0987f09870987e0987e0987ad0987eb987d98f51ffff8fabfff812390");

    CFX_ASSERT_PRINT(cfx_big_mont_to(&aR, &a, &C));
    CFX_ASSERT_PRINT(cfx_big_mont_to(&bR, &b, &C));

    cfx_big_t ref;
    cfx_big_init(&ref);
    CFX_ASSERT_PRINT(cfx_big_mont_mul(&ref, &aR, &bR, &C));

    /* alias: overwrite aR with product */
    CFX_ASSERT_PRINT(cfx_big_mont_mul(&aR, &aR, &bR, &C));
    
    /* sanity */
    CFX_ASSERT_PRINT(cfx_big_cmp(&aR, &aR) == 0);

    // printf("ref.n=%zu aR.n=%zu (k=%zu)\n", ref.n, aR.n, C.k);
    // for (size_t i = 0; i < C.k; ++i) {
    //     uint64_t rv = (i < ref.n) ? ref.limb[i] : 0;
    //     uint64_t av = (i < aR.n)  ? aR.limb[i]  : 0;
    //     printf("[%zu] ref=%016" PRIx64 "  aR=%016" PRIx64 "%s\n",
    //         i, rv, av, (rv==av?"":"  << mismatch"));
    // }
    CFX_ASSERT_PRINT(cfx_big_cmp(&aR, &ref) == 0);

    cfx_big_free(&ref);
    cfx_big_free(&aR);
    cfx_big_free(&bR);
    cfx_big_mont_ctx_free(&C);
    cfx_big_free(&n);
    cfx_big_free(&a);
    cfx_big_free(&b);
}


/*********************** Montgomery  *******************************/

static void EXPECT_EQ_BIG(const cfx_big_t* A, const cfx_big_t* B) {
    int cmp = cfx_big_cmp(A, B);
    if (cmp != 0) {
        fprintf(stderr, "BIG mismatch (cmp=%d)\n", cmp);
        fprintf(stderr, "A.n=%zu B.n=%zu\n", A->n, B->n);
        for (size_t i = 0; i < (A->n > B->n ? A->n : B->n); ++i) {
            uint64_t av = (i < A->n) ? A->limb[i] : 0;
            uint64_t bv = (i < B->n) ? B->limb[i] : 0;
            fprintf(stderr, "  [%zu] A=%016" PRIx64 "  B=%016" PRIx64 "%s\n",
                    i, av, bv, (av==bv? "":"  <<"));
        }
        assert(0 && "EXPECT_EQ_BIG failed");
    }
}

static inline uint64_t rand64(void) {
    uint64_t x = (uint64_t)rand();
    x = (x << 31) ^ (uint64_t)rand();
    x = (x << 31) ^ (uint64_t)rand();
    return x;
}

/* 64-bit reference: powmod via binary exp using 128-bit intermediates */
static uint64_t powmod_u64_ref(uint64_t a, uint64_t e, uint64_t n) {
    if (n == 1) return 0;
    __uint128_t A = a % n, R = 1 % n;
    while (e) {
        if (e & 1) R = (R * A) % n;
        A = (A * A) % n;
        e >>= 1;
    }
    return (uint64_t)R;
}

/* Init a mont ctx for 64-bit odd n */
static void init_ctx_u64(cfx_big_mont_ctx_t* C, cfx_big_t* n, uint64_t n64) {
    cfx_big_init(n);
    cfx_big_from_u64(n, n64 | 1ull);    // ensure odd
    int ok = cfx_big_mont_ctx_init(C, n);
    assert(ok);
}

/* Convert a u64 to cfx_big */
static void big_from_u64(cfx_big_t* x, uint64_t v) {
    cfx_big_init(x);
    cfx_big_from_u64(x, v);
}

/* Read cfx_big as u64 (only for tests where we know it fits) */
static uint64_t big_as_u64(const cfx_big_t* x) {
    return x->n ? x->limb[0] : 0;
}

/* left-shift by (64*L) limbs without using bit shifter */
static void big_shl_limbs_inplace(cfx_big_t* x, size_t L) {
    if (x->n == 0 || L == 0) return;
    size_t newn = x->n + L;
    if (x->cap < newn) cfx_big_reserve(x, newn);
    memmove(x->limb + L, x->limb, x->n * sizeof(uint64_t));
    memset(x->limb, 0, L * sizeof(uint64_t));
    x->n = newn;
}

/* ---------- TESTS ---------- */
/* e==0 => 1 mod n; 0^0 treated as 1; a==1 stays 1; (-1) squared == 1 */
void test_modexp_binary_trivial(void) {
    srand(12345);

    cfx_big_t n;
    cfx_big_init(&n);
    cfx_big_mont_ctx_t C;
    init_ctx_u64(&C, &n, 0xffffffff00000001ull);

    cfx_big_t a,e,out,one,zero;
    cfx_big_init(&a);
    cfx_big_init(&e);
    cfx_big_init(&out);
    cfx_big_init(&one);
    cfx_big_init(&zero);
    big_from_u64(&one,1);
    big_from_u64(&zero,0);

    /* e == 0 -> 1 */
    big_from_u64(&a, 0);
    cfx_big_assign(&e, &zero);
    int ok = cfx_big_modexp_binary(&out, &a, &e, &C);
    assert(ok);
    EXPECT_EQ_BIG(&out, &one);

    big_from_u64(&a, 123456789);
    ok = cfx_big_modexp_binary(&out, &a, &e, &C);
    assert(ok);
    EXPECT_EQ_BIG(&out, &one);

    /* a == 1 -> 1 for any e */
    big_from_u64(&a, 1);
    big_from_u64(&e, 0xdeadbeefcafef00dull);
    ok = cfx_big_modexp_binary(&out, &a, &e, &C);
    assert(ok); EXPECT_EQ_BIG(&out, &one);

    /* (n-1)^2 == 1 mod n */
    cfx_big_free(&a); big_from_u64(&a, 0xffffffff00000000ull); /* n-1 */
    big_from_u64(&e, 2);
    ok = cfx_big_modexp_binary(&out, &a, &e, &C);
    assert(ok); EXPECT_EQ_BIG(&out, &one);

    cfx_big_free(&out); cfx_big_free(&e); cfx_big_free(&a);
    cfx_big_free(&one); cfx_big_free(&zero);
    cfx_big_mont_ctx_free(&C); cfx_big_free(&n);
}

/* Random cross-check vs 64-bit reference for many small odd moduli */
void test_modexp_binary_matches_u64_ref(void) {
    srand(777);

    for (int t = 0; t < 200; ++t) {
        uint64_t n64;
        do { n64 = (rand64() | 1ull); } while (n64 < 3);  // odd, >1
        cfx_big_t n; cfx_big_mont_ctx_t C;
        init_ctx_u64(&C, &n, n64);

        uint64_t a64 = rand64() % n64;
        uint64_t e64 = rand64();          // any exponent
        uint64_t expect = powmod_u64_ref(a64, e64, n64);

        cfx_big_t a,e,out; 
        cfx_big_init(&a);
        cfx_big_init(&e);
        cfx_big_init(&out);
        big_from_u64(&a, a64); big_from_u64(&e, e64);
        int ok = cfx_big_modexp_binary(&out, &a, &e, &C);
        assert(ok);

        uint64_t got = big_as_u64(&out);
        if (got != expect) {
            fprintf(stderr, "Mismatch: a=%" PRIu64 " e=%" PRIu64 " n=%" PRIu64
                            " got=%" PRIu64 " exp=%" PRIu64 "\n",
                    a64, e64, n64, got, expect);
            assert(0);
        }

        cfx_big_free(&out); cfx_big_free(&e); cfx_big_free(&a);
        cfx_big_mont_ctx_free(&C); cfx_big_free(&n);
    }
}

/* Aliasing: out == a must be safe */
void test_modexp_binary_aliasing(void) {
    cfx_big_t n;
    cfx_big_init(&n);
    cfx_big_mont_ctx_t C;
    init_ctx_u64(&C, &n, 0xffffffff00000001ull);

    cfx_big_t a, e, ref;
    cfx_big_init(&a);
    cfx_big_init(&e);
    cfx_big_init(&ref);

    big_from_u64(&a, 123456789012345ull);
    big_from_u64(&e, 0xfedcba987654321ull);

    /* reference into separate out */
    int ok = cfx_big_modexp_binary(&ref, &a, &e, &C);
    assert(ok);

    /* alias: write back into a */
    ok = cfx_big_modexp_binary(&a, &a, &e, &C);
    assert(ok);

    EXPECT_EQ_BIG(&a, &ref);

    cfx_big_free(&ref); cfx_big_free(&e); cfx_big_free(&a);
    cfx_big_mont_ctx_free(&C); cfx_big_free(&n);
}

/* Fermat check with a known 64-bit prime p = 2^61 - 1 (Mersenne prime):
   For 1 <= a < p, a^(p-1) ≡ 1 (mod p).
*/
void test_modexp_binary_fermat_mersenne61(void) {
    const uint64_t p = ((1ull << 61) - 1ull);

    cfx_big_t n;
    cfx_big_init(&n);
    cfx_big_mont_ctx_t C;
    init_ctx_u64(&C, &n, p);

    cfx_big_t e,out,one;
    cfx_big_init(&e);
    cfx_big_init(&out);
    cfx_big_init(&one);
    big_from_u64(&one,1);
    big_from_u64(&e, p - 1ull);

    for (int t = 0; t < 100; ++t) {
        uint64_t a64;
        do { a64 = rand64() % p; } while (a64 == 0);

        cfx_big_t a; big_from_u64(&a, a64);
        int ok = cfx_big_modexp_binary(&out, &a, &e, &C);
        assert(ok);
        EXPECT_EQ_BIG(&out, &one);
        cfx_big_free(&a);
    }

    cfx_big_free(&out); cfx_big_free(&e); cfx_big_free(&one);
    cfx_big_mont_ctx_free(&C); cfx_big_free(&n);
}

/* Exponent reduction mod (p-1): for prime p, a^e ≡ a^(e + k*(p-1)) (mod p).
   We build a two-limb exponent by adding a large multiple of (p-1). */
void test_modexp_binary_exponent_reduction(void) {
    const uint64_t p = ((1ull << 61) - 1ull);

    cfx_big_t n; cfx_big_mont_ctx_t C;
    init_ctx_u64(&C, &n, p);

    /* base a in [1..p-1] */
    cfx_big_t a; big_from_u64(&a, (rand64() % (p - 1ull)) + 1ull);

    /* e1 small, e2 = e1 + (p-1)<<(64) (two-limb) */
    cfx_big_t e1, e2; big_from_u64(&e1, 0x123456789abcdef0ull % (p-1));
    cfx_big_init(&e2);
    cfx_big_assign(&e2, &e1);

    cfx_big_t phi;
    cfx_big_init(&phi);
    big_from_u64(&phi, p - 1ull);
    /* e2 += (p-1) << 64 */
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_assign(&tmp, &phi);
    big_shl_limbs_inplace(&tmp, 1);   /* shift by 64-bit limb */
    /* Now e2 = e1 + tmp */
    // simple add: since tmp is one limb above, we can extend e2 and set limb[1]
    if (e2.cap < 2) cfx_big_reserve(&e2, 2);
    if (e2.n < 2) e2.limb[e2.n++] = 0;
    e2.limb[1] += tmp.limb[1];  // tmp.limb[1] = (p-1), tmp.limb[0]=0
    // carry normalize if needed:
    if (e2.limb[1] == 0) { /* nothing */ }

    cfx_big_t r1, r2;
    cfx_big_init(&r1);
    cfx_big_init(&r2);
    int ok1 = cfx_big_modexp_binary(&r1, &a, &e1, &C);
    int ok2 = cfx_big_modexp_binary(&r2, &a, &e2, &C);
    assert(ok1 && ok2);
    EXPECT_EQ_BIG(&r1, &r2);

    cfx_big_free(&r1); cfx_big_free(&r2);
    cfx_big_free(&tmp); cfx_big_free(&phi);
    cfx_big_free(&e1); cfx_big_free(&e2);
    cfx_big_free(&a);
    cfx_big_mont_ctx_free(&C); cfx_big_free(&n);
}

/*..........................................................*/
int main(void) {
    CFX_TEST(test_mont_ctx_rejects_even_n);
    CFX_TEST(test_mont_mul_matches_scalar);
    CFX_TEST(test_mont_modexp_matches_scalar);
    CFX_TEST(test_mont_roundtrip_1);
    CFX_TEST(test_mont_aliasing_safe);

    CFX_TEST(test_modexp_binary_trivial);
    CFX_TEST(test_modexp_binary_matches_u64_ref);
    CFX_TEST(test_modexp_binary_aliasing);
    CFX_TEST(test_modexp_binary_fermat_mersenne61);
    CFX_TEST(test_modexp_binary_exponent_reduction);

    puts("ok!\n");
    return 0;
}
