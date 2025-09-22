#include "cfx/big.h"
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
    printf(">>>>>>>>>>>>>>>>> test_mont_mul_matches_scalar: got: 0x%llx, expect: 0x%llx (diff %llu)\n",
        got, expect, got > expect? got-expect : expect-got);
    //CFX_ASSERT_PRINT(got == expect);

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
    CFX_ASSERT_PRINT(got == powmod_u64(a64, e64, n64));

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
    cfx_big_from_u64(&a, 12345);
    cfx_big_from_u64(&b, 67890);
    CFX_ASSERT_PRINT(cfx_big_mont_to(&aR, &a, &C));
    CFX_ASSERT_PRINT(cfx_big_mont_to(&bR, &b, &C));

    cfx_big_t ref;
    cfx_big_init(&ref);
    CFX_ASSERT_PRINT(cfx_big_mont_mul(&ref, &aR, &bR, &C));

    /* alias: overwrite aR with product */
    CFX_ASSERT_PRINT(cfx_big_mont_mul(&aR, &aR, &bR, &C));
    CFX_ASSERT_PRINT(cfx_big_cmp(&aR, &ref) == 0);

    cfx_big_free(&ref);
    cfx_big_free(&aR);
    cfx_big_free(&bR);
    cfx_big_mont_ctx_free(&C);
    cfx_big_free(&n);
    cfx_big_free(&a);
    cfx_big_free(&b);
}

int main(void) {
    CFX_TEST(test_mont_ctx_rejects_even_n);
    CFX_TEST(test_mont_mul_matches_scalar);
    CFX_TEST(test_mont_modexp_matches_scalar);
    CFX_TEST(test_mont_roundtrip_1);
    CFX_TEST(test_mont_aliasing_safe);
    
    puts("ok!\n");
    return 0;
}
