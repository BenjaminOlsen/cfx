/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* test_big_exp.c - Tests for big integer exponentiation */

#include "test_common.h"
#include "cfx/primes.h"

static void test_exp_edge_cases(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out);

    /* 0^0 := 1 */
    cfx_big_from_limb(&n, 0);
    cfx_big_from_limb(&p, 0);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    /* 0^k = 0 (k>0) */
    cfx_big_from_limb(&n, 0);
    cfx_big_from_limb(&p, 5);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 0));

    /* 1^p = 1 */
    cfx_big_from_limb(&n, 1);
    cfx_big_from_limb(&p, 1234567);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    /* n^0 = 1 */
    cfx_big_from_limb(&n, 42);
    cfx_big_from_limb(&p, 0);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    cfx_big_free(&out);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_small_values(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out);

    /* n^1 = n */
    cfx_big_from_limb(&n, 123456789);
    cfx_big_from_limb(&p, 1);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 1 && out.limb[0] == 123456789);

    /* 2^10 = 1024 */
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 10);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, UINT64_C(1024)));

    /* 3^5 = 243 */
    cfx_big_from_limb(&n, 3);
    cfx_big_from_limb(&p, 5);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(cfx_big_eq_u64(&out, UINT64_C(243)));

    cfx_big_free(&out);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_powers_of_two_boundaries(void) {
    cfx_big_t n, p, out;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out);

#if CFX_LIMB_BITS == 64
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 64);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    expect_limb_pattern(&out, 2, 0, 1, 0);

    cfx_big_from_limb(&p, 128);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 3);
    expect_limb_pattern(&out, 3, 1, 0, 0);

    cfx_big_from_limb(&p, 127);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    CFX_ASSERT(out.limb[0] == 0);
    CFX_ASSERT(out.limb[1] == (UINT64_C(1) << 63));
#else
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 32);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    expect_limb_pattern(&out, 2, 0, 1, 0);

    cfx_big_from_limb(&p, 64);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 3);
    expect_limb_pattern(&out, 3, 1, 0, 0);

    cfx_big_from_limb(&p, 63);
    cfx_big_exp(&out, &n, &p);
    CFX_ASSERT(out.n == 2);
    CFX_ASSERT(out.limb[0] == 0);
    CFX_ASSERT(out.limb[1] == (UINT32_C(1) << 31));
#endif

    cfx_big_free(&out);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_aliasing(void) {
    cfx_big_t n, p;
    cfx_big_init(&n);
    cfx_big_init(&p);

    cfx_big_from_limb(&n, 7);
    cfx_big_from_limb(&p, 6);
    cfx_big_exp(&n, &n, &p);
    CFX_ASSERT(n.n == 1);
    CFX_ASSERT(n.limb[0] == 117649);

#if CFX_LIMB_BITS == 64
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 80);
    cfx_big_exp(&p, &n, &p);
    CFX_ASSERT(p.n == 2);
    CFX_ASSERT(p.limb[0] == 0);
    CFX_ASSERT(p.limb[1] == (UINT64_C(1) << 16));
#else
    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 48);
    cfx_big_exp(&p, &n, &p);
    CFX_ASSERT(p.n == 2);
    CFX_ASSERT(p.limb[0] == 0);
    CFX_ASSERT(p.limb[1] == (UINT32_C(1) << 16));
#endif

    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_exp_compare_with_naive_mul(void) {
    cfx_big_t n, p, out1, out2;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&out1);
    cfx_big_init(&out2);

    cfx_limb_t bases[] = {2,3,5,10,17,1234567};
    cfx_limb_t exps[]  = {2,3,4,5,8,16};

    for (size_t i=0;i<sizeof(bases)/sizeof(bases[0]);++i) {
        for (size_t j=0;j<sizeof(exps)/sizeof(exps[0]);++j) {
            cfx_big_from_limb(&n, bases[i]);
            cfx_big_from_limb(&p, exps[j]);

            cfx_big_exp(&out1, &n, &p);

            cfx_big_from_limb(&out2, 1);
            cfx_limb_t e = exps[j];
            while (e--) {
                cfx_big_mul_auto(&out2, &n);
            }

            CFX_ASSERT(cfx_big_cmp(&out1, &out2) == 0);
        }
    }

    cfx_big_free(&out2);
    cfx_big_free(&out1);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_big_prime(void) {
    cfx_big_t b1, b2;
    cfx_big_init(&b1);
    cfx_big_init(&b2);

    unsigned cnt = 0;
    for (size_t i = 0; i < cfx_primes_len; i += 100) {
        ++cnt;
        cfx_limb_t val = (cfx_limb_t)cfx_primes[i];
        cfx_big_assign_sm(&b1, val);
        cfx_big_copy(&b2, &b1);
        cfx_big_add_sm(&b2, 1);
        CFX_ASSERT(cfx_big_is_prime(&b1) == 1);
        if(b2.limb[0] != 3) CFX_ASSERT(cfx_big_is_prime(&b2) == 0);
    }

    cfx_big_from_dec(&b1, "1361129467683753853853498429727072845819");
    CFX_ASSERT(cfx_big_is_prime(&b1) == 1);

    cfx_big_from_dec(&b1, "18446744073709551557");
    CFX_ASSERT(cfx_big_is_prime(&b1) == 1);

    cfx_big_free(&b2);
    cfx_big_free(&b1);
}

static void test_primepow(void) {
    cfx_big_t out;
    cfx_big_init(&out);

    cfx_big_pow_sm(&out, 2, 10);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1024));

    cfx_big_pow_sm(&out, 3, 5);
    CFX_ASSERT(cfx_big_eq_u64(&out, 243));

    cfx_big_pow_sm(&out, 7, 0);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    cfx_big_pow_sm(&out, 5, 1);
    CFX_ASSERT(cfx_big_eq_u64(&out, 5));

    cfx_big_pow_sm(&out, 13, 7);
    CFX_ASSERT(cfx_big_eq_u64(&out, 62748517));

    /* compare with cfx_big_exp */
    cfx_big_t base, exp, ref;
    cfx_big_init(&base);
    cfx_big_init(&exp);
    cfx_big_init(&ref);

    cfx_big_from_limb(&base, 17);
    cfx_big_from_limb(&exp, 11);
    cfx_big_exp(&ref, &base, &exp);
    cfx_big_pow_sm(&out, 17, 11);
    CFX_ASSERT(cfx_big_cmp(&out, &ref) == 0);

    cfx_big_free(&ref);
    cfx_big_free(&exp);
    cfx_big_free(&base);
    cfx_big_free(&out);
}

static void test_exp_u64(void) {
    cfx_big_t n, out;
    cfx_big_init(&n);
    cfx_big_init(&out);

    cfx_big_from_limb(&n, 2);
    cfx_big_exp_u64(&out, &n, 10);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1024));

    cfx_big_from_limb(&n, 7);
    cfx_big_exp_u64(&out, &n, 6);
    CFX_ASSERT(cfx_big_eq_u64(&out, 117649));

    cfx_big_from_limb(&n, 42);
    cfx_big_exp_u64(&out, &n, 0);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    cfx_big_from_limb(&n, 12345);
    cfx_big_exp_u64(&out, &n, 1);
    CFX_ASSERT(cfx_big_eq_u64(&out, 12345));

    cfx_big_t p, ref;
    cfx_big_init(&p);
    cfx_big_init(&ref);

    cfx_big_from_limb(&n, 13);
    cfx_big_from_limb(&p, 8);
    cfx_big_exp(&ref, &n, &p);
    cfx_big_exp_u64(&out, &n, 8);
    CFX_ASSERT(cfx_big_cmp(&out, &ref) == 0);

    cfx_big_free(&ref);
    cfx_big_free(&p);
    cfx_big_free(&out);
    cfx_big_free(&n);
}

static void test_mod_exp(void) {
    cfx_big_t n, p, m, out;
    cfx_big_init(&n);
    cfx_big_init(&p);
    cfx_big_init(&m);
    cfx_big_init(&out);

    cfx_big_from_limb(&n, 2);
    cfx_big_from_limb(&p, 10);
    cfx_big_from_limb(&m, 1000);
    cfx_big_mod_exp(&out, &n, &p, &m);
    CFX_ASSERT(cfx_big_eq_u64(&out, 24));  /* 1024 mod 1000 */

    cfx_big_from_limb(&n, 3);
    cfx_big_from_limb(&p, 7);
    cfx_big_from_limb(&m, 100);
    cfx_big_mod_exp(&out, &n, &p, &m);
    CFX_ASSERT(cfx_big_eq_u64(&out, 87));  /* 2187 mod 100 */

    /* fermat's little theorem: 7^12 = 1 mod 13 */
    cfx_big_from_limb(&n, 7);
    cfx_big_from_limb(&p, 12);
    cfx_big_from_limb(&m, 13);
    cfx_big_mod_exp(&out, &n, &p, &m);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    cfx_big_from_limb(&n, 123);
    cfx_big_from_limb(&p, 0);
    cfx_big_from_limb(&m, 7);
    cfx_big_mod_exp(&out, &n, &p, &m);
    CFX_ASSERT(cfx_big_eq_u64(&out, 1));

    cfx_big_from_limb(&n, 0);
    cfx_big_from_limb(&p, 5);
    cfx_big_from_limb(&m, 17);
    cfx_big_mod_exp(&out, &n, &p, &m);
    CFX_ASSERT(cfx_big_eq_u64(&out, 0));

    cfx_big_from_hex(&n, "deadbeef");
    cfx_big_from_limb(&p, 100);
    cfx_big_from_hex(&m, "1000000007");
    cfx_big_mod_exp(&out, &n, &p, &m);
    CFX_ASSERT(cfx_big_cmp(&out, &m) < 0);

    cfx_big_free(&out);
    cfx_big_free(&m);
    cfx_big_free(&p);
    cfx_big_free(&n);
}

static void test_mod_exp_rsa_style(void) {
    /* rsa round-trip: (m^e)^d mod n = m, with p=61 q=53 */
    cfx_big_t msg, e, d, n, cipher, decrypted;
    cfx_big_init(&msg);
    cfx_big_init(&e);
    cfx_big_init(&d);
    cfx_big_init(&n);
    cfx_big_init(&cipher);
    cfx_big_init(&decrypted);

    cfx_big_from_limb(&n, 3233);
    cfx_big_from_limb(&e, 17);
    cfx_big_from_limb(&d, 2753);
    cfx_big_from_limb(&msg, 42);

    cfx_big_mod_exp(&cipher, &msg, &e, &n);
    cfx_big_mod_exp(&decrypted, &cipher, &d, &n);
    CFX_ASSERT(cfx_big_cmp(&msg, &decrypted) == 0);

    cfx_big_free(&decrypted);
    cfx_big_free(&cipher);
    cfx_big_free(&n);
    cfx_big_free(&d);
    cfx_big_free(&e);
    cfx_big_free(&msg);
}

int main(void) {
    CFX_TEST(test_exp_edge_cases);
    CFX_TEST(test_exp_small_values);
    CFX_TEST(test_exp_powers_of_two_boundaries);
    CFX_TEST(test_exp_aliasing);
    CFX_TEST(test_exp_compare_with_naive_mul);
    CFX_TEST(test_big_prime);
    CFX_TEST(test_primepow);
    CFX_TEST(test_exp_u64);
    CFX_TEST(test_mod_exp);
    CFX_TEST(test_mod_exp_rsa_style);
    puts("OK");
    return 0;
}
