#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <assert.h>

#include "cfx/arith.h"
#include "cfx/macros.h"

#ifndef CFX_ASSERT
#define CFX_ASSERT(expr) assert(expr)
#endif

/* Wide type for reference math */
#if (CFX_LIMB_BITS == 32)
    typedef uint64_t cfx_wide_t;
    #define CFX_HAVE_WIDE_TEST 1
#elif (CFX_LIMB_BITS == 64) && CFX_HAS_UINT128
    typedef __uint128_t cfx_wide_t;
    #define CFX_HAVE_WIDE_TEST 1
#endif

static cfx_acc_t make_acc(cfx_limb_t hi, cfx_limb_t lo) {
#if defined(CFX_ACC_NATIVE)
    cfx_acc_t a = ((cfx_acc_t)hi << CFX_LIMB_BITS) | (cfx_acc_t)lo;
    return a;
#else
    cfx_acc_t a;
    a.hi = hi;
    a.lo = lo;
    return a;
#endif
}

static void dump_acc(const char *label, const cfx_acc_t *a) {
    (void)label;
    (void)a;
#if 0
    printf("%s: hi=" CFX_PRI0xLIMB " lo=" CFX_PRI0xLIMB "\n",
           label, cfx_acc_hi(*a), cfx_acc_lo(*a));
#endif
}



static void test_acc_zero_and_from_lo(void) {
    cfx_acc_t a = cfx_acc_zero();
    CFX_ASSERT(cfx_acc_lo(a) == 0);
    CFX_ASSERT(cfx_acc_hi(a) == 0);

    cfx_limb_t vals[] = {0, 1, 7, 12345u, (cfx_limb_t)CFX_LIMB_MAX};
    for (size_t i = 0; i < sizeof(vals)/sizeof(vals[0]); ++i) {
        cfx_acc_t b = cfx_acc_from_lo(vals[i]);
        CFX_ASSERT(cfx_acc_lo(b) == vals[i]);
        CFX_ASSERT(cfx_acc_hi(b) == 0);
    }
}

static void test_acc_lo_hi_roundtrip(void) {
    /* A few arbitrary patterns */
    cfx_limb_t his[] = {
        0,
        1,
        (cfx_limb_t)0x12345678u,
        (cfx_limb_t)(CFX_LIMB_MAX >> 1),
        (cfx_limb_t)CFX_LIMB_MAX
    };
    cfx_limb_t los[] = {
        0,
        1,
        (cfx_limb_t)0x87654321u,
        (cfx_limb_t)(CFX_LIMB_MAX >> 1),
        (cfx_limb_t)CFX_LIMB_MAX
    };

    for (size_t i = 0; i < sizeof(his)/sizeof(his[0]); ++i) {
        cfx_acc_t a = make_acc(his[i], los[i]);
        CFX_ASSERT(cfx_acc_hi(a) == his[i]);
        CFX_ASSERT(cfx_acc_lo(a) == los[i]);
    }
}

static void test_mul_wide_basic(void) {
    struct {
        cfx_limb_t x, y;
    } cases[] = {
        {0, 0},
        {0, 1},
        {1, 0},
        {1, 1},
        {2, 3},
        {7, 9},
        {12345u, 6789u},
        { (cfx_limb_t)0xffffffffu, (cfx_limb_t)0xffffffffu },
        { (cfx_limb_t)CFX_LIMB_MAX, (cfx_limb_t)1 },
        { (cfx_limb_t)CFX_LIMB_MAX, (cfx_limb_t)CFX_LIMB_MAX }
    };

    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
        cfx_limb_t x = cases[i].x;
        cfx_limb_t y = cases[i].y;
        cfx_limb_t hi, lo;
        cfx_mul_wide(x, y, &hi, &lo);

#if CFX_LIMB_BITS == 32
        cfx_wide_t exp = (cfx_wide_t)x * (cfx_wide_t)y;
        cfx_limb_t exp_lo = (cfx_limb_t)exp;
        cfx_limb_t exp_hi = (cfx_limb_t)(exp >> 32);
        CFX_ASSERT(lo == exp_lo);
        CFX_ASSERT(hi == exp_hi);
#elif (CFX_LIMB_BITS == 64) && CFX_HAS_UINT128
        cfx_wide_t exp = (cfx_wide_t)x * (cfx_wide_t)y;
        cfx_limb_t exp_lo = (cfx_limb_t)exp;
        cfx_limb_t exp_hi = (cfx_limb_t)(exp >> 64);
        CFX_ASSERT(lo == exp_lo);
        CFX_ASSERT(hi == exp_hi);
#else
        /* No 128-bit; at least check low limb matches x*y mod 2^64 */
        cfx_limb_t exp_lo = (cfx_limb_t)((uint64_t)x * (uint64_t)y);
        CFX_ASSERT(lo == exp_lo);
        (void)hi; /* high limb still tested indirectly via other functions */
#endif
    }
}

static void test_acc_mul_vs_mul_wide(void) {
    cfx_limb_t xs[] = {0, 1, 2, 3, 17, 12345u, (cfx_limb_t)CFX_LIMB_MAX};
    cfx_limb_t ys[] = {0, 1, 5, 9, 23, 99991u, (cfx_limb_t)CFX_LIMB_MAX};

    for (size_t i = 0; i < sizeof(xs)/sizeof(xs[0]); ++i) {
        for (size_t j = 0; j < sizeof(ys)/sizeof(ys[0]); ++j) {
            cfx_limb_t x = xs[i];
            cfx_limb_t y = ys[j];

            cfx_limb_t hi, lo;
            cfx_mul_wide(x, y, &hi, &lo);

            cfx_acc_t a;
            cfx_acc_mul(&a, x, y);

            CFX_ASSERT(cfx_acc_lo(a) == lo);
            CFX_ASSERT(cfx_acc_hi(a) == hi);
        }
    }
}

static void test_acc_add_lo_hi(void) {
    cfx_acc_t a = make_acc((cfx_limb_t)0x12u, (cfx_limb_t)0x34u);

    /* add_lo no carry */
    cfx_acc_t b = a;
    cfx_limb_t add = 1;
    cfx_limb_t old_lo = cfx_acc_lo(b);
    cfx_limb_t old_hi = cfx_acc_hi(b);
    cfx_acc_add_lo(&b, add);
    cfx_limb_t exp_lo = old_lo + add;
    cfx_limb_t exp_hi = old_hi + (exp_lo < old_lo ? 1 : 0);
    CFX_ASSERT(cfx_acc_lo(b) == exp_lo);
    CFX_ASSERT(cfx_acc_hi(b) == exp_hi);

    /* add_lo with carry */
    b = make_acc(5, (cfx_limb_t)CFX_LIMB_MAX);
    add = 1;
    old_lo = cfx_acc_lo(b);
    old_hi = cfx_acc_hi(b);
    cfx_acc_add_lo(&b, add);
    exp_lo = old_lo + add;
    exp_hi = old_hi + (exp_lo < old_lo ? 1 : 0);
    CFX_ASSERT(cfx_acc_lo(b) == exp_lo);
    CFX_ASSERT(cfx_acc_hi(b) == exp_hi);

    /* add_hi is just hi += add_hi */
    b = make_acc(7, 42);
    cfx_acc_add_hi(&b, 3);
    CFX_ASSERT(cfx_acc_hi(b) == 10);
    CFX_ASSERT(cfx_acc_lo(b) == 42);
}

static void test_acc_mac(void) {
    cfx_acc_t acc = cfx_acc_zero();

    /* acc += 2*3 */
    cfx_acc_mac(&acc, 2, 3);
    CFX_ASSERT(cfx_acc_lo(acc) == 6);
    CFX_ASSERT(cfx_acc_hi(acc) == 0);

    /* acc += large product and compare via wide type when available */
    cfx_limb_t x = (cfx_limb_t)123456789u;
    cfx_limb_t y = (cfx_limb_t)987654321u;
    cfx_acc_mac(&acc, x, y);

#ifdef CFX_HAVE_WIDE_TEST
    cfx_wide_t ref = (cfx_wide_t)2 * 3 + (cfx_wide_t)x * (cfx_wide_t)y;
    cfx_limb_t exp_lo = (cfx_limb_t)ref;
    cfx_limb_t exp_hi = (cfx_limb_t)(ref >> CFX_LIMB_BITS);
    CFX_ASSERT(cfx_acc_lo(acc) == exp_lo);
    CFX_ASSERT(cfx_acc_hi(acc) == exp_hi);
#endif
}

static void test_acc_cmp(void) {
    cfx_acc_t a = make_acc(1, 2);
    cfx_acc_t b = make_acc(1, 3);
    cfx_acc_t c = make_acc(2, 0);
    cfx_acc_t d = make_acc(1, 2);

    CFX_ASSERT(cfx_acc_lt(&a, &b));
    CFX_ASSERT(cfx_acc_le(&a, &b));
    CFX_ASSERT(!cfx_acc_eq(&a, &b));

    CFX_ASSERT(cfx_acc_lt(&b, &c));
    CFX_ASSERT(cfx_acc_le(&b, &c));

    CFX_ASSERT(cfx_acc_eq(&a, &d));
    CFX_ASSERT(!cfx_acc_lt(&a, &d));
    CFX_ASSERT(cfx_acc_le(&a, &d));
}

static void test_acc_mul_eq(void) {
    /* Simple cases, including aliasing (a *= a). */

    /* case 1: a = 1, b = 1 */
    cfx_acc_t a = make_acc(0, 1);
    cfx_acc_t b = make_acc(0, 1);
    cfx_acc_mul_eq(&a, &b);
    CFX_ASSERT(cfx_acc_hi(a) == 0);
    CFX_ASSERT(cfx_acc_lo(a) == 1);

#ifdef CFX_HAVE_WIDE_TEST
    /* case 2: random-ish values, compare against reference wide product */
    a = make_acc(0x1234u, 0x5678u);
    b = make_acc(0x9abcu, 0xdef0u);

    cfx_wide_t A = ((cfx_wide_t)cfx_acc_hi(a) << CFX_LIMB_BITS) |
                   (cfx_wide_t)cfx_acc_lo(a);
    cfx_wide_t B = ((cfx_wide_t)cfx_acc_hi(b) << CFX_LIMB_BITS) |
                   (cfx_wide_t)cfx_acc_lo(b);
    cfx_wide_t P = A * B;

    cfx_acc_mul_eq(&a, &b);

    cfx_limb_t exp_lo = (cfx_limb_t)P;
    cfx_limb_t exp_hi = (cfx_limb_t)(P >> CFX_LIMB_BITS);

    CFX_ASSERT(cfx_acc_lo(a) == exp_lo);
    CFX_ASSERT(cfx_acc_hi(a) == exp_hi);

    /* case 3: squaring (a *= a) */
    a = make_acc(0x1111u, 0x2222u);
    A = ((cfx_wide_t)cfx_acc_hi(a) << CFX_LIMB_BITS) |
        (cfx_wide_t)cfx_acc_lo(a);
    P = A * A;
    cfx_acc_mul_eq(&a, &a);  /* alias */
    exp_lo = (cfx_limb_t)P;
    exp_hi = (cfx_limb_t)(P >> CFX_LIMB_BITS);
    CFX_ASSERT(cfx_acc_lo(a) == exp_lo);
    CFX_ASSERT(cfx_acc_hi(a) == exp_hi);
#endif
}

static void test_acc_divrem(void) {
    /* Basic cases including u < v, u == v, general case. */

    /* case 1: u < v */
    cfx_acc_t u = make_acc(0, 10);
    cfx_acc_t v = make_acc(0, 20);
    cfx_acc_t q, r;
    cfx_acc_divrem(&q, &r, &u, &v);
    CFX_ASSERT(cfx_acc_hi(q) == 0 && cfx_acc_lo(q) == 0);
    CFX_ASSERT(cfx_acc_hi(r) == cfx_acc_hi(u));
    CFX_ASSERT(cfx_acc_lo(r) == cfx_acc_lo(u));

    /* case 2: u == v */
    u = make_acc(0, 1234);
    v = make_acc(0, 1234);
    cfx_acc_divrem(&q, &r, &u, &v);
    CFX_ASSERT(cfx_acc_hi(q) == 0 && cfx_acc_lo(q) == 1);
    CFX_ASSERT(cfx_acc_hi(r) == 0 && cfx_acc_lo(r) == 0);

    /* case 3: general values with wide reference if available */
#ifdef CFX_HAVE_WIDE_TEST
    struct {
        cfx_limb_t u_hi, u_lo;
        cfx_limb_t v_hi, v_lo;
    } cases[] = {
        {0, 100, 0, 7},
        {0, 1000, 0, 10},
        {1, 0, 0, 3},
        {0x10u, 0, 0, 0x123u},
        {0x1234u, 0x5678u, 0, 0x9u},
        {0x12345678u, 0x9abcdef0u, 0x1u, 0x2345u}
    };

    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
        u = make_acc(cases[i].u_hi, cases[i].u_lo);
        v = make_acc(cases[i].v_hi, cases[i].v_lo);

        if (cases[i].v_hi == 0 && cases[i].v_lo == 0)
            continue; /* skip div-by-zero here */

        cfx_acc_divrem(&q, &r, &u, &v);

        cfx_wide_t U = ((cfx_wide_t)cases[i].u_hi << CFX_LIMB_BITS) |
                       (cfx_wide_t)cases[i].u_lo;
        cfx_wide_t V = ((cfx_wide_t)cases[i].v_hi << CFX_LIMB_BITS) |
                       (cfx_wide_t)cases[i].v_lo;

        cfx_wide_t Q = U / V;
        cfx_wide_t R = U % V;

        cfx_limb_t exp_q_lo = (cfx_limb_t)Q;
        cfx_limb_t exp_q_hi = (cfx_limb_t)(Q >> CFX_LIMB_BITS);
        cfx_limb_t exp_r_lo = (cfx_limb_t)R;
        cfx_limb_t exp_r_hi = (cfx_limb_t)(R >> CFX_LIMB_BITS);

        CFX_ASSERT(cfx_acc_lo(q) == exp_q_lo);
        CFX_ASSERT(cfx_acc_hi(q) == exp_q_hi);
        CFX_ASSERT(cfx_acc_lo(r) == exp_r_lo);
        CFX_ASSERT(cfx_acc_hi(r) == exp_r_hi);
    }
#else
    /* No wide reference type available:
       check basic invariants: r < v and u ~ q*v + r (mod B^2). */

    u = make_acc(0, 1000);
    v = make_acc(0, 37);
    cfx_acc_divrem(&q, &r, &u, &v);

    CFX_ASSERT(cfx_acc_lt(&r, &v));

    /* Check u = q*v + r modulo base^2 via mul/add */
    cfx_acc_t prod = v;
    cfx_acc_t qq = q;

    /* naive multiply qq * v via double-and-add on acc (small sanity check) */
    cfx_acc_t acc_zero = cfx_acc_zero();
    prod = acc_zero;
    int total_bits = 2 * CFX_LIMB_BITS;
    for (int i = 0; i < total_bits; ++i) {
        /* if bit i of qq is set, prod += v << i (only for small tests) */
        cfx_limb_t bit;
        if (i >= CFX_LIMB_BITS)
            bit = (cfx_acc_hi(qq) >> (i - CFX_LIMB_BITS)) & (cfx_limb_t)1;
        else
            bit = (cfx_acc_lo(qq) >> i) & (cfx_limb_t)1;

        if (bit) {
            /* v shifted by i bits within 2 limbs – only safe for small i.
               We only rely on this in this limited test case. */
            cfx_acc_t tmp = v;
            for (int s = 0; s < i; ++s) {
                cfx_limb_t carry = cfx_acc_lo(tmp) >> (CFX_LIMB_BITS - 1);
                cfx_limb_t lo2   = cfx_acc_lo(tmp) << 1;
                cfx_limb_t hi2   = (cfx_acc_hi(tmp) << 1) | carry;
                tmp = make_acc(hi2, lo2);
            }
            /* prod += tmp */
            cfx_acc_add_lo(&prod, cfx_acc_lo(tmp));
            cfx_acc_add_hi(&prod, cfx_acc_hi(tmp));
        }
    }

    /* prod += r */
    cfx_acc_add_lo(&prod, cfx_acc_lo(r));
    cfx_acc_add_hi(&prod, cfx_acc_hi(r));

    /* prod should match u modulo base^2 (here all small enough) */
    CFX_ASSERT(cfx_acc_lo(prod) == cfx_acc_lo(u));
    CFX_ASSERT(cfx_acc_hi(prod) == cfx_acc_hi(u));
#endif
}

int main(void) {
    CFX_TEST(test_acc_zero_and_from_lo);
    CFX_TEST(test_acc_lo_hi_roundtrip);
    CFX_TEST(test_mul_wide_basic);
    CFX_TEST(test_acc_mul_vs_mul_wide);
    CFX_TEST(test_acc_add_lo_hi);
    CFX_TEST(test_acc_mac);
    CFX_TEST(test_acc_cmp);
    CFX_TEST(test_acc_mul_eq);
    CFX_TEST(test_acc_divrem);
    return 0;
}
