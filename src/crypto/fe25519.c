/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/fe25519.h"
#include "cfx/compat.h"

/*
 * radix 2^51 for GF(2^255 - 19)
 * each limb uses 51 bits, mask = 2^51 - 1 = 0x7FFFFFFFFFFFF
 * p = 2^255 - 19, in limbs: (2^51 - 19, 2^51 - 1, 2^51 - 1, 2^51 - 1, 2^51 - 1)
 */

#define MASK51 0x7FFFFFFFFFFFFULL  /* 2^51 - 1 */

#ifdef _MSC_VER
#define uint128_t      cfx_uint128_t
#define mul64          cfx_mul64
#define add128         cfx_add128
#define shr128_51      cfx_shr128_51
#define lo128          cfx_lo128
#define u128_from_u64  cfx_u128_from_u64
#endif

const fe25519_t cfx_fe25519_zero = {{ 0, 0, 0, 0, 0 }};
const fe25519_t cfx_fe25519_one  = {{ 1, 0, 0, 0, 0 }};

void cfx_fe25519_0(fe25519_t* h) {
    h->v[0] = 0;
    h->v[1] = 0;
    h->v[2] = 0;
    h->v[3] = 0;
    h->v[4] = 0;
}

void cfx_fe25519_1(fe25519_t* h) {
    h->v[0] = 1;
    h->v[1] = 0;
    h->v[2] = 0;
    h->v[3] = 0;
    h->v[4] = 0;
}

void cfx_fe25519_copy(fe25519_t* h, const fe25519_t* f) {
    h->v[0] = f->v[0];
    h->v[1] = f->v[1];
    h->v[2] = f->v[2];
    h->v[3] = f->v[3];
    h->v[4] = f->v[4];
}

/* h = f + g, no reduction needed - headroom absorbs the growth */
void cfx_fe25519_add(fe25519_t* h, const fe25519_t* f, const fe25519_t* g) {
    h->v[0] = f->v[0] + g->v[0];
    h->v[1] = f->v[1] + g->v[1];
    h->v[2] = f->v[2] + g->v[2];
    h->v[3] = f->v[3] + g->v[3];
    h->v[4] = f->v[4] + g->v[4];
}

/* h = f - g, add 2p first to ensure non-negative result */
void cfx_fe25519_sub(fe25519_t* h, const fe25519_t* f, const fe25519_t* g) {
    const uint64_t two_p0 = 0xFFFFFFFFFFFDA;  /* 2*(2^51 - 19) */
    const uint64_t two_p1 = 0xFFFFFFFFFFFFE;  /* 2*(2^51 - 1) */

    h->v[0] = (f->v[0] + two_p0) - g->v[0];
    h->v[1] = (f->v[1] + two_p1) - g->v[1];
    h->v[2] = (f->v[2] + two_p1) - g->v[2];
    h->v[3] = (f->v[3] + two_p1) - g->v[3];
    h->v[4] = (f->v[4] + two_p1) - g->v[4];
}

/* h = -f */
void cfx_fe25519_neg(fe25519_t* h, const fe25519_t* f) {
    fe25519_t zero;
    cfx_fe25519_0(&zero);
    cfx_fe25519_sub(h, &zero, f);
}

/* carry propagation: reduce each limb to ~51 bits, final carry wraps * 19 */
void cfx_fe25519_carry(fe25519_t* h) {
    uint64_t c;

    c = h->v[0] >> 51;
    h->v[0] &= MASK51;
    h->v[1] += c;

    c = h->v[1] >> 51;
    h->v[1] &= MASK51;
    h->v[2] += c;

    c = h->v[2] >> 51;
    h->v[2] &= MASK51;
    h->v[3] += c;

    c = h->v[3] >> 51;
    h->v[3] &= MASK51;
    h->v[4] += c;

    c = h->v[4] >> 51;
    h->v[4] &= MASK51;
    h->v[0] += c * 19;  /* 2^255 ≡ 19 */

    /* one more in case c*19 overflowed */
    c = h->v[0] >> 51;
    h->v[0] &= MASK51;
    h->v[1] += c;
}

/* h = f * g, schoolbook 5x5 with inline reduction, 128-bit accumulators */
void cfx_fe25519_mul(fe25519_t* h, const fe25519_t* f, const fe25519_t* g) {
    uint64_t f0 = f->v[0], f1 = f->v[1], f2 = f->v[2], f3 = f->v[3], f4 = f->v[4];
    uint64_t g0 = g->v[0], g1 = g->v[1], g2 = g->v[2], g3 = g->v[3], g4 = g->v[4];

    /* precompute g[i] * 19 for reduction */
    uint64_t g1_19 = g1 * 19;
    uint64_t g2_19 = g2 * 19;
    uint64_t g3_19 = g3 * 19;
    uint64_t g4_19 = g4 * 19;

    /*
     * accumulate products - h[i] collects terms where exponents sum to i (mod 5)
     * high terms (i+j >= 5) reduced via *19
     *
     * h0 = f0*g0 + f1*g4*19 + f2*g3*19 + f3*g2*19 + f4*g1*19
     * h1 = f0*g1 + f1*g0    + f2*g4*19 + f3*g3*19 + f4*g2*19
     * h2 = f0*g2 + f1*g1    + f2*g0    + f3*g4*19 + f4*g3*19
     * h3 = f0*g3 + f1*g2    + f2*g1    + f3*g0    + f4*g4*19
     * h4 = f0*g4 + f1*g3    + f2*g2    + f3*g1    + f4*g0
     */

#ifdef __SIZEOF_INT128__
    /* use 128-bit arithmetic where available */
    __uint128_t h0, h1, h2, h3, h4;

    h0 = (__uint128_t)f0 * g0
       + (__uint128_t)f1 * g4_19
       + (__uint128_t)f2 * g3_19
       + (__uint128_t)f3 * g2_19
       + (__uint128_t)f4 * g1_19;

    h1 = (__uint128_t)f0 * g1
       + (__uint128_t)f1 * g0
       + (__uint128_t)f2 * g4_19
       + (__uint128_t)f3 * g3_19
       + (__uint128_t)f4 * g2_19;

    h2 = (__uint128_t)f0 * g2
       + (__uint128_t)f1 * g1
       + (__uint128_t)f2 * g0
       + (__uint128_t)f3 * g4_19
       + (__uint128_t)f4 * g3_19;

    h3 = (__uint128_t)f0 * g3
       + (__uint128_t)f1 * g2
       + (__uint128_t)f2 * g1
       + (__uint128_t)f3 * g0
       + (__uint128_t)f4 * g4_19;

    h4 = (__uint128_t)f0 * g4
       + (__uint128_t)f1 * g3
       + (__uint128_t)f2 * g2
       + (__uint128_t)f3 * g1
       + (__uint128_t)f4 * g0;

    /* carry propagation */
    uint64_t c;

    c = (uint64_t)(h0 >> 51);
    h->v[0] = (uint64_t)h0 & MASK51;
    h1 += c;

    c = (uint64_t)(h1 >> 51);
    h->v[1] = (uint64_t)h1 & MASK51;
    h2 += c;

    c = (uint64_t)(h2 >> 51);
    h->v[2] = (uint64_t)h2 & MASK51;
    h3 += c;

    c = (uint64_t)(h3 >> 51);
    h->v[3] = (uint64_t)h3 & MASK51;
    h4 += c;

    c = (uint64_t)(h4 >> 51);
    h->v[4] = (uint64_t)h4 & MASK51;
    h->v[0] += c * 19;

    c = h->v[0] >> 51;
    h->v[0] &= MASK51;
    h->v[1] += c;

#else
    /* MSVC fallback using _umul128 */
    uint128_t h0, h1, h2, h3, h4;

    h0 = mul64(f0, g0);
    h0 = add128(h0, mul64(f1, g4_19));
    h0 = add128(h0, mul64(f2, g3_19));
    h0 = add128(h0, mul64(f3, g2_19));
    h0 = add128(h0, mul64(f4, g1_19));

    h1 = mul64(f0, g1);
    h1 = add128(h1, mul64(f1, g0));
    h1 = add128(h1, mul64(f2, g4_19));
    h1 = add128(h1, mul64(f3, g3_19));
    h1 = add128(h1, mul64(f4, g2_19));

    h2 = mul64(f0, g2);
    h2 = add128(h2, mul64(f1, g1));
    h2 = add128(h2, mul64(f2, g0));
    h2 = add128(h2, mul64(f3, g4_19));
    h2 = add128(h2, mul64(f4, g3_19));

    h3 = mul64(f0, g3);
    h3 = add128(h3, mul64(f1, g2));
    h3 = add128(h3, mul64(f2, g1));
    h3 = add128(h3, mul64(f3, g0));
    h3 = add128(h3, mul64(f4, g4_19));

    h4 = mul64(f0, g4);
    h4 = add128(h4, mul64(f1, g3));
    h4 = add128(h4, mul64(f2, g2));
    h4 = add128(h4, mul64(f3, g1));
    h4 = add128(h4, mul64(f4, g0));

    /* carry propagation */
    uint64_t c;

    c = shr128_51(h0);
    h->v[0] = lo128(h0) & MASK51;
    h1 = add128(h1, u128_from_u64(c));

    c = shr128_51(h1);
    h->v[1] = lo128(h1) & MASK51;
    h2 = add128(h2, u128_from_u64(c));

    c = shr128_51(h2);
    h->v[2] = lo128(h2) & MASK51;
    h3 = add128(h3, u128_from_u64(c));

    c = shr128_51(h3);
    h->v[3] = lo128(h3) & MASK51;
    h4 = add128(h4, u128_from_u64(c));

    c = shr128_51(h4);
    h->v[4] = lo128(h4) & MASK51;
    h->v[0] += c * 19;

    c = h->v[0] >> 51;
    h->v[0] &= MASK51;
    h->v[1] += c;
#endif
}

/* h = f^2, exploits symmetry: f[i]*f[j] appears twice when i != j */
void cfx_fe25519_sqr(fe25519_t* h, const fe25519_t* f) {
    uint64_t f0 = f->v[0], f1 = f->v[1], f2 = f->v[2], f3 = f->v[3], f4 = f->v[4];

    /* double the cross terms */
    uint64_t f0_2 = f0 * 2;
    uint64_t f1_2 = f1 * 2;
    uint64_t f3_2 = f3 * 2;

    /* precompute f[i] * 19 for reduction */
    uint64_t f1_38 = f1 * 38;  /* 2 * 19 */
    uint64_t f2_19 = f2 * 19;
    uint64_t f3_38 = f3 * 38;
    uint64_t f4_19 = f4 * 19;

#ifdef __SIZEOF_INT128__
    __uint128_t h0, h1, h2, h3, h4;

    /*
     * h0 = f0^2 + 2*19*f1*f4 + 2*19*f2*f3
     * h1 = 2*f0*f1 + 2*19*f2*f4 + 19*f3^2
     * h2 = 2*f0*f2 + f1^2 + 2*19*f3*f4
     * h3 = 2*f0*f3 + 2*f1*f2 + 19*f4^2
     * h4 = 2*f0*f4 + 2*f1*f3 + f2^2
     */
    uint64_t f3_19 = f3 * 19;

    h0 = (__uint128_t)f0   * f0
       + (__uint128_t)f1_38 * f4
       + (__uint128_t)f2_19 * f3_2;

    h1 = (__uint128_t)f0_2 * f1
       + (__uint128_t)f2_19 * f4 * 2
       + (__uint128_t)f3_19 * f3;

    h2 = (__uint128_t)f0_2 * f2
       + (__uint128_t)f1   * f1
       + (__uint128_t)f3_38 * f4;

    h3 = (__uint128_t)f0_2 * f3
       + (__uint128_t)f1_2 * f2
       + (__uint128_t)f4_19 * f4;

    h4 = (__uint128_t)f0_2 * f4
       + (__uint128_t)f1_2 * f3
       + (__uint128_t)f2   * f2;

    /* carry propagation */
    uint64_t c;

    c = (uint64_t)(h0 >> 51);
    h->v[0] = (uint64_t)h0 & MASK51;
    h1 += c;

    c = (uint64_t)(h1 >> 51);
    h->v[1] = (uint64_t)h1 & MASK51;
    h2 += c;

    c = (uint64_t)(h2 >> 51);
    h->v[2] = (uint64_t)h2 & MASK51;
    h3 += c;

    c = (uint64_t)(h3 >> 51);
    h->v[3] = (uint64_t)h3 & MASK51;
    h4 += c;

    c = (uint64_t)(h4 >> 51);
    h->v[4] = (uint64_t)h4 & MASK51;
    h->v[0] += c * 19;

    c = h->v[0] >> 51;
    h->v[0] &= MASK51;
    h->v[1] += c;

#else
    /* MSVC fallback */
    uint64_t f3_19 = f3 * 19;
    uint128_t h0, h1, h2, h3, h4;

    h0 = mul64(f0, f0);
    h0 = add128(h0, mul64(f1_38, f4));
    h0 = add128(h0, mul64(f2_19, f3_2));

    h1 = mul64(f0_2, f1);
    h1 = add128(h1, mul64(f2_19, f4 * 2));
    h1 = add128(h1, mul64(f3_19, f3));

    h2 = mul64(f0_2, f2);
    h2 = add128(h2, mul64(f1, f1));
    h2 = add128(h2, mul64(f3_38, f4));

    h3 = mul64(f0_2, f3);
    h3 = add128(h3, mul64(f1_2, f2));
    h3 = add128(h3, mul64(f4_19, f4));

    h4 = mul64(f0_2, f4);
    h4 = add128(h4, mul64(f1_2, f3));
    h4 = add128(h4, mul64(f2, f2));

    uint64_t c;

    c = shr128_51(h0);
    h->v[0] = lo128(h0) & MASK51;
    h1 = add128(h1, u128_from_u64(c));

    c = shr128_51(h1);
    h->v[1] = lo128(h1) & MASK51;
    h2 = add128(h2, u128_from_u64(c));

    c = shr128_51(h2);
    h->v[2] = lo128(h2) & MASK51;
    h3 = add128(h3, u128_from_u64(c));

    c = shr128_51(h3);
    h->v[3] = lo128(h3) & MASK51;
    h4 = add128(h4, u128_from_u64(c));

    c = shr128_51(h4);
    h->v[4] = lo128(h4) & MASK51;
    h->v[0] += c * 19;

    c = h->v[0] >> 51;
    h->v[0] &= MASK51;
    h->v[1] += c;
#endif
}

/* h = f^(2^n) - repeated squaring */
void cfx_fe25519_sqr_n(fe25519_t* h, const fe25519_t* f, unsigned n) {
    cfx_fe25519_sqr(h, f);
    for (unsigned i = 1; i < n; i++) {
        cfx_fe25519_sqr(h, h);
    }
}

/* h = f * 121666, used in montgomery ladder: (A+2)/4 = 121666 for curve25519 */
void cfx_fe25519_mul121666(fe25519_t* h, const fe25519_t* f) {
#ifdef __SIZEOF_INT128__
    __uint128_t t;
    uint64_t c;

    t = (__uint128_t)f->v[0] * 121666;
    h->v[0] = (uint64_t)t & MASK51;
    c = (uint64_t)(t >> 51);

    t = (__uint128_t)f->v[1] * 121666 + c;
    h->v[1] = (uint64_t)t & MASK51;
    c = (uint64_t)(t >> 51);

    t = (__uint128_t)f->v[2] * 121666 + c;
    h->v[2] = (uint64_t)t & MASK51;
    c = (uint64_t)(t >> 51);

    t = (__uint128_t)f->v[3] * 121666 + c;
    h->v[3] = (uint64_t)t & MASK51;
    c = (uint64_t)(t >> 51);

    t = (__uint128_t)f->v[4] * 121666 + c;
    h->v[4] = (uint64_t)t & MASK51;
    c = (uint64_t)(t >> 51);

    h->v[0] += c * 19;
#else
    /* MSVC fallback */
    uint128_t t;
    uint64_t c;

    t = mul64(f->v[0], 121666);
    h->v[0] = lo128(t) & MASK51;
    c = shr128_51(t);

    t = add128(mul64(f->v[1], 121666), u128_from_u64(c));
    h->v[1] = lo128(t) & MASK51;
    c = shr128_51(t);

    t = add128(mul64(f->v[2], 121666), u128_from_u64(c));
    h->v[2] = lo128(t) & MASK51;
    c = shr128_51(t);

    t = add128(mul64(f->v[3], 121666), u128_from_u64(c));
    h->v[3] = lo128(t) & MASK51;
    c = shr128_51(t);

    t = add128(mul64(f->v[4], 121666), u128_from_u64(c));
    h->v[4] = lo128(t) & MASK51;
    c = shr128_51(t);

    h->v[0] += c * 19;
#endif
}

/* full reduction to canonical [0, p) range, required before comparison or serialization */
void cfx_fe25519_reduce(fe25519_t* h) {
    cfx_fe25519_carry(h);
    cfx_fe25519_carry(h);  /* second round ensures all limbs < 2^51 */

    /* now check if h >= p and conditionally subtract p */
    /* p = (2^51 - 19, 2^51 - 1, 2^51 - 1, 2^51 - 1, 2^51 - 1) */

    /* compute h - p, will underflow if h < p */
    int64_t t0 = (int64_t)h->v[0] - (int64_t)(MASK51 - 18);  /* -(2^51 - 19) */
    int64_t t1 = (int64_t)h->v[1] - (int64_t)MASK51;
    int64_t t2 = (int64_t)h->v[2] - (int64_t)MASK51;
    int64_t t3 = (int64_t)h->v[3] - (int64_t)MASK51;
    int64_t t4 = (int64_t)h->v[4] - (int64_t)MASK51;

    /* propagate borrows */
    t1 += t0 >> 51; t0 &= MASK51;
    t2 += t1 >> 51; t1 &= MASK51;
    t3 += t2 >> 51; t2 &= MASK51;
    t4 += t3 >> 51; t3 &= MASK51;

    /* if t4 >= 0, we had h >= p, use the reduced value */
    uint64_t mask = (uint64_t)0 - (uint64_t)(t4 >= 0);

    h->v[0] = (h->v[0] & ~mask) | ((uint64_t)t0 & mask);
    h->v[1] = (h->v[1] & ~mask) | ((uint64_t)t1 & mask);
    h->v[2] = (h->v[2] & ~mask) | ((uint64_t)t2 & mask);
    h->v[3] = (h->v[3] & ~mask) | ((uint64_t)t3 & mask);
    h->v[4] = (h->v[4] & ~mask) | ((uint64_t)t4 & mask);
}

/* constant-time conditional swap */
void cfx_fe25519_cswap(fe25519_t* f, fe25519_t* g, unsigned int b) {
    uint64_t mask = (uint64_t)0 - (uint64_t)(b != 0);
    uint64_t t;

    t = mask & (f->v[0] ^ g->v[0]); f->v[0] ^= t; g->v[0] ^= t;
    t = mask & (f->v[1] ^ g->v[1]); f->v[1] ^= t; g->v[1] ^= t;
    t = mask & (f->v[2] ^ g->v[2]); f->v[2] ^= t; g->v[2] ^= t;
    t = mask & (f->v[3] ^ g->v[3]); f->v[3] ^= t; g->v[3] ^= t;
    t = mask & (f->v[4] ^ g->v[4]); f->v[4] ^= t; g->v[4] ^= t;
}

/* constant-time conditional move: if b, h = f */
void cfx_fe25519_cmov(fe25519_t* h, const fe25519_t* f, unsigned int b) {
    uint64_t mask = (uint64_t)0 - (uint64_t)(b != 0);

    h->v[0] = (h->v[0] & ~mask) | (f->v[0] & mask);
    h->v[1] = (h->v[1] & ~mask) | (f->v[1] & mask);
    h->v[2] = (h->v[2] & ~mask) | (f->v[2] & mask);
    h->v[3] = (h->v[3] & ~mask) | (f->v[3] & mask);
    h->v[4] = (h->v[4] & ~mask) | (f->v[4] & mask);
}

/* returns 1 if f == 0, else 0 */
int cfx_fe25519_iszero(const fe25519_t* f) {
    fe25519_t t;
    cfx_fe25519_copy(&t, f);
    cfx_fe25519_reduce(&t);

    uint64_t r = t.v[0] | t.v[1] | t.v[2] | t.v[3] | t.v[4];

    /* constant-time: r == 0 iff all limbs are 0 */
    r = (r | ((uint64_t)0 - r)) >> 63;  /* 0 if r==0, 1 otherwise */
    return (int)(1 - r);
}

/* returns 1 if f is negative (LSB of reduced form) */
int cfx_fe25519_isnegative(const fe25519_t* f) {
    fe25519_t t;
    cfx_fe25519_copy(&t, f);
    cfx_fe25519_reduce(&t);
    return (int)(t.v[0] & 1);
}

/* returns 1 if f == g */
int cfx_fe25519_eq(const fe25519_t* f, const fe25519_t* g) {
    fe25519_t d;
    cfx_fe25519_sub(&d, f, g);
    return cfx_fe25519_iszero(&d);
}

/* deserialize 32 bytes (little-endian) to field element, clears bit 255 per RFC 7748 */
void cfx_fe25519_frombytes(fe25519_t* h, const uint8_t s[32]) {
    uint64_t h0 = (uint64_t)s[0]
                | ((uint64_t)s[1] << 8)
                | ((uint64_t)s[2] << 16)
                | ((uint64_t)s[3] << 24)
                | ((uint64_t)s[4] << 32)
                | ((uint64_t)s[5] << 40)
                | ((uint64_t)(s[6] & 0x07) << 48);

    uint64_t h1 = ((uint64_t)s[6] >> 3)
                | ((uint64_t)s[7] << 5)
                | ((uint64_t)s[8] << 13)
                | ((uint64_t)s[9] << 21)
                | ((uint64_t)s[10] << 29)
                | ((uint64_t)s[11] << 37)
                | ((uint64_t)(s[12] & 0x3F) << 45);

    uint64_t h2 = ((uint64_t)s[12] >> 6)
                | ((uint64_t)s[13] << 2)
                | ((uint64_t)s[14] << 10)
                | ((uint64_t)s[15] << 18)
                | ((uint64_t)s[16] << 26)
                | ((uint64_t)s[17] << 34)
                | ((uint64_t)s[18] << 42)
                | ((uint64_t)(s[19] & 0x01) << 50);

    uint64_t h3 = ((uint64_t)s[19] >> 1)
                | ((uint64_t)s[20] << 7)
                | ((uint64_t)s[21] << 15)
                | ((uint64_t)s[22] << 23)
                | ((uint64_t)s[23] << 31)
                | ((uint64_t)s[24] << 39)
                | ((uint64_t)(s[25] & 0x0F) << 47);

    uint64_t h4 = ((uint64_t)s[25] >> 4)
                | ((uint64_t)s[26] << 4)
                | ((uint64_t)s[27] << 12)
                | ((uint64_t)s[28] << 20)
                | ((uint64_t)s[29] << 28)
                | ((uint64_t)s[30] << 36)
                | ((uint64_t)(s[31] & 0x7F) << 44);  /* clear bit 255 */

    h->v[0] = h0;
    h->v[1] = h1;
    h->v[2] = h2;
    h->v[3] = h3;
    h->v[4] = h4;
}

/* serialize to 32 bytes (little-endian), fully reduces first */
void cfx_fe25519_tobytes(uint8_t s[32], const fe25519_t* h) {
    fe25519_t t;
    cfx_fe25519_copy(&t, h);
    cfx_fe25519_reduce(&t);

    uint64_t h0 = t.v[0], h1 = t.v[1], h2 = t.v[2], h3 = t.v[3], h4 = t.v[4];

    /* pack 5 x 51-bit limbs into 32 bytes */
    s[0]  = (uint8_t)(h0);
    s[1]  = (uint8_t)(h0 >> 8);
    s[2]  = (uint8_t)(h0 >> 16);
    s[3]  = (uint8_t)(h0 >> 24);
    s[4]  = (uint8_t)(h0 >> 32);
    s[5]  = (uint8_t)(h0 >> 40);
    s[6]  = (uint8_t)((h0 >> 48) | (h1 << 3));
    s[7]  = (uint8_t)(h1 >> 5);
    s[8]  = (uint8_t)(h1 >> 13);
    s[9]  = (uint8_t)(h1 >> 21);
    s[10] = (uint8_t)(h1 >> 29);
    s[11] = (uint8_t)(h1 >> 37);
    s[12] = (uint8_t)((h1 >> 45) | (h2 << 6));
    s[13] = (uint8_t)(h2 >> 2);
    s[14] = (uint8_t)(h2 >> 10);
    s[15] = (uint8_t)(h2 >> 18);
    s[16] = (uint8_t)(h2 >> 26);
    s[17] = (uint8_t)(h2 >> 34);
    s[18] = (uint8_t)(h2 >> 42);
    s[19] = (uint8_t)((h2 >> 50) | (h3 << 1));
    s[20] = (uint8_t)(h3 >> 7);
    s[21] = (uint8_t)(h3 >> 15);
    s[22] = (uint8_t)(h3 >> 23);
    s[23] = (uint8_t)(h3 >> 31);
    s[24] = (uint8_t)(h3 >> 39);
    s[25] = (uint8_t)((h3 >> 47) | (h4 << 4));
    s[26] = (uint8_t)(h4 >> 4);
    s[27] = (uint8_t)(h4 >> 12);
    s[28] = (uint8_t)(h4 >> 20);
    s[29] = (uint8_t)(h4 >> 28);
    s[30] = (uint8_t)(h4 >> 36);
    s[31] = (uint8_t)(h4 >> 44);
}

/* h = 1/f via fermat: f^(-1) = f^(p-2), p-2 = 2^255 - 21, computed via addition chain */
void cfx_fe25519_inv(fe25519_t* h, const fe25519_t* f) {
    fe25519_t t0, t1, t2, t3;

    /* addition chain for p-2 = 2^255 - 21 */

    cfx_fe25519_sqr(&t0, f);           /* t0 = f^2 */
    cfx_fe25519_sqr(&t1, &t0);         /* t1 = f^4 */
    cfx_fe25519_sqr(&t1, &t1);         /* t1 = f^8 */
    cfx_fe25519_mul(&t1, f, &t1);      /* t1 = f^9 */
    cfx_fe25519_mul(&t0, &t0, &t1);    /* t0 = f^11 */
    cfx_fe25519_sqr(&t2, &t0);         /* t2 = f^22 */
    cfx_fe25519_mul(&t1, &t1, &t2);    /* t1 = f^31 = f^(2^5 - 1) */

    cfx_fe25519_sqr_n(&t2, &t1, 5);    /* t2 = f^(2^10 - 32) */
    cfx_fe25519_mul(&t1, &t1, &t2);    /* t1 = f^(2^10 - 1) */

    cfx_fe25519_sqr_n(&t2, &t1, 10);   /* t2 = f^(2^20 - 1024) */
    cfx_fe25519_mul(&t2, &t1, &t2);    /* t2 = f^(2^20 - 1) */

    cfx_fe25519_sqr_n(&t3, &t2, 20);   /* t3 = f^(2^40 - 2^20) */
    cfx_fe25519_mul(&t2, &t2, &t3);    /* t2 = f^(2^40 - 1) */

    cfx_fe25519_sqr_n(&t2, &t2, 10);   /* t2 = f^(2^50 - 2^10) */
    cfx_fe25519_mul(&t1, &t1, &t2);    /* t1 = f^(2^50 - 1) */

    cfx_fe25519_sqr_n(&t2, &t1, 50);   /* t2 = f^(2^100 - 2^50) */
    cfx_fe25519_mul(&t2, &t1, &t2);    /* t2 = f^(2^100 - 1) */

    cfx_fe25519_sqr_n(&t3, &t2, 100);  /* t3 = f^(2^200 - 2^100) */
    cfx_fe25519_mul(&t2, &t2, &t3);    /* t2 = f^(2^200 - 1) */

    cfx_fe25519_sqr_n(&t2, &t2, 50);   /* t2 = f^(2^250 - 2^50) */
    cfx_fe25519_mul(&t1, &t1, &t2);    /* t1 = f^(2^250 - 1) */

    cfx_fe25519_sqr_n(&t1, &t1, 5);    /* t1 = f^(2^255 - 32) */
    cfx_fe25519_mul(h, &t0, &t1);      /* h = f^(2^255 - 21) = f^(p-2) */
}

/* h = f^((p-5)/8) = f^(2^252 - 3), used for square roots */
void cfx_fe25519_pow22523(fe25519_t* h, const fe25519_t* f) {
    fe25519_t t0, t1, t2;

    cfx_fe25519_sqr(&t0, f);
    cfx_fe25519_sqr(&t1, &t0);
    cfx_fe25519_sqr(&t1, &t1);
    cfx_fe25519_mul(&t1, f, &t1);
    cfx_fe25519_mul(&t0, &t0, &t1);
    cfx_fe25519_sqr(&t0, &t0);
    cfx_fe25519_mul(&t0, &t1, &t0);

    cfx_fe25519_sqr_n(&t1, &t0, 5);
    cfx_fe25519_mul(&t0, &t0, &t1);

    cfx_fe25519_sqr_n(&t1, &t0, 10);
    cfx_fe25519_mul(&t1, &t0, &t1);

    cfx_fe25519_sqr_n(&t2, &t1, 20);
    cfx_fe25519_mul(&t1, &t1, &t2);

    cfx_fe25519_sqr_n(&t1, &t1, 10);
    cfx_fe25519_mul(&t0, &t0, &t1);

    cfx_fe25519_sqr_n(&t1, &t0, 50);
    cfx_fe25519_mul(&t1, &t0, &t1);

    cfx_fe25519_sqr_n(&t2, &t1, 100);
    cfx_fe25519_mul(&t1, &t1, &t2);

    cfx_fe25519_sqr_n(&t1, &t1, 50);
    cfx_fe25519_mul(&t0, &t0, &t1);

    cfx_fe25519_sqr_n(&t0, &t0, 2);
    cfx_fe25519_mul(h, &t0, f);
}
