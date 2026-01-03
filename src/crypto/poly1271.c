/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/poly1271.h"
#include "cfx/memory.h"

#include <string.h>

/*
 * Poly1271: polynomial MAC over GF(2^127 - 1)
 *
 * The Mersenne prime M127 = 2^127 - 1 allows very fast reduction:
 * if x = high * 2^127 + low, then x ≡ high + low (mod M127)
 *
 * We represent 127-bit values as two 64-bit limbs:
 *   value = lo + hi * 2^64, where hi has at most 63 bits
 */

#define M127_HI 0x7FFFFFFFFFFFFFFFULL
#define M127_LO 0xFFFFFFFFFFFFFFFFULL

static inline uint64_t load64_le(const uint8_t* p) {
    return (uint64_t)p[0]       | ((uint64_t)p[1] << 8)  |
           ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}

static inline void store64_le(uint8_t* p, uint64_t v) {
    p[0] = (uint8_t)(v);       p[1] = (uint8_t)(v >> 8);
    p[2] = (uint8_t)(v >> 16); p[3] = (uint8_t)(v >> 24);
    p[4] = (uint8_t)(v >> 32); p[5] = (uint8_t)(v >> 40);
    p[6] = (uint8_t)(v >> 48); p[7] = (uint8_t)(v >> 56);
}

#if defined(__SIZEOF_INT128__)
typedef unsigned __int128 uint128_t;

static inline void mul64(uint64_t a, uint64_t b, uint64_t* lo, uint64_t* hi) {
    uint128_t r = (uint128_t)a * b;
    *lo = (uint64_t)r;
    *hi = (uint64_t)(r >> 64);
}

#elif defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>

static inline void mul64(uint64_t a, uint64_t b, uint64_t* lo, uint64_t* hi) {
    *lo = _umul128(a, b, hi);
}

#else
static inline void mul64(uint64_t a, uint64_t b, uint64_t* lo, uint64_t* hi) {
    uint32_t a0 = (uint32_t)a, a1 = (uint32_t)(a >> 32);
    uint32_t b0 = (uint32_t)b, b1 = (uint32_t)(b >> 32);

    uint64_t p00 = (uint64_t)a0 * b0;
    uint64_t p01 = (uint64_t)a0 * b1;
    uint64_t p10 = (uint64_t)a1 * b0;
    uint64_t p11 = (uint64_t)a1 * b1;

    uint64_t mid = p01 + (p00 >> 32) + (uint32_t)p10;
    *lo = (p00 & 0xFFFFFFFF) | (mid << 32);
    *hi = p11 + (mid >> 32) + (p10 >> 32);
}
#endif

/*
 * All carry detection uses direct comparisons like (result < operand)
 * which compile to constant-time code on modern compilers. (fingers crossed)
 */

/* full reduction mod 2^127 - 1 (constant-time) */
static void reduce_m127_full(uint64_t acc[3]) {
    uint64_t lo, hi, carry, overflow, old_lo;

    uint64_t high_bits = (acc[1] >> 63) | (acc[2] << 1);
    uint64_t high_hi = acc[2] >> 63;

    lo = acc[0];
    hi = acc[1] & 0x7FFFFFFFFFFFFFFFULL;

    /* fold high bits: 2^127 ≡ 1 (mod M127) */
    old_lo = lo;
    lo += high_bits;
    carry = (lo < old_lo);
    hi += high_hi + carry;

    overflow = hi >> 63;
    hi &= 0x7FFFFFFFFFFFFFFFULL;
    old_lo = lo;
    lo += overflow;
    carry = (lo < old_lo);
    hi += carry;

    /* second overflow possible but rare */
    overflow = hi >> 63;
    hi &= 0x7FFFFFFFFFFFFFFFULL;
    lo += overflow;

    acc[0] = lo;
    acc[1] = hi;
    acc[2] = 0;
}

/* partial reduction: fold acc[2] into lower limbs (constant-time) */
static inline void reduce_m127_partial(uint64_t acc[3]) {
    /* 2^128 ≡ 2 (mod M127), so acc[2] * 2^128 ≡ acc[2] * 2 */
    uint64_t fold = acc[2] << 1;
    uint64_t old_lo = acc[0];
    acc[0] += fold;
    /* carry detection: acc[0] < old_lo after add means overflow */
    uint64_t carry = (acc[0] < old_lo);  /* 0 or 1, constant-time comparison */
    acc[1] += carry;
    acc[2] = 0;
}

/*
 * Multiply 3-limb accumulator by 2-limb r.
 * Input: acc can have up to ~136 bits (acc[2] small from lazy reduction)
 * Output: partially reduced, acc[2] can be non-zero
 */
static void mul_r_lazy(uint64_t acc[3], const uint64_t r[2]) {
    uint64_t r0 = r[0], r1 = r[1];
    uint64_t a0 = acc[0], a1 = acc[1], a2 = acc[2];

    uint64_t p0_lo, p0_hi, p1_lo, p1_hi, p2_lo, p2_hi, p3_lo, p3_hi;
    uint64_t p4_lo, p4_hi, p5_lo, p5_hi;

    /* 2-limb × 2-limb = 4 products */
    mul64(a0, r0, &p0_lo, &p0_hi);
    mul64(a0, r1, &p1_lo, &p1_hi);
    mul64(a1, r0, &p2_lo, &p2_hi);
    mul64(a1, r1, &p3_lo, &p3_hi);

    /* a2 × r (a2 is small, typically 0-8) */
    mul64(a2, r0, &p4_lo, &p4_hi);
    mul64(a2, r1, &p5_lo, &p5_hi);

    uint64_t res0, res1, res2, res3, res4;
    uint64_t carry;

    res0 = p0_lo;

    res1 = p0_hi;
    carry = 0;
    res1 += p1_lo; carry += (res1 < p1_lo);
    res1 += p2_lo; carry += (res1 < p2_lo);

    res2 = p1_hi + carry;
    carry = (res2 < p1_hi);
    res2 += p2_hi; carry += (res2 < p2_hi);
    res2 += p3_lo; carry += (res2 < p3_lo);
    res2 += p4_lo; carry += (res2 < p4_lo);  /* a2*r0 low */

    res3 = p3_hi + carry;
    carry = (res3 < p3_hi);
    res3 += p4_hi; carry += (res3 < p4_hi);  /* a2*r0 high */
    res3 += p5_lo; carry += (res3 < p5_lo);  /* a2*r1 low */

    res4 = p5_hi + carry;  /* a2*r1 high */

    /* fold bits [127..] back using 2^127 ≡ 1 (mod M127) */
    uint64_t high0 = (res1 >> 63) | (res2 << 1);
    uint64_t high1 = (res2 >> 63) | (res3 << 1);
    uint64_t high2 = (res3 >> 63) | (res4 << 1);
    uint64_t high3 = res4 >> 63;

    acc[0] = res0;
    acc[1] = res1 & 0x7FFFFFFFFFFFFFFFULL;

    /* add folded high bits */
    carry = 0;
    acc[0] += high0; carry = (acc[0] < high0);
    acc[1] += high1 + carry; carry = (acc[1] < high1 + carry);
    acc[2] = high2 + high3 + carry;

    /* always do partial reduction (constant-time) */
    reduce_m127_partial(acc);
}

/* multiply with full reduction (used periodically and at end) */
static void mul_r_full(uint64_t acc[3], const uint64_t r[2]) {
    mul_r_lazy(acc, r);
    reduce_m127_full(acc);
}

/*
 * Multiply 2-limb value by 2-limb r, return 3-limb result (partially reduced).
 * Used for computing m * r in 2-block processing.
 */
static void mul_2limb(uint64_t out[3], const uint64_t a[2], const uint64_t r[2]) {
    uint64_t a0 = a[0], a1 = a[1];
    uint64_t r0 = r[0], r1 = r[1];

    uint64_t p0_lo, p0_hi, p1_lo, p1_hi, p2_lo, p2_hi, p3_lo, p3_hi;

    mul64(a0, r0, &p0_lo, &p0_hi);
    mul64(a0, r1, &p1_lo, &p1_hi);
    mul64(a1, r0, &p2_lo, &p2_hi);
    mul64(a1, r1, &p3_lo, &p3_hi);

    uint64_t res0, res1, res2, res3;
    uint64_t carry;

    res0 = p0_lo;

    res1 = p0_hi;
    carry = 0;
    res1 += p1_lo; carry += (res1 < p1_lo);
    res1 += p2_lo; carry += (res1 < p2_lo);

    res2 = p1_hi + carry;
    carry = (res2 < p1_hi);
    res2 += p2_hi; carry += (res2 < p2_hi);
    res2 += p3_lo; carry += (res2 < p3_lo);

    res3 = p3_hi + carry;

    /* fold bits [127..] back */
    uint64_t high0 = (res1 >> 63) | (res2 << 1);
    uint64_t high1 = (res2 >> 63) | (res3 << 1);
    uint64_t high2 = res3 >> 63;

    out[0] = res0;
    out[1] = res1 & 0x7FFFFFFFFFFFFFFFULL;

    carry = 0;
    out[0] += high0; carry = (out[0] < high0);
    out[1] += high1 + carry; carry = (out[1] < high1 + carry);
    out[2] = high2 + carry;
}

/*
 * Process 2 blocks at once: acc = (acc + m1) * r² + m2 * r
 * This allows better ILP as the two multiplications are independent.
 */
static void process_2blocks(uint64_t acc[3], const uint8_t* m1, const uint8_t* m2,
                            const uint64_t r[2], const uint64_t r2[2]) {
    /* load m1 */
    uint64_t b1_0 = load64_le(m1);
    uint64_t b1_1 = (load64_le(m1 + 7) >> 8) | (0x01ULL << 56);

    /* load m2 */
    uint64_t b2_0 = load64_le(m2);
    uint64_t b2_1 = (load64_le(m2 + 7) >> 8) | (0x01ULL << 56);

    /* t1 = acc + m1 */
    uint64_t t1[3];
    uint64_t carry;
    t1[0] = acc[0] + b1_0;
    carry = (t1[0] < acc[0]);
    t1[1] = acc[1] + b1_1 + carry;
    carry = (t1[1] < acc[1] + carry);
    t1[2] = acc[2] + carry;

    /* prod1 = t1 * r² (need 3-limb × 2-limb multiply) */
    uint64_t prod1[3];
    {
        uint64_t a0 = t1[0], a1 = t1[1], a2 = t1[2];
        uint64_t rr0 = r2[0], rr1 = r2[1];

        uint64_t p0_lo, p0_hi, p1_lo, p1_hi, p2_lo, p2_hi, p3_lo, p3_hi;
        uint64_t p4_lo, p4_hi, p5_lo, p5_hi;

        mul64(a0, rr0, &p0_lo, &p0_hi);
        mul64(a0, rr1, &p1_lo, &p1_hi);
        mul64(a1, rr0, &p2_lo, &p2_hi);
        mul64(a1, rr1, &p3_lo, &p3_hi);
        mul64(a2, rr0, &p4_lo, &p4_hi);
        mul64(a2, rr1, &p5_lo, &p5_hi);

        uint64_t res0, res1, res2, res3, res4;

        res0 = p0_lo;

        res1 = p0_hi;
        carry = 0;
        res1 += p1_lo; carry += (res1 < p1_lo);
        res1 += p2_lo; carry += (res1 < p2_lo);

        res2 = p1_hi + carry;
        carry = (res2 < p1_hi);
        res2 += p2_hi; carry += (res2 < p2_hi);
        res2 += p3_lo; carry += (res2 < p3_lo);
        res2 += p4_lo; carry += (res2 < p4_lo);

        res3 = p3_hi + carry;
        carry = (res3 < p3_hi);
        res3 += p4_hi; carry += (res3 < p4_hi);
        res3 += p5_lo; carry += (res3 < p5_lo);

        res4 = p5_hi + carry;

        uint64_t high0 = (res1 >> 63) | (res2 << 1);
        uint64_t high1 = (res2 >> 63) | (res3 << 1);
        uint64_t high2 = (res3 >> 63) | (res4 << 1);
        uint64_t high3 = res4 >> 63;

        prod1[0] = res0;
        prod1[1] = res1 & 0x7FFFFFFFFFFFFFFFULL;

        carry = 0;
        prod1[0] += high0; carry = (prod1[0] < high0);
        prod1[1] += high1 + carry; carry = (prod1[1] < high1 + carry);
        prod1[2] = high2 + high3 + carry;
    }

    /* prod2 = m2 * r (2-limb × 2-limb) */
    uint64_t t2[2] = {b2_0, b2_1};
    uint64_t prod2[3];
    mul_2limb(prod2, t2, r);

    /* acc = prod1 + prod2 */
    carry = 0;
    acc[0] = prod1[0] + prod2[0];
    carry = (acc[0] < prod1[0]);
    acc[1] = prod1[1] + prod2[1] + carry;
    carry = (acc[1] < prod1[1] + carry);
    acc[2] = prod1[2] + prod2[2] + carry;

    /* always do partial reduction (constant-time) */
    reduce_m127_partial(acc);
}

/*
 * Process 4 blocks at once: acc = (acc + m1)*r⁴ + m2*r³ + m3*r² + m4*r
 * Four independent multiplies for maximum ILP (constant-time).
 */
static void process_4blocks(uint64_t acc[3], const uint8_t* m1, const uint8_t* m2,
                            const uint8_t* m3, const uint8_t* m4,
                            const uint64_t r[2], const uint64_t r2[2],
                            const uint64_t r3[2], const uint64_t r4[2]) {
    /* load all 4 blocks */
    uint64_t b1_0 = load64_le(m1);
    uint64_t b1_1 = (load64_le(m1 + 7) >> 8) | (0x01ULL << 56);
    uint64_t b2_0 = load64_le(m2);
    uint64_t b2_1 = (load64_le(m2 + 7) >> 8) | (0x01ULL << 56);
    uint64_t b3_0 = load64_le(m3);
    uint64_t b3_1 = (load64_le(m3 + 7) >> 8) | (0x01ULL << 56);
    uint64_t b4_0 = load64_le(m4);
    uint64_t b4_1 = (load64_le(m4 + 7) >> 8) | (0x01ULL << 56);

    /* t1 = acc + m1 */
    uint64_t t1[3];
    uint64_t carry;
    t1[0] = acc[0] + b1_0;
    carry = (t1[0] < acc[0]);
    t1[1] = acc[1] + b1_1 + carry;
    carry = (t1[1] < acc[1] + carry);
    t1[2] = acc[2] + carry;

    /* prod1 = t1 * r⁴ (3-limb × 2-limb) */
    uint64_t prod1[3];
    {
        uint64_t a0 = t1[0], a1 = t1[1], a2 = t1[2];
        uint64_t rr0 = r4[0], rr1 = r4[1];

        uint64_t p0_lo, p0_hi, p1_lo, p1_hi, p2_lo, p2_hi, p3_lo, p3_hi;
        uint64_t p4_lo, p4_hi, p5_lo, p5_hi;

        mul64(a0, rr0, &p0_lo, &p0_hi);
        mul64(a0, rr1, &p1_lo, &p1_hi);
        mul64(a1, rr0, &p2_lo, &p2_hi);
        mul64(a1, rr1, &p3_lo, &p3_hi);
        mul64(a2, rr0, &p4_lo, &p4_hi);
        mul64(a2, rr1, &p5_lo, &p5_hi);

        uint64_t res0, res1, res2, res3, res4;

        res0 = p0_lo;

        res1 = p0_hi;
        carry = 0;
        res1 += p1_lo; carry += (res1 < p1_lo);
        res1 += p2_lo; carry += (res1 < p2_lo);

        res2 = p1_hi + carry;
        carry = (res2 < p1_hi);
        res2 += p2_hi; carry += (res2 < p2_hi);
        res2 += p3_lo; carry += (res2 < p3_lo);
        res2 += p4_lo; carry += (res2 < p4_lo);

        res3 = p3_hi + carry;
        carry = (res3 < p3_hi);
        res3 += p4_hi; carry += (res3 < p4_hi);
        res3 += p5_lo; carry += (res3 < p5_lo);

        res4 = p5_hi + carry;

        uint64_t high0 = (res1 >> 63) | (res2 << 1);
        uint64_t high1 = (res2 >> 63) | (res3 << 1);
        uint64_t high2 = (res3 >> 63) | (res4 << 1);
        uint64_t high3 = res4 >> 63;

        prod1[0] = res0;
        prod1[1] = res1 & 0x7FFFFFFFFFFFFFFFULL;

        carry = 0;
        prod1[0] += high0; carry = (prod1[0] < high0);
        prod1[1] += high1 + carry; carry = (prod1[1] < high1 + carry);
        prod1[2] = high2 + high3 + carry;
    }

    /* prod2 = m2 * r³ */
    uint64_t t2[2] = {b2_0, b2_1};
    uint64_t prod2[3];
    mul_2limb(prod2, t2, r3);

    /* prod3 = m3 * r² */
    uint64_t t3[2] = {b3_0, b3_1};
    uint64_t prod3[3];
    mul_2limb(prod3, t3, r2);

    /* prod4 = m4 * r */
    uint64_t t4[2] = {b4_0, b4_1};
    uint64_t prod4[3];
    mul_2limb(prod4, t4, r);

    /* acc = prod1 + prod2 + prod3 + prod4 */
    carry = 0;
    acc[0] = prod1[0] + prod2[0];
    carry = (acc[0] < prod1[0]);
    acc[1] = prod1[1] + prod2[1] + carry;
    carry = (acc[1] < prod1[1] + carry);
    acc[2] = prod1[2] + prod2[2] + carry;

    carry = 0;
    acc[0] += prod3[0];
    carry = (acc[0] < prod3[0]);
    acc[1] += prod3[1] + carry;
    carry = (acc[1] < prod3[1] + carry);
    acc[2] += prod3[2] + carry;

    carry = 0;
    acc[0] += prod4[0];
    carry = (acc[0] < prod4[0]);
    acc[1] += prod4[1] + carry;
    carry = (acc[1] < prod4[1] + carry);
    acc[2] += prod4[2] + carry;

    /* always do partial reduction (constant-time) */
    reduce_m127_partial(acc);
}

/* add full 15-byte block with 0x01 padding (optimized, no memcpy) */
static inline void add_block_full(uint64_t acc[3], const uint8_t* block) {
    uint64_t m0 = load64_le(block);
    /* load bytes 8-14 into low 56 bits, put 0x01 in high byte */
    uint64_t m1 = (load64_le(block + 7) >> 8) | (0x01ULL << 56);

    uint64_t carry;
    acc[0] += m0;
    carry = (acc[0] < m0);
    acc[1] += m1 + carry;
    carry = (acc[1] < m1 + carry);
    acc[2] += carry;
}

/* add partial block with 0x01 padding (uses memcpy) */
static void add_block_partial(uint64_t acc[3], const uint8_t* block, size_t len) {
    uint8_t padded[16] = {0};
    memcpy(padded, block, len);
    padded[len] = 0x01;

    uint64_t m0 = load64_le(padded);
    uint64_t m1 = load64_le(padded + 8);

    uint64_t carry;
    acc[0] += m0;
    carry = (acc[0] < m0);
    acc[1] += m1 + carry;
    carry = (acc[1] < m1 + carry);
    acc[2] += carry;
}

/*
 * Clamp r for side-channel resistance (similar to Poly1305)
 *
 * We clear certain bits to ensure constant-time behavior:
 * - Top bit (bit 127) must be 0 (r < 2^127)
 * - Clear top 4 bits of bytes 3, 7, 11, 15 for alignment
 * - Clear bottom 2 bits of bytes 4, 8, 12 to reduce key space
 *
 * This reduces the effective key space but ensures no timing variations.
 */
static void clamp_r(uint8_t r[16]) {
    r[3]  &= 0x0F;
    r[7]  &= 0x0F;
    r[11] &= 0x0F;
    r[15] &= 0x07;  /* also clears bit 127 */

    r[4]  &= 0xFC;
    r[8]  &= 0xFC;
    r[12] &= 0xFC;
}

/* reduce fully every N blocks (tunable) */
#define LAZY_REDUCE_INTERVAL 8

/* multiply two 2-limb values mod M127, result in out[2] */
static void mul_mod_m127(uint64_t out[2], const uint64_t a[2], const uint64_t b[2]) {
    uint64_t a0 = a[0], a1 = a[1];
    uint64_t b0 = b[0], b1 = b[1];

    uint64_t p0_lo, p0_hi, p1_lo, p1_hi, p2_lo, p2_hi, p3_lo, p3_hi;

    mul64(a0, b0, &p0_lo, &p0_hi);
    mul64(a0, b1, &p1_lo, &p1_hi);
    mul64(a1, b0, &p2_lo, &p2_hi);
    mul64(a1, b1, &p3_lo, &p3_hi);

    uint64_t res0, res1, res2, res3;
    uint64_t carry;

    res0 = p0_lo;

    res1 = p0_hi;
    carry = 0;
    res1 += p1_lo; carry += (res1 < p1_lo);
    res1 += p2_lo; carry += (res1 < p2_lo);

    res2 = p1_hi + carry;
    carry = (res2 < p1_hi);
    res2 += p2_hi; carry += (res2 < p2_hi);
    res2 += p3_lo; carry += (res2 < p3_lo);

    res3 = p3_hi + carry;

    /* fold bits [127..] back */
    uint64_t high0 = (res1 >> 63) | (res2 << 1);
    uint64_t high1 = (res2 >> 63) | (res3 << 1);

    out[0] = res0;
    out[1] = res1 & 0x7FFFFFFFFFFFFFFFULL;

    carry = 0;
    out[0] += high0; carry = (out[0] < high0);
    out[1] += high1 + carry; carry = (out[1] < high1 + carry);

    /* final reduction */
    uint64_t overflow = (out[1] >> 63) + carry;
    out[1] &= 0x7FFFFFFFFFFFFFFFULL;
    uint64_t old_out0 = out[0];
    out[0] += overflow;
    uint64_t carry2 = (out[0] < old_out0);
    out[1] += carry2;

    /* second pass if needed */
    overflow = out[1] >> 63;
    out[1] &= 0x7FFFFFFFFFFFFFFFULL;
    out[0] += overflow;
}

/* compute r² mod M127 for 2-block processing */
static void square_r(uint64_t r2[2], const uint64_t r[2]) {
    uint64_t r0 = r[0], r1 = r[1];

    uint64_t p0_lo, p0_hi, p1_lo, p1_hi, p2_lo, p2_hi;

    /* r² = r0² + 2*r0*r1*2^64 + r1²*2^128 */
    mul64(r0, r0, &p0_lo, &p0_hi);  /* r0² */
    mul64(r0, r1, &p1_lo, &p1_hi);  /* r0*r1 (will be doubled) */
    mul64(r1, r1, &p2_lo, &p2_hi);  /* r1² */

    /* double the cross term */
    uint64_t p1_lo2 = p1_lo << 1;
    uint64_t p1_hi2 = (p1_hi << 1) | (p1_lo >> 63);

    /* sum into 256-bit result */
    uint64_t res0, res1, res2, res3;
    uint64_t carry;

    res0 = p0_lo;

    res1 = p0_hi;
    carry = 0;
    res1 += p1_lo2; carry = (res1 < p1_lo2);

    res2 = p1_hi2 + carry;
    carry = (res2 < p1_hi2);
    res2 += p2_lo; carry += (res2 < p2_lo);

    res3 = p2_hi + carry;

    /* reduce mod M127: fold bits [127..] back */
    uint64_t high0 = (res1 >> 63) | (res2 << 1);
    uint64_t high1 = (res2 >> 63) | (res3 << 1);

    r2[0] = res0;
    r2[1] = res1 & 0x7FFFFFFFFFFFFFFFULL;

    carry = 0;
    r2[0] += high0; carry = (r2[0] < high0);
    r2[1] += high1 + carry; carry = (r2[1] < high1 + carry);

    /* note: res3 >> 63 is always 0 for r < 2^127 squared, so we skip high2 */

    /* final reduction (constant-time) */
    uint64_t overflow = (r2[1] >> 63) + carry;
    r2[1] &= 0x7FFFFFFFFFFFFFFFULL;
    uint64_t old_r2_0 = r2[0];
    r2[0] += overflow;
    uint64_t carry2 = (r2[0] < old_r2_0);
    r2[1] += carry2;

    /* second pass if needed (constant-time) */
    overflow = r2[1] >> 63;
    r2[1] &= 0x7FFFFFFFFFFFFFFFULL;
    r2[0] += overflow;
}

void cfx_poly1271_init(cfx_poly1271_ctx_t* ctx, const uint8_t key[32]) {
    uint8_t r_clamped[16];
    memcpy(r_clamped, key, 16);
    clamp_r(r_clamped);

    ctx->r[0] = load64_le(r_clamped);
    ctx->r[1] = load64_le(r_clamped + 8);

    /* precompute r², r³, r⁴ for multi-block processing */
    square_r(ctx->r2, ctx->r);
    mul_mod_m127(ctx->r3, ctx->r2, ctx->r);   /* r³ = r² * r */
    mul_mod_m127(ctx->r4, ctx->r2, ctx->r2);  /* r⁴ = r² * r² */

    ctx->s[0] = load64_le(key + 16);
    ctx->s[1] = load64_le(key + 24);

    ctx->acc[0] = 0;
    ctx->acc[1] = 0;
    ctx->acc[2] = 0;

    ctx->buflen = 0;
    ctx->blkcnt = 0;

    CFX_MEMZERO_S(r_clamped, sizeof(r_clamped));
}

void cfx_poly1271_update(cfx_poly1271_ctx_t* ctx, const uint8_t* msg, size_t len) {
    /* handle buffered data first */
    if (ctx->buflen > 0) {
        size_t need = CFX_POLY1271_BLOCK_SIZE - ctx->buflen;
        if (len < need) {
            memcpy(ctx->buf + ctx->buflen, msg, len);
            ctx->buflen += (uint8_t)len;
            return;
        }
        memcpy(ctx->buf + ctx->buflen, msg, need);
        add_block_full(ctx->acc, ctx->buf);

        ctx->blkcnt++;
        if (ctx->blkcnt >= LAZY_REDUCE_INTERVAL) {
            mul_r_full(ctx->acc, ctx->r);
            ctx->blkcnt = 0;
        } else {
            mul_r_lazy(ctx->acc, ctx->r);
        }

        msg += need;
        len -= need;
        ctx->buflen = 0;
    }

    /* process 4 blocks at a time using r⁴ */
    while (len >= 4 * CFX_POLY1271_BLOCK_SIZE) {
        process_4blocks(ctx->acc, msg,
                        msg + CFX_POLY1271_BLOCK_SIZE,
                        msg + 2 * CFX_POLY1271_BLOCK_SIZE,
                        msg + 3 * CFX_POLY1271_BLOCK_SIZE,
                        ctx->r, ctx->r2, ctx->r3, ctx->r4);
        ctx->blkcnt += 4;

        /* full reduction periodically */
        if (ctx->blkcnt >= LAZY_REDUCE_INTERVAL) {
            reduce_m127_full(ctx->acc);
            ctx->blkcnt = 0;
        }

        msg += 4 * CFX_POLY1271_BLOCK_SIZE;
        len -= 4 * CFX_POLY1271_BLOCK_SIZE;
    }

    /* process remaining 2 blocks using r² */
    if (len >= 2 * CFX_POLY1271_BLOCK_SIZE) {
        process_2blocks(ctx->acc, msg, msg + CFX_POLY1271_BLOCK_SIZE,
                        ctx->r, ctx->r2);
        ctx->blkcnt += 2;

        if (ctx->blkcnt >= LAZY_REDUCE_INTERVAL) {
            reduce_m127_full(ctx->acc);
            ctx->blkcnt = 0;
        }

        msg += 2 * CFX_POLY1271_BLOCK_SIZE;
        len -= 2 * CFX_POLY1271_BLOCK_SIZE;
    }

    /* handle remaining single block */
    if (len >= CFX_POLY1271_BLOCK_SIZE) {
        add_block_full(ctx->acc, msg);

        ctx->blkcnt++;
        if (ctx->blkcnt >= LAZY_REDUCE_INTERVAL) {
            mul_r_full(ctx->acc, ctx->r);
            ctx->blkcnt = 0;
        } else {
            mul_r_lazy(ctx->acc, ctx->r);
        }

        msg += CFX_POLY1271_BLOCK_SIZE;
        len -= CFX_POLY1271_BLOCK_SIZE;
    }

    /* buffer any remaining bytes */
    if (len > 0) {
        memcpy(ctx->buf, msg, len);
        ctx->buflen = (uint8_t)len;
    }
}

void cfx_poly1271_finish(cfx_poly1271_ctx_t* ctx, uint8_t tag[16]) {
    if (ctx->buflen > 0) {
        add_block_partial(ctx->acc, ctx->buf, ctx->buflen);
        mul_r_full(ctx->acc, ctx->r);
    } else {
        /* ensure full reduction even if no final block */
        reduce_m127_full(ctx->acc);
    }

    /* fully reduce: if acc >= M127, subtract M127 (constant-time) */
    /* acc >= M127 iff (acc[1] > M127_HI) || (acc[1] == M127_HI && acc[0] == M127_LO) */
    /* M127 = (M127_HI, M127_LO) = (0x7FFF..., 0xFFFF...) */
    /* Since acc[1] <= M127_HI after full reduction, we only need to check equality */
    uint64_t eq_hi = ~(ctx->acc[1] ^ M127_HI);  /* all 1s if equal */
    eq_hi &= eq_hi >> 32; eq_hi &= eq_hi >> 16; eq_hi &= eq_hi >> 8;
    eq_hi &= eq_hi >> 4; eq_hi &= eq_hi >> 2; eq_hi &= eq_hi >> 1;
    uint64_t hi_eq_mask = 0 - (eq_hi & 1);  /* all 1s if acc[1] == M127_HI */

    uint64_t eq_lo = ~(ctx->acc[0] ^ M127_LO);
    eq_lo &= eq_lo >> 32; eq_lo &= eq_lo >> 16; eq_lo &= eq_lo >> 8;
    eq_lo &= eq_lo >> 4; eq_lo &= eq_lo >> 2; eq_lo &= eq_lo >> 1;
    uint64_t lo_eq_mask = 0 - (eq_lo & 1);  /* all 1s if acc[0] == M127_LO */

    /* need_reduce if both are equal (acc == M127) */
    uint64_t need_reduce = hi_eq_mask & lo_eq_mask & 1;

    /* conditionally add 1 (which wraps M127 to 0) */
    uint64_t old_acc0 = ctx->acc[0];
    ctx->acc[0] += need_reduce;
    uint64_t carry = (ctx->acc[0] < old_acc0);
    ctx->acc[1] += carry;
    ctx->acc[1] &= 0x7FFFFFFFFFFFFFFFULL;

    /* tag = (acc + s) mod 2^128 */
    uint64_t t0 = ctx->acc[0] + ctx->s[0];
    carry = (t0 < ctx->acc[0]);
    uint64_t t1 = ctx->acc[1] + ctx->s[1] + carry;

    store64_le(tag, t0);
    store64_le(tag + 8, t1);

    CFX_MEMZERO_S(ctx, sizeof(*ctx));
}

void cfx_poly1271(uint8_t tag[16], const uint8_t* msg, size_t len, const uint8_t key[32]) {
    /* dedicated one-shot: skip buffer management overhead */
    uint8_t r_clamped[16];
    memcpy(r_clamped, key, 16);
    clamp_r(r_clamped);

    uint64_t r[2], r2[2], r3[2], r4[2], s[2];
    r[0] = load64_le(r_clamped);
    r[1] = load64_le(r_clamped + 8);

    square_r(r2, r);
    mul_mod_m127(r3, r2, r);
    mul_mod_m127(r4, r2, r2);

    s[0] = load64_le(key + 16);
    s[1] = load64_le(key + 24);

    uint64_t acc[3] = {0, 0, 0};
    uint8_t blkcnt = 0;

    /* process 4 blocks at a time */
    while (len >= 4 * CFX_POLY1271_BLOCK_SIZE) {
        process_4blocks(acc, msg,
                        msg + CFX_POLY1271_BLOCK_SIZE,
                        msg + 2 * CFX_POLY1271_BLOCK_SIZE,
                        msg + 3 * CFX_POLY1271_BLOCK_SIZE,
                        r, r2, r3, r4);
        blkcnt += 4;
        if (blkcnt >= LAZY_REDUCE_INTERVAL) {
            reduce_m127_full(acc);
            blkcnt = 0;
        }
        msg += 4 * CFX_POLY1271_BLOCK_SIZE;
        len -= 4 * CFX_POLY1271_BLOCK_SIZE;
    }

    /* process remaining 2 blocks */
    if (len >= 2 * CFX_POLY1271_BLOCK_SIZE) {
        process_2blocks(acc, msg, msg + CFX_POLY1271_BLOCK_SIZE, r, r2);
        blkcnt += 2;
        if (blkcnt >= LAZY_REDUCE_INTERVAL) {
            reduce_m127_full(acc);
            blkcnt = 0;
        }
        msg += 2 * CFX_POLY1271_BLOCK_SIZE;
        len -= 2 * CFX_POLY1271_BLOCK_SIZE;
    }

    /* process remaining single block */
    if (len >= CFX_POLY1271_BLOCK_SIZE) {
        add_block_full(acc, msg);
        mul_r_lazy(acc, r);
        msg += CFX_POLY1271_BLOCK_SIZE;
        len -= CFX_POLY1271_BLOCK_SIZE;
    }

    /* final partial block */
    if (len > 0) {
        add_block_partial(acc, msg, len);
        mul_r_full(acc, r);
    } else {
        reduce_m127_full(acc);
    }

    /* fully reduce: if acc == M127, wrap to 0 (constant-time) */
    uint64_t eq_hi = ~(acc[1] ^ M127_HI);
    eq_hi &= eq_hi >> 32; eq_hi &= eq_hi >> 16; eq_hi &= eq_hi >> 8;
    eq_hi &= eq_hi >> 4; eq_hi &= eq_hi >> 2; eq_hi &= eq_hi >> 1;
    uint64_t hi_eq_mask = 0 - (eq_hi & 1);

    uint64_t eq_lo = ~(acc[0] ^ M127_LO);
    eq_lo &= eq_lo >> 32; eq_lo &= eq_lo >> 16; eq_lo &= eq_lo >> 8;
    eq_lo &= eq_lo >> 4; eq_lo &= eq_lo >> 2; eq_lo &= eq_lo >> 1;
    uint64_t lo_eq_mask = 0 - (eq_lo & 1);

    uint64_t need_reduce = hi_eq_mask & lo_eq_mask & 1;
    uint64_t old_acc0 = acc[0];
    acc[0] += need_reduce;
    uint64_t carry = (acc[0] < old_acc0);
    acc[1] += carry;
    acc[1] &= 0x7FFFFFFFFFFFFFFFULL;

    /* tag = (acc + s) mod 2^128 */
    uint64_t t0 = acc[0] + s[0];
    carry = (t0 < acc[0]);
    uint64_t t1 = acc[1] + s[1] + carry;

    store64_le(tag, t0);
    store64_le(tag + 8, t1);

    CFX_MEMZERO_S(r_clamped, sizeof(r_clamped));
    CFX_MEMZERO_S(r, sizeof(r));
    CFX_MEMZERO_S(r2, sizeof(r2));
    CFX_MEMZERO_S(r3, sizeof(r3));
    CFX_MEMZERO_S(r4, sizeof(r4));
    CFX_MEMZERO_S(s, sizeof(s));
    CFX_MEMZERO_S(acc, sizeof(acc));
}

/* ========================================================================== */
/* AVX2 implementation - radix-2^26, 4-way parallel                           */
/* ========================================================================== */

#if CFX_HAVE_AVX2

#include "cfx/arch.h"

#define AVX2_MASK26 ((1ULL << 26) - 1)
#define AVX2_MASK23 ((1ULL << 23) - 1)

/* radix-26 conversions */

static inline void avx2_decode_block_to_limbs(uint64_t limbs[5], const uint8_t* block) {
    uint64_t lo = load64_le(block);
    uint64_t hi = (uint64_t)block[8] | ((uint64_t)block[9] << 8) |
                  ((uint64_t)block[10] << 16) | ((uint64_t)block[11] << 24) |
                  ((uint64_t)block[12] << 32) | ((uint64_t)block[13] << 40) |
                  ((uint64_t)block[14] << 48) | (1ULL << 56); /* delimiter */

    limbs[0] = lo & AVX2_MASK26;
    limbs[1] = (lo >> 26) & AVX2_MASK26;
    limbs[2] = ((lo >> 52) | (hi << 12)) & AVX2_MASK26;
    limbs[3] = (hi >> 14) & AVX2_MASK26;
    limbs[4] = (hi >> 40) & AVX2_MASK26;
}

/* decode 4 consecutive 15-byte blocks directly to vector format using gather */
static inline void avx2_decode_4blocks_to_vec(__m256i out[5], const uint8_t* msg) {
    const __m256i lo_idx = _mm256_set_epi64x(45, 30, 15, 0);
    const __m256i hi_idx = _mm256_set_epi64x(53, 38, 23, 8);

    __m256i lo = _mm256_i64gather_epi64((const long long*)msg, lo_idx, 1);
    __m256i hi = _mm256_i64gather_epi64((const long long*)msg, hi_idx, 1);

    const __m256i mask56 = _mm256_set1_epi64x(0x00FFFFFFFFFFFFFFULL);
    const __m256i delim = _mm256_set1_epi64x(1ULL << 56);
    hi = _mm256_or_si256(_mm256_and_si256(hi, mask56), delim);

    const __m256i mask26 = _mm256_set1_epi64x(AVX2_MASK26);

    out[0] = _mm256_and_si256(lo, mask26);
    out[1] = _mm256_and_si256(_mm256_srli_epi64(lo, 26), mask26);
    __m256i lo_hi = _mm256_srli_epi64(lo, 52);
    __m256i hi_lo = _mm256_slli_epi64(hi, 12);
    out[2] = _mm256_and_si256(_mm256_or_si256(lo_hi, hi_lo), mask26);
    out[3] = _mm256_and_si256(_mm256_srli_epi64(hi, 14), mask26);
    out[4] = _mm256_and_si256(_mm256_srli_epi64(hi, 40), mask26);
}

static inline void avx2_decode_partial_to_limbs(uint64_t limbs[5], const uint8_t* block, size_t len) {
    uint8_t pad[16] = {0};
    memcpy(pad, block, len);
    pad[len] = 0x01;

    uint64_t lo = load64_le(pad);
    uint64_t hi = load64_le(pad + 8);

    limbs[0] = lo & AVX2_MASK26;
    limbs[1] = (lo >> 26) & AVX2_MASK26;
    limbs[2] = ((lo >> 52) | (hi << 12)) & AVX2_MASK26;
    limbs[3] = (hi >> 14) & AVX2_MASK26;
    limbs[4] = (hi >> 40) & AVX2_MASK26;
}

static inline void avx2_convert_to_limbs(uint64_t limbs[5], const uint64_t v[2]) {
    uint64_t lo = v[0], hi = v[1];
    limbs[0] = lo & AVX2_MASK26;
    limbs[1] = (lo >> 26) & AVX2_MASK26;
    limbs[2] = ((lo >> 52) | (hi << 12)) & AVX2_MASK26;
    limbs[3] = (hi >> 14) & AVX2_MASK26;
    limbs[4] = (hi >> 40) & AVX2_MASK26;
}

static inline void avx2_convert_from_limbs(uint64_t v[2], const uint64_t limbs[5]) {
    v[0] = limbs[0] | (limbs[1] << 26) | (limbs[2] << 52);
    v[1] = (limbs[2] >> 12) | (limbs[3] << 14) | (limbs[4] << 40);
}

/* scalar radix-26 arithmetic */

static inline void avx2_carry26_scalar(uint64_t limbs[5]) {
    uint64_t c;
    c = limbs[0] >> 26; limbs[0] &= AVX2_MASK26; limbs[1] += c;
    c = limbs[1] >> 26; limbs[1] &= AVX2_MASK26; limbs[2] += c;
    c = limbs[2] >> 26; limbs[2] &= AVX2_MASK26; limbs[3] += c;
    c = limbs[3] >> 26; limbs[3] &= AVX2_MASK26; limbs[4] += c;
}

static inline void avx2_reduce_p_scalar(uint64_t limbs[5]) {
    uint64_t c = limbs[4] >> 23;
    limbs[4] &= AVX2_MASK23;
    limbs[0] += c;
    avx2_carry26_scalar(limbs);
    c = limbs[4] >> 23;
    limbs[4] &= AVX2_MASK23;
    limbs[0] += c;
    avx2_carry26_scalar(limbs);
}

static void avx2_mul_limbs_scalar(uint64_t out[5], const uint64_t a[5], const uint64_t b[5]) {
    uint64_t d0 = a[0] * b[0];
    uint64_t d1 = a[0] * b[1] + a[1] * b[0];
    uint64_t d2 = a[0] * b[2] + a[1] * b[1] + a[2] * b[0];
    uint64_t d3 = a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0];
    uint64_t d4 = a[0] * b[4] + a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + a[4] * b[0];
    uint64_t d5 = a[1] * b[4] + a[2] * b[3] + a[3] * b[2] + a[4] * b[1];
    uint64_t d6 = a[2] * b[4] + a[3] * b[3] + a[4] * b[2];
    uint64_t d7 = a[3] * b[4] + a[4] * b[3];
    uint64_t d8 = a[4] * b[4];

    d0 += d5 << 3;
    d1 += d6 << 3;
    d2 += d7 << 3;
    d3 += d8 << 3;

    uint64_t c;
    c = d0 >> 26; d0 &= AVX2_MASK26; d1 += c;
    c = d1 >> 26; d1 &= AVX2_MASK26; d2 += c;
    c = d2 >> 26; d2 &= AVX2_MASK26; d3 += c;
    c = d3 >> 26; d3 &= AVX2_MASK26; d4 += c;

    out[0] = d0; out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4;
    avx2_reduce_p_scalar(out);
}

static void avx2_square_limbs_scalar(uint64_t out[5], const uint64_t a[5]) {
    uint64_t d0 = a[0] * a[0];
    uint64_t d1 = 2 * a[0] * a[1];
    uint64_t d2 = 2 * a[0] * a[2] + a[1] * a[1];
    uint64_t d3 = 2 * a[0] * a[3] + 2 * a[1] * a[2];
    uint64_t d4 = 2 * a[0] * a[4] + 2 * a[1] * a[3] + a[2] * a[2];
    uint64_t d5 = 2 * a[1] * a[4] + 2 * a[2] * a[3];
    uint64_t d6 = 2 * a[2] * a[4] + a[3] * a[3];
    uint64_t d7 = 2 * a[3] * a[4];
    uint64_t d8 = a[4] * a[4];

    d0 += d5 << 3;
    d1 += d6 << 3;
    d2 += d7 << 3;
    d3 += d8 << 3;

    uint64_t c;
    c = d0 >> 26; d0 &= AVX2_MASK26; d1 += c;
    c = d1 >> 26; d1 &= AVX2_MASK26; d2 += c;
    c = d2 >> 26; d2 &= AVX2_MASK26; d3 += c;
    c = d3 >> 26; d3 &= AVX2_MASK26; d4 += c;

    out[0] = d0; out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4;
    avx2_reduce_p_scalar(out);
}

/* multiply 4 parallel accumulators by r^4 */
static inline void avx2_mul_r4(__m256i acc[5], const __m256i r4v[5]) {
    const __m256i r0 = r4v[0], r1 = r4v[1], r2 = r4v[2], r3 = r4v[3], r4 = r4v[4];
    __m256i a0 = acc[0], a1 = acc[1], a2 = acc[2], a3 = acc[3], a4 = acc[4];

    __m256i d0 = _mm256_mul_epu32(a0, r0);
    __m256i d1 = _mm256_add_epi64(_mm256_mul_epu32(a0, r1), _mm256_mul_epu32(a1, r0));
    __m256i d2 = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(a0, r2), _mm256_mul_epu32(a1, r1)),
                                  _mm256_mul_epu32(a2, r0));
    __m256i d3 = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(a0, r3), _mm256_mul_epu32(a1, r2)),
                                  _mm256_add_epi64(_mm256_mul_epu32(a2, r1), _mm256_mul_epu32(a3, r0)));
    __m256i d4 = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(a0, r4), _mm256_mul_epu32(a1, r3)),
                                  _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(a2, r2), _mm256_mul_epu32(a3, r1)),
                                                   _mm256_mul_epu32(a4, r0)));
    __m256i d5 = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(a1, r4), _mm256_mul_epu32(a2, r3)),
                                  _mm256_add_epi64(_mm256_mul_epu32(a3, r2), _mm256_mul_epu32(a4, r1)));
    __m256i d6 = _mm256_add_epi64(_mm256_add_epi64(_mm256_mul_epu32(a2, r4), _mm256_mul_epu32(a3, r3)),
                                  _mm256_mul_epu32(a4, r2));
    __m256i d7 = _mm256_add_epi64(_mm256_mul_epu32(a3, r4), _mm256_mul_epu32(a4, r3));
    __m256i d8 = _mm256_mul_epu32(a4, r4);

    d0 = _mm256_add_epi64(d0, _mm256_slli_epi64(d5, 3));
    d1 = _mm256_add_epi64(d1, _mm256_slli_epi64(d6, 3));
    d2 = _mm256_add_epi64(d2, _mm256_slli_epi64(d7, 3));
    d3 = _mm256_add_epi64(d3, _mm256_slli_epi64(d8, 3));

    __m256i mask26 = _mm256_set1_epi64x(AVX2_MASK26);
    __m256i mask23 = _mm256_set1_epi64x(AVX2_MASK23);
    __m256i c;

    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);
    c = _mm256_srli_epi64(d1, 26); d1 = _mm256_and_si256(d1, mask26); d2 = _mm256_add_epi64(d2, c);
    c = _mm256_srli_epi64(d2, 26); d2 = _mm256_and_si256(d2, mask26); d3 = _mm256_add_epi64(d3, c);
    c = _mm256_srli_epi64(d3, 26); d3 = _mm256_and_si256(d3, mask26); d4 = _mm256_add_epi64(d4, c);
    c = _mm256_srli_epi64(d4, 23); d4 = _mm256_and_si256(d4, mask23); d0 = _mm256_add_epi64(d0, c);
    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);

    acc[0] = d0; acc[1] = d1; acc[2] = d2; acc[3] = d3; acc[4] = d4;
}

static inline void avx2_add_4blocks_vec(__m256i acc[5], const __m256i blocks[5]) {
    acc[0] = _mm256_add_epi64(acc[0], blocks[0]);
    acc[1] = _mm256_add_epi64(acc[1], blocks[1]);
    acc[2] = _mm256_add_epi64(acc[2], blocks[2]);
    acc[3] = _mm256_add_epi64(acc[3], blocks[3]);
    acc[4] = _mm256_add_epi64(acc[4], blocks[4]);
}

static inline uint64_t avx2_hsum_epi64(__m256i v) {
    __m256i hi = _mm256_permute4x64_epi64(v, 0x4E);
    __m256i sum1 = _mm256_add_epi64(v, hi);
    __m128i lo128 = _mm256_castsi256_si128(sum1);
    __m128i hi128 = _mm_unpackhi_epi64(lo128, lo128);
    __m128i sum2 = _mm_add_epi64(lo128, hi128);
    return (uint64_t)_mm_cvtsi128_si64(sum2);
}

static void avx2_combine_4lanes(uint64_t out[5], const __m256i acc[5], const __m256i rv[5]) {
    const __m256i a0 = acc[0], a1 = acc[1], a2 = acc[2], a3 = acc[3], a4 = acc[4];
    const __m256i r0 = rv[0], r1 = rv[1], r2 = rv[2], r3 = rv[3], r4 = rv[4];

    __m256i d0 = _mm256_mul_epu32(a0, r0);
    __m256i d1 = _mm256_mul_epu32(a0, r1);
    __m256i d2 = _mm256_mul_epu32(a0, r2);
    __m256i d3 = _mm256_mul_epu32(a0, r3);
    __m256i d4 = _mm256_mul_epu32(a0, r4);
    __m256i d5 = _mm256_mul_epu32(a1, r4);
    __m256i d6 = _mm256_mul_epu32(a2, r4);
    __m256i d7 = _mm256_mul_epu32(a3, r4);
    __m256i d8 = _mm256_mul_epu32(a4, r4);

    d1 = _mm256_add_epi64(d1, _mm256_mul_epu32(a1, r0));
    d2 = _mm256_add_epi64(d2, _mm256_mul_epu32(a1, r1));
    d3 = _mm256_add_epi64(d3, _mm256_mul_epu32(a1, r2));
    d4 = _mm256_add_epi64(d4, _mm256_mul_epu32(a1, r3));
    d5 = _mm256_add_epi64(d5, _mm256_mul_epu32(a2, r3));
    d6 = _mm256_add_epi64(d6, _mm256_mul_epu32(a3, r3));
    d7 = _mm256_add_epi64(d7, _mm256_mul_epu32(a4, r3));

    d2 = _mm256_add_epi64(d2, _mm256_mul_epu32(a2, r0));
    d3 = _mm256_add_epi64(d3, _mm256_mul_epu32(a2, r1));
    d4 = _mm256_add_epi64(d4, _mm256_mul_epu32(a2, r2));
    d5 = _mm256_add_epi64(d5, _mm256_mul_epu32(a3, r2));
    d6 = _mm256_add_epi64(d6, _mm256_mul_epu32(a4, r2));

    d3 = _mm256_add_epi64(d3, _mm256_mul_epu32(a3, r0));
    d4 = _mm256_add_epi64(d4, _mm256_mul_epu32(a3, r1));
    d5 = _mm256_add_epi64(d5, _mm256_mul_epu32(a4, r1));

    d4 = _mm256_add_epi64(d4, _mm256_mul_epu32(a4, r0));

    d0 = _mm256_add_epi64(d0, _mm256_slli_epi64(d5, 3));
    d1 = _mm256_add_epi64(d1, _mm256_slli_epi64(d6, 3));
    d2 = _mm256_add_epi64(d2, _mm256_slli_epi64(d7, 3));
    d3 = _mm256_add_epi64(d3, _mm256_slli_epi64(d8, 3));

    __m256i mask26 = _mm256_set1_epi64x(AVX2_MASK26);
    __m256i mask23 = _mm256_set1_epi64x(AVX2_MASK23);
    __m256i c;

    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);
    c = _mm256_srli_epi64(d1, 26); d1 = _mm256_and_si256(d1, mask26); d2 = _mm256_add_epi64(d2, c);
    c = _mm256_srli_epi64(d2, 26); d2 = _mm256_and_si256(d2, mask26); d3 = _mm256_add_epi64(d3, c);
    c = _mm256_srli_epi64(d3, 26); d3 = _mm256_and_si256(d3, mask26); d4 = _mm256_add_epi64(d4, c);
    c = _mm256_srli_epi64(d4, 23); d4 = _mm256_and_si256(d4, mask23); d0 = _mm256_add_epi64(d0, c);
    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);

    out[0] = avx2_hsum_epi64(d0);
    out[1] = avx2_hsum_epi64(d1);
    out[2] = avx2_hsum_epi64(d2);
    out[3] = avx2_hsum_epi64(d3);
    out[4] = avx2_hsum_epi64(d4);

    avx2_carry26_scalar(out);
    avx2_reduce_p_scalar(out);
}

static void avx2_process_nblocks(uint64_t acc[5], const uint8_t* msg, int nrounds,
                                  const __m256i r4v[5], const __m256i rv[5]) {
    __m256i pacc[5];
    for (int i = 0; i < 5; i++)
        pacc[i] = _mm256_set_epi64x(0, 0, 0, (int64_t)acc[i]);

    for (int round = 0; round < nrounds; round++) {
        __m256i blocks[5];
        avx2_decode_4blocks_to_vec(blocks, msg + round * 60);
        avx2_add_4blocks_vec(pacc, blocks);

        if (round < nrounds - 1)
            avx2_mul_r4(pacc, r4v);
    }

    avx2_combine_4lanes(acc, pacc, rv);
}

static void avx2_process_1block_scalar(uint64_t acc[5], const uint8_t* m, const uint64_t r[5]) {
    uint64_t b[5];
    avx2_decode_block_to_limbs(b, m);
    for (int i = 0; i < 5; i++) acc[i] += b[i];
    avx2_carry26_scalar(acc);
    avx2_mul_limbs_scalar(acc, acc, r);
}

static void avx2_finalize(uint64_t acc[5], const uint64_t s[2], uint8_t tag[16]) {
    avx2_reduce_p_scalar(acc);

    uint64_t v[2];
    avx2_convert_from_limbs(v, acc);

    uint64_t eq_lo = ~(v[0] ^ M127_LO);
    eq_lo &= eq_lo >> 32; eq_lo &= eq_lo >> 16; eq_lo &= eq_lo >> 8;
    eq_lo &= eq_lo >> 4;  eq_lo &= eq_lo >> 2;  eq_lo &= eq_lo >> 1;

    uint64_t eq_hi = ~(v[1] ^ M127_HI);
    eq_hi &= eq_hi >> 32; eq_hi &= eq_hi >> 16; eq_hi &= eq_hi >> 8;
    eq_hi &= eq_hi >> 4;  eq_hi &= eq_hi >> 2;  eq_hi &= eq_hi >> 1;

    uint64_t wrap = (uint64_t)0 - (eq_lo & eq_hi & 1);
    v[0] += wrap & 1;
    uint64_t c = (v[0] == 0) & (wrap & 1);
    v[1] = (v[1] + c) & M127_HI;

    uint64_t t0 = v[0] + s[0];
    c = (t0 < v[0]);
    uint64_t t1 = v[1] + s[1] + c;

    store64_le(tag, t0);
    store64_le(tag + 8, t1);
}

/* public API */

void cfx_poly1271_avx2_init(cfx_poly1271_avx2_ctx_t* ctx, const uint8_t key[32]) {
    uint8_t rc[16];
    memcpy(rc, key, 16);
    clamp_r(rc);

    uint64_t r64[2] = {load64_le(rc), load64_le(rc + 8)};
    avx2_convert_to_limbs(ctx->r, r64);

    avx2_square_limbs_scalar(ctx->r2, ctx->r);
    avx2_mul_limbs_scalar(ctx->r3, ctx->r2, ctx->r);
    avx2_mul_limbs_scalar(ctx->r4, ctx->r2, ctx->r2);

    for (int i = 0; i < 5; i++) {
        ctx->rv[i] = _mm256_set_epi64x((int64_t)ctx->r[i], (int64_t)ctx->r2[i],
                                        (int64_t)ctx->r3[i], (int64_t)ctx->r4[i]);
        ctx->r4v[i] = _mm256_set1_epi64x((int64_t)ctx->r4[i]);
    }

    ctx->s[0] = load64_le(key + 16);
    ctx->s[1] = load64_le(key + 24);

    for (int i = 0; i < 5; i++) ctx->acc[i] = 0;
    ctx->buflen = 0;

    CFX_MEMZERO_S(rc, sizeof(rc));
}

void cfx_poly1271_avx2_update(cfx_poly1271_avx2_ctx_t* ctx, const uint8_t* msg, size_t len) {
    if (ctx->buflen) {
        size_t need = 15 - ctx->buflen;
        if (len < need) {
            memcpy(ctx->buf + ctx->buflen, msg, len);
            ctx->buflen += (uint8_t)len;
            return;
        }
        memcpy(ctx->buf + ctx->buflen, msg, need);
        avx2_process_1block_scalar(ctx->acc, ctx->buf, ctx->r);
        msg += need;
        len -= need;
        ctx->buflen = 0;
    }

    while (len >= 1920) {
        avx2_process_nblocks(ctx->acc, msg, 32, ctx->r4v, ctx->rv);
        msg += 1920; len -= 1920;
    }
    while (len >= 960) {
        avx2_process_nblocks(ctx->acc, msg, 16, ctx->r4v, ctx->rv);
        msg += 960; len -= 960;
    }
    while (len >= 480) {
        avx2_process_nblocks(ctx->acc, msg, 8, ctx->r4v, ctx->rv);
        msg += 480; len -= 480;
    }
    while (len >= 240) {
        avx2_process_nblocks(ctx->acc, msg, 4, ctx->r4v, ctx->rv);
        msg += 240; len -= 240;
    }
    while (len >= 60) {
        avx2_process_nblocks(ctx->acc, msg, 1, ctx->r4v, ctx->rv);
        msg += 60; len -= 60;
    }

    while (len >= 15) {
        avx2_process_1block_scalar(ctx->acc, msg, ctx->r);
        msg += 15;
        len -= 15;
    }

    if (len) {
        memcpy(ctx->buf, msg, len);
        ctx->buflen = (uint8_t)len;
    }
}

void cfx_poly1271_avx2_finish(cfx_poly1271_avx2_ctx_t* ctx, uint8_t tag[16]) {
    if (ctx->buflen) {
        uint64_t b[5];
        avx2_decode_partial_to_limbs(b, ctx->buf, ctx->buflen);
        for (int i = 0; i < 5; i++) ctx->acc[i] += b[i];
        avx2_carry26_scalar(ctx->acc);
        avx2_mul_limbs_scalar(ctx->acc, ctx->acc, ctx->r);
    }
    avx2_finalize(ctx->acc, ctx->s, tag);
    CFX_MEMZERO_S(ctx, sizeof(*ctx));
}

void cfx_poly1271_avx2(uint8_t tag[16], const uint8_t* msg, size_t len, const uint8_t key[32]) {
    cfx_poly1271_avx2_ctx_t ctx;
    cfx_poly1271_avx2_init(&ctx, key);
    cfx_poly1271_avx2_update(&ctx, msg, len);
    cfx_poly1271_avx2_finish(&ctx, tag);
}

int cfx_poly1271_avx2_verify(const uint8_t tag[16], const uint8_t* msg, size_t len, const uint8_t key[32]) {
    uint8_t computed[16];
    cfx_poly1271_avx2(computed, msg, len, key);
    uint8_t diff = 0;
    for (int i = 0; i < 16; i++) diff |= computed[i] ^ tag[i];
    CFX_MEMZERO_S(computed, sizeof(computed));
    return diff ? -1 : 0;
}

#endif /* CFX_HAVE_AVX2 */
