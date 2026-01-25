/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/*
 * Poly1271 AVX2 Implementation
 *
 * This implementation uses AVX2 SIMD instructions to process 4 message blocks
 * in parallel.
 *
 * KEY DIFFERENCES FROM SCALAR:
 * - Uses radix-2^26 representation (5 limbs of 26 bits each)
 * - Processes 4 independent block×r^k products simultaneously
 * - Uses vectorized gather for efficient block loading
 *
 * Based on the reference implementation from github.com/poly1271/poly1271
 */

#include "cfx/poly1271.h"
#include "cfx/memory.h"

#if CFX_HAVE_AVX2

#include <immintrin.h>
#include <string.h>

/*============================================================================
 * CONSTANTS
 *============================================================================*/

/*
 * Radix-2^26 representation:
 * A 127-bit value is stored as 5 limbs of ~26 bits each:
 *   value = limb0 + limb1*2^26 + limb2*2^52 + limb3*2^78 + limb4*2^104
 *
 * Why radix-2^26?
 * - Fits in 32-bit lanes for SIMD multiply
 * - 26*5 = 130 bits > 127 bits needed
 * - Allows lazy carry propagation without overflow
 */
#define MASK26 ((1ULL << 26) - 1)  /* 26 bits: 0x03FFFFFF */
#define MASK23 ((1ULL << 23) - 1)  /* 23 bits for top limb (127 - 4*26 = 23) */

/*============================================================================
 * RADIX CONVERSION
 *============================================================================*/

/*
 * Decode a 15-byte message block into 5 radix-2^26 limbs.
 *
 * The 0x01 delimiter is placed at bit 120 (byte 14, then shifted).
 * This ensures non-zero blocks and encodes length.
 */
static inline void decode_block_to_limbs(uint64_t limbs[5], const uint8_t *block) {
    /* Load as two 64-bit values (little-endian) */
    uint64_t lo = (uint64_t)block[0] | ((uint64_t)block[1] << 8) |
                  ((uint64_t)block[2] << 16) | ((uint64_t)block[3] << 24) |
                  ((uint64_t)block[4] << 32) | ((uint64_t)block[5] << 40) |
                  ((uint64_t)block[6] << 48) | ((uint64_t)block[7] << 56);
    uint64_t hi = (uint64_t)block[8] | ((uint64_t)block[9] << 8) |
                  ((uint64_t)block[10] << 16) | ((uint64_t)block[11] << 24) |
                  ((uint64_t)block[12] << 32) | ((uint64_t)block[13] << 40) |
                  ((uint64_t)block[14] << 48) | (1ULL << 56); /* delimiter at bit 120 */

    /* Split 128 bits into 5 x 26-bit limbs (with 23 bits in top limb) */
    limbs[0] = lo & MASK26;                          /* bits 0-25 */
    limbs[1] = (lo >> 26) & MASK26;                  /* bits 26-51 */
    limbs[2] = ((lo >> 52) | (hi << 12)) & MASK26;   /* bits 52-77 (crosses boundary) */
    limbs[3] = (hi >> 14) & MASK26;                  /* bits 78-103 */
    limbs[4] = (hi >> 40) & MASK26;                  /* bits 104-127 (only 23 bits used) */
}

/*
 * Decode 4 consecutive 15-byte blocks directly into SIMD vectors using gather.
 *
 * Each output vector contains one limb from each of the 4 blocks:
 *   out[0] = [block0.limb0, block1.limb0, block2.limb0, block3.limb0]
 */
static inline void decode_4blocks_to_vec(__m256i out[5], const uint8_t *msg) {
    /*
     * Block layout: blocks are at offsets 0, 15, 30, 45 bytes.
     * We gather the low 8 bytes and high 7 bytes of each block.
     */
    const __m256i lo_idx = _mm256_set_epi64x(45, 30, 15, 0);   /* offsets for lo[0..7] */
    const __m256i hi_idx = _mm256_set_epi64x(53, 38, 23, 8);   /* offsets for hi[8..14] */

    /* Gather 8 bytes from each of 4 positions */
    __m256i lo = _mm256_i64gather_epi64((const long long *)msg, lo_idx, 1);
    __m256i hi = _mm256_i64gather_epi64((const long long *)msg, hi_idx, 1);

    /* Mask hi to 56 bits (7 bytes) and add delimiter at bit 56 */
    const __m256i mask56 = _mm256_set1_epi64x(0x00FFFFFFFFFFFFFFULL);
    const __m256i delim = _mm256_set1_epi64x(1ULL << 56);
    hi = _mm256_or_si256(_mm256_and_si256(hi, mask56), delim);

    /* Convert to radix-2^26 (same logic as scalar, but vectorized) */
    const __m256i mask26 = _mm256_set1_epi64x(MASK26);

    out[0] = _mm256_and_si256(lo, mask26);                                          /* bits 0-25 */
    out[1] = _mm256_and_si256(_mm256_srli_epi64(lo, 26), mask26);                    /* bits 26-51 */
    __m256i lo_hi = _mm256_srli_epi64(lo, 52);
    __m256i hi_lo = _mm256_slli_epi64(hi, 12);
    out[2] = _mm256_and_si256(_mm256_or_si256(lo_hi, hi_lo), mask26);                /* bits 52-77 */
    out[3] = _mm256_and_si256(_mm256_srli_epi64(hi, 14), mask26);                    /* bits 78-103 */
    out[4] = _mm256_and_si256(_mm256_srli_epi64(hi, 40), mask26);                    /* bits 104-127 */
}

/*
 * Decode a partial block (< 15 bytes) into limbs.
 */
static inline void decode_partial_to_limbs(uint64_t limbs[5], const uint8_t *block, size_t len) {
    uint8_t pad[16] = {0};
    memcpy(pad, block, len);
    pad[len] = 0x01;  /* delimiter immediately after message */

    uint64_t lo = (uint64_t)pad[0] | ((uint64_t)pad[1] << 8) |
                  ((uint64_t)pad[2] << 16) | ((uint64_t)pad[3] << 24) |
                  ((uint64_t)pad[4] << 32) | ((uint64_t)pad[5] << 40) |
                  ((uint64_t)pad[6] << 48) | ((uint64_t)pad[7] << 56);
    uint64_t hi = (uint64_t)pad[8] | ((uint64_t)pad[9] << 8) |
                  ((uint64_t)pad[10] << 16) | ((uint64_t)pad[11] << 24) |
                  ((uint64_t)pad[12] << 32) | ((uint64_t)pad[13] << 40) |
                  ((uint64_t)pad[14] << 48) | ((uint64_t)pad[15] << 56);

    limbs[0] = lo & MASK26;
    limbs[1] = (lo >> 26) & MASK26;
    limbs[2] = ((lo >> 52) | (hi << 12)) & MASK26;
    limbs[3] = (hi >> 14) & MASK26;
    limbs[4] = (hi >> 40) & MASK26;
}

/*
 * Convert between radix-2^26 (5 limbs) and radix-2^64 (2 limbs) representations.
 */
static inline void convert_to_limbs(uint64_t limbs[5], const uint64_t v[2]) {
    uint64_t lo = v[0], hi = v[1];
    limbs[0] = lo & MASK26;
    limbs[1] = (lo >> 26) & MASK26;
    limbs[2] = ((lo >> 52) | (hi << 12)) & MASK26;
    limbs[3] = (hi >> 14) & MASK26;
    limbs[4] = (hi >> 40) & MASK26;
}

static inline void convert_from_limbs(uint64_t v[2], const uint64_t limbs[5]) {
    v[0] = limbs[0] | (limbs[1] << 26) | (limbs[2] << 52);
    v[1] = (limbs[2] >> 12) | (limbs[3] << 14) | (limbs[4] << 40);
}

/*============================================================================
 * SCALAR RADIX-2^26 ARITHMETIC
 *============================================================================
 * Used for single-block processing and key setup.
 */

/*
 * Propagate carries between limbs.
 */
static inline void carry26_scalar(uint64_t limbs[5]) {
    uint64_t c;
    c = limbs[0] >> 26; limbs[0] &= MASK26; limbs[1] += c;
    c = limbs[1] >> 26; limbs[1] &= MASK26; limbs[2] += c;
    c = limbs[2] >> 26; limbs[2] &= MASK26; limbs[3] += c;
    c = limbs[3] >> 26; limbs[3] &= MASK26; limbs[4] += c;
}

/*
 * Reduce mod 2^127 - 1.
 *
 * In radix-2^26, bit 127 is at position 127 - 4*26 = 23 within limb4.
 * Since 2^127 = 1 (mod p), we add these bits back to limb0.
 */
static inline void reduce_p_scalar(uint64_t limbs[5]) {
    uint64_t c = limbs[4] >> 23;  /* bits >= 127 */
    limbs[4] &= MASK23;           /* keep only bits 0-22 */
    limbs[0] += c;                /* add back (since 2^127 = 1) */
    carry26_scalar(limbs);

    /* May need one more reduction if carry propagated to limb4 */
    c = limbs[4] >> 23;
    limbs[4] &= MASK23;
    limbs[0] += c;
    carry26_scalar(limbs);
}

/*
 * Schoolbook multiply in radix-2^26, mod p.
 */
static void mul_limbs_scalar(uint64_t out[5], const uint64_t a[5], const uint64_t b[5]) {
    /* Schoolbook: 25 multiplications */
    uint64_t d0 = a[0] * b[0];
    uint64_t d1 = a[0] * b[1] + a[1] * b[0];
    uint64_t d2 = a[0] * b[2] + a[1] * b[1] + a[2] * b[0];
    uint64_t d3 = a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0];
    uint64_t d4 = a[0] * b[4] + a[1] * b[3] + a[2] * b[2] + a[3] * b[1] + a[4] * b[0];
    uint64_t d5 = a[1] * b[4] + a[2] * b[3] + a[3] * b[2] + a[4] * b[1];
    uint64_t d6 = a[2] * b[4] + a[3] * b[3] + a[4] * b[2];
    uint64_t d7 = a[3] * b[4] + a[4] * b[3];
    uint64_t d8 = a[4] * b[4];

    /* Fold high limbs: 2^130 = 8 (mod p), so d5*2^130 = d5*8, etc. */
    d0 += d5 << 3;
    d1 += d6 << 3;
    d2 += d7 << 3;
    d3 += d8 << 3;

    /* Carry propagation */
    uint64_t c;
    c = d0 >> 26; d0 &= MASK26; d1 += c;
    c = d1 >> 26; d1 &= MASK26; d2 += c;
    c = d2 >> 26; d2 &= MASK26; d3 += c;
    c = d3 >> 26; d3 &= MASK26; d4 += c;

    out[0] = d0; out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4;
    reduce_p_scalar(out);
}

/*
 * Square in radix-2^26, mod p.
 */
static void square_limbs_scalar(uint64_t out[5], const uint64_t a[5]) {
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
    c = d0 >> 26; d0 &= MASK26; d1 += c;
    c = d1 >> 26; d1 &= MASK26; d2 += c;
    c = d2 >> 26; d2 &= MASK26; d3 += c;
    c = d3 >> 26; d3 &= MASK26; d4 += c;

    out[0] = d0; out[1] = d1; out[2] = d2; out[3] = d3; out[4] = d4;
    reduce_p_scalar(out);
}

/*============================================================================
 * UTILITY FUNCTIONS
 *============================================================================*/

static inline uint64_t load64_le(const uint8_t *p) {
    return (uint64_t)p[0] | ((uint64_t)p[1] << 8) | ((uint64_t)p[2] << 16) |
           ((uint64_t)p[3] << 24) | ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}

static inline void store64_le(uint8_t *p, uint64_t v) {
    for (int i = 0; i < 8; i++) p[i] = (uint8_t)(v >> (i * 8));
}

/*
 * Key clamping: same as scalar implementation.
 */
static void clamp_r(uint8_t r[16]) {
    r[3] &= 0x0F; r[7] &= 0x0F; r[11] &= 0x0F; r[15] &= 0x07;
    r[4] &= 0xFC; r[8] &= 0xFC; r[12] &= 0xFC;
}

/*============================================================================
 * AVX2 VECTORIZED ARITHMETIC
 *============================================================================
 *
 * Each __m256i register holds 4 x 64-bit values, one from each block.
 * We compute 4 independent products simultaneously.
 */

/*
 * Multiply 4 parallel accumulators by r^4 (broadcast to all lanes).
 */
static inline void mul_r4_avx2(__m256i acc[5], const __m256i r4v[5]) {
    const __m256i r0 = r4v[0], r1 = r4v[1], r2 = r4v[2], r3 = r4v[3], r4 = r4v[4];
    __m256i a0 = acc[0], a1 = acc[1], a2 = acc[2], a3 = acc[3], a4 = acc[4];

    /* Schoolbook multiply (vectorized): 25 SIMD multiplications */
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

    /* Fold high limbs: 2^130 = 8 (mod p) */
    d0 = _mm256_add_epi64(d0, _mm256_slli_epi64(d5, 3));
    d1 = _mm256_add_epi64(d1, _mm256_slli_epi64(d6, 3));
    d2 = _mm256_add_epi64(d2, _mm256_slli_epi64(d7, 3));
    d3 = _mm256_add_epi64(d3, _mm256_slli_epi64(d8, 3));

    /* Carry propagation (vectorized) */
    __m256i mask26 = _mm256_set1_epi64x(MASK26);
    __m256i mask23 = _mm256_set1_epi64x(MASK23);
    __m256i c;

    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);
    c = _mm256_srli_epi64(d1, 26); d1 = _mm256_and_si256(d1, mask26); d2 = _mm256_add_epi64(d2, c);
    c = _mm256_srli_epi64(d2, 26); d2 = _mm256_and_si256(d2, mask26); d3 = _mm256_add_epi64(d3, c);
    c = _mm256_srli_epi64(d3, 26); d3 = _mm256_and_si256(d3, mask26); d4 = _mm256_add_epi64(d4, c);
    c = _mm256_srli_epi64(d4, 23); d4 = _mm256_and_si256(d4, mask23); d0 = _mm256_add_epi64(d0, c);
    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);

    acc[0] = d0; acc[1] = d1; acc[2] = d2; acc[3] = d3; acc[4] = d4;
}

/*
 * Add 4 blocks (already in vector format) to the accumulator.
 */
static inline void add_4blocks_vec_avx2(__m256i acc[5], const __m256i blocks[5]) {
    acc[0] = _mm256_add_epi64(acc[0], blocks[0]);
    acc[1] = _mm256_add_epi64(acc[1], blocks[1]);
    acc[2] = _mm256_add_epi64(acc[2], blocks[2]);
    acc[3] = _mm256_add_epi64(acc[3], blocks[3]);
    acc[4] = _mm256_add_epi64(acc[4], blocks[4]);
}

/*
 * Horizontal sum: add all 4 lanes of a __m256i to produce a single uint64_t.
 */
static inline uint64_t hsum_epi64(__m256i v) {
    /* v = [A, B, C, D] */
    __m256i hi = _mm256_permute4x64_epi64(v, 0x4E);  /* [C, D, A, B] */
    __m256i sum1 = _mm256_add_epi64(v, hi);          /* [A+C, B+D, ...] */
    __m128i lo128 = _mm256_castsi256_si128(sum1);    /* [A+C, B+D] */
    __m128i hi128 = _mm_unpackhi_epi64(lo128, lo128);/* [B+D, B+D] */
    __m128i sum2 = _mm_add_epi64(lo128, hi128);      /* [A+B+C+D, ...] */
    return (uint64_t)_mm_cvtsi128_si64(sum2);
}

/*
 * Combine 4 lanes into a single accumulator value.
 *
 * result = lane0*r^4 + lane1*r^3 + lane2*r^2 + lane3*r
 *
 * rv contains [r, r^2, r^3, r^4] per limb, and we multiply then sum horizontally.
 */
static void combine_4lanes_avx2(uint64_t out[5], const __m256i acc[5], const __m256i rv[5]) {
    /* Schoolbook multiply with different r powers per lane */
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

    /* Fold high limbs */
    d0 = _mm256_add_epi64(d0, _mm256_slli_epi64(d5, 3));
    d1 = _mm256_add_epi64(d1, _mm256_slli_epi64(d6, 3));
    d2 = _mm256_add_epi64(d2, _mm256_slli_epi64(d7, 3));
    d3 = _mm256_add_epi64(d3, _mm256_slli_epi64(d8, 3));

    /* Carry propagation */
    __m256i mask26 = _mm256_set1_epi64x(MASK26);
    __m256i mask23 = _mm256_set1_epi64x(MASK23);
    __m256i c;

    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);
    c = _mm256_srli_epi64(d1, 26); d1 = _mm256_and_si256(d1, mask26); d2 = _mm256_add_epi64(d2, c);
    c = _mm256_srli_epi64(d2, 26); d2 = _mm256_and_si256(d2, mask26); d3 = _mm256_add_epi64(d3, c);
    c = _mm256_srli_epi64(d3, 26); d3 = _mm256_and_si256(d3, mask26); d4 = _mm256_add_epi64(d4, c);
    c = _mm256_srli_epi64(d4, 23); d4 = _mm256_and_si256(d4, mask23); d0 = _mm256_add_epi64(d0, c);
    c = _mm256_srli_epi64(d0, 26); d0 = _mm256_and_si256(d0, mask26); d1 = _mm256_add_epi64(d1, c);

    /* Horizontal sum: combine all 4 lanes */
    out[0] = hsum_epi64(d0);
    out[1] = hsum_epi64(d1);
    out[2] = hsum_epi64(d2);
    out[3] = hsum_epi64(d3);
    out[4] = hsum_epi64(d4);

    carry26_scalar(out);
    reduce_p_scalar(out);
}

/*============================================================================
 * BLOCK PROCESSING
 *============================================================================*/

/*
 * Process n rounds of 4 blocks each using SIMD.
 */
static void process_nblocks_avx2(uint64_t acc[5], const uint8_t *msg, int nrounds,
    const __m256i r4v[5], const __m256i rv[5]) {
    __m256i pacc[5];
    for (int i = 0; i < 5; i++)
        pacc[i] = _mm256_set_epi64x(0, 0, 0, (int64_t)acc[i]);

    for (int round = 0; round < nrounds; round++) {
        __m256i blocks[5];
        decode_4blocks_to_vec(blocks, msg + round * 60);  /* 60 = 4 * 15 bytes */
        add_4blocks_vec_avx2(pacc, blocks);

        if (round < nrounds - 1)
            mul_r4_avx2(pacc, r4v);
    }

    combine_4lanes_avx2(acc, pacc, rv);
}

/* Convenience wrappers for different batch sizes */
static inline void process_4blocks_avx2(uint64_t acc[5], const uint8_t *msg,
    const __m256i r4v[5], const __m256i rv[5]) {
    process_nblocks_avx2(acc, msg, 1, r4v, rv);
}

static inline void process_16blocks_avx2(uint64_t acc[5], const uint8_t *msg,
    const __m256i r4v[5], const __m256i rv[5]) {
    process_nblocks_avx2(acc, msg, 4, r4v, rv);
}

static inline void process_64blocks_avx2(uint64_t acc[5], const uint8_t *msg,
    const __m256i r4v[5], const __m256i rv[5]) {
    process_nblocks_avx2(acc, msg, 16, r4v, rv);
}

/*
 * Process a single block using scalar arithmetic.
 */
static void process_1block_scalar(uint64_t acc[5], const uint8_t *m, const uint64_t r[5]) {
    uint64_t b[5];
    decode_block_to_limbs(b, m);
    for (int i = 0; i < 5; i++) acc[i] += b[i];
    carry26_scalar(acc);
    mul_limbs_scalar(acc, acc, r);
}

/*============================================================================
 * FINALIZATION
 *============================================================================*/

/*
 * Finalize: reduce, check for p, add s, output tag.
 */
static void finalize_avx2(uint64_t acc[5], const uint64_t s[2], uint8_t tag[16]) {
    reduce_p_scalar(acc);

    uint64_t v[2];
    convert_from_limbs(v, acc);

    /* Constant-time check if acc == p (2^127 - 1) */
    uint64_t eq_lo = ~(v[0] ^ 0xFFFFFFFFFFFFFFFFULL);
    eq_lo &= eq_lo >> 32; eq_lo &= eq_lo >> 16; eq_lo &= eq_lo >> 8;
    eq_lo &= eq_lo >> 4;  eq_lo &= eq_lo >> 2;  eq_lo &= eq_lo >> 1;

    uint64_t eq_hi = ~(v[1] ^ 0x7FFFFFFFFFFFFFFFULL);
    eq_hi &= eq_hi >> 32; eq_hi &= eq_hi >> 16; eq_hi &= eq_hi >> 8;
    eq_hi &= eq_hi >> 4;  eq_hi &= eq_hi >> 2;  eq_hi &= eq_hi >> 1;

    uint64_t wrap = (uint64_t)0 - (eq_lo & eq_hi & 1);
    v[0] += wrap & 1;
    uint64_t c = (v[0] == 0) & (wrap & 1);
    v[1] = (v[1] + c) & 0x7FFFFFFFFFFFFFFFULL;

    /* Add s mod 2^128 */
    uint64_t t0 = v[0] + s[0];
    c = (t0 < v[0]);
    uint64_t t1 = v[1] + s[1] + c;

    store64_le(tag, t0);
    store64_le(tag + 8, t1);
}

/*============================================================================
 * PUBLIC API
 *============================================================================*/

void cfx_poly1271_avx2_init(cfx_poly1271_avx2_ctx_t *ctx, const uint8_t key[32]) {
    uint8_t rc[16];
    memcpy(rc, key, 16);
    clamp_r(rc);

    uint64_t r64[2] = {load64_le(rc), load64_le(rc + 8)};
    convert_to_limbs(ctx->r, r64);

    /* Precompute powers of r */
    square_limbs_scalar(ctx->r2, ctx->r);
    mul_limbs_scalar(ctx->r3, ctx->r2, ctx->r);
    mul_limbs_scalar(ctx->r4, ctx->r2, ctx->r2);

    /*
     * Set up vector constants:
     * - rv[i] contains [r^4, r^3, r^2, r] for limb i (for combine_4lanes)
     * - r4v[i] contains [r^4, r^4, r^4, r^4] for limb i (for mul_r4)
     */
    for (int i = 0; i < 5; i++) {
        ctx->rv[i] = _mm256_set_epi64x((int64_t)ctx->r[i], (int64_t)ctx->r2[i],
            (int64_t)ctx->r3[i], (int64_t)ctx->r4[i]);
        ctx->r4v[i] = _mm256_set1_epi64x((int64_t)ctx->r4[i]);
    }

    ctx->s[0] = load64_le(key + 16);
    ctx->s[1] = load64_le(key + 24);

    for (int i = 0; i < 5; i++) ctx->acc[i] = 0;
    ctx->buflen = 0;

    CFX_MEMZERO_S(rc, 16);
}

void cfx_poly1271_avx2_update(cfx_poly1271_avx2_ctx_t *ctx, const uint8_t *msg, size_t len) {
    /* Handle buffered partial block */
    if (ctx->buflen) {
        size_t need = 15 - ctx->buflen;
        if (len < need) {
            memcpy(ctx->buf + ctx->buflen, msg, len);
            ctx->buflen += (uint8_t)len;
            return;
        }
        memcpy(ctx->buf + ctx->buflen, msg, need);
        process_1block_scalar(ctx->acc, ctx->buf, ctx->r);
        msg += need;
        len -= need;
        ctx->buflen = 0;
    }

    /* Process in large batches for efficiency */
    while (len >= 960) {   /* 64 blocks */
        process_64blocks_avx2(ctx->acc, msg, ctx->r4v, ctx->rv);
        msg += 960; len -= 960;
    }
    while (len >= 240) {   /* 16 blocks */
        process_16blocks_avx2(ctx->acc, msg, ctx->r4v, ctx->rv);
        msg += 240; len -= 240;
    }
    while (len >= 60) {    /* 4 blocks */
        process_4blocks_avx2(ctx->acc, msg, ctx->r4v, ctx->rv);
        msg += 60; len -= 60;
    }

    /* Process remaining blocks one at a time */
    while (len >= 15) {
        process_1block_scalar(ctx->acc, msg, ctx->r);
        msg += 15;
        len -= 15;
    }

    /* Buffer any remaining partial block */
    if (len) {
        memcpy(ctx->buf, msg, len);
        ctx->buflen = (uint8_t)len;
    }
}

void cfx_poly1271_avx2_finish(cfx_poly1271_avx2_ctx_t *ctx, uint8_t tag[16]) {
    if (ctx->buflen) {
        uint64_t b[5];
        decode_partial_to_limbs(b, ctx->buf, ctx->buflen);
        for (int i = 0; i < 5; i++) ctx->acc[i] += b[i];
        carry26_scalar(ctx->acc);
        mul_limbs_scalar(ctx->acc, ctx->acc, ctx->r);
    }
    finalize_avx2(ctx->acc, ctx->s, tag);
    CFX_MEMZERO_S(ctx, sizeof(*ctx));
}

void cfx_poly1271_avx2(uint8_t tag[16], const uint8_t *msg, size_t len, const uint8_t key[32]) {
    cfx_poly1271_avx2_ctx_t ctx;
    cfx_poly1271_avx2_init(&ctx, key);
    cfx_poly1271_avx2_update(&ctx, msg, len);
    cfx_poly1271_avx2_finish(&ctx, tag);
}

int cfx_poly1271_avx2_verify(const uint8_t tag[16], const uint8_t *msg, size_t len, const uint8_t key[32]) {
    uint8_t computed[16];
    cfx_poly1271_avx2(computed, msg, len, key);
    uint8_t diff = 0;
    for (int i = 0; i < 16; i++) diff |= computed[i] ^ tag[i];
    CFX_MEMZERO_S(computed, 16);
    return diff ? -1 : 0;
}

#endif /* CFX_HAVE_AVX2 */
