/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * AVX2-optimized ChaCha20 block8 implementation.
 *
 * Generates 8 ChaCha20 blocks in parallel using 256-bit AVX2 vectors.
 * Each __m256i holds one word from all 8 blocks (SoA layout).
 *
 * Requires: -mavx2 compiler flag
 *
 * Self-guarding: only compiles when CFX_TARGET_X86_64_AVX2 or AVX512 is defined.
 */

#if defined(CFX_TARGET_X86_64_AVX2) || defined(CFX_TARGET_X86_64_AVX512)

#include "chacha20_backend.h"
#include "cfx/chacha20.h"

#include <immintrin.h>
#include <string.h>

#define _EXPA 0x61707865u
#define _ND_3 0x3320646eu
#define _2_BY 0x79622d32u
#define _TE_K 0x6b206574u

static inline __m256i cfx_mm256_rotl32(__m256i x, int n) {
    return _mm256_or_si256(_mm256_slli_epi32(x, n),
                           _mm256_srli_epi32(x, 32 - n));
}

#define MM256_QR(a,b,c,d) do {          \
    a = _mm256_add_epi32(a, b);         \
    d = _mm256_xor_si256(d, a);         \
    d = cfx_mm256_rotl32(d, 16);        \
    c = _mm256_add_epi32(c, d);         \
    b = _mm256_xor_si256(b, c);         \
    b = cfx_mm256_rotl32(b, 12);        \
    a = _mm256_add_epi32(a, b);         \
    d = _mm256_xor_si256(d, a);         \
    d = cfx_mm256_rotl32(d, 8);         \
    c = _mm256_add_epi32(c, d);         \
    b = _mm256_xor_si256(b, c);         \
    b = cfx_mm256_rotl32(b, 7);         \
} while (0)

/*
 * Transpose from SoA (x[w] = word w for all 8 blocks) to AoS (out[b] = block b).
 */
#define LANE32(v, idx) (uint32_t)_mm_extract_epi32(v, idx)

static inline void transpose_16x8_to_blocks(const __m256i x[16], uint32_t out[8][16]) {
    __m128i lo[16], hi[16];

    for (int w = 0; w < 16; ++w) {
        lo[w] = _mm256_extracti128_si256(x[w], 0); /* lanes 0..3 */
        hi[w] = _mm256_extracti128_si256(x[w], 1); /* lanes 4..7 */
    }

    uint32_t *dst;

#define EXTRACT_BLOCK(ARR, lane_idx, out_idx)                           \
    do {                                                                \
        dst = out[(out_idx)];                                           \
        for (int w = 0; w < 16; ++w) {                                  \
            dst[w] = LANE32((ARR)[w], (lane_idx));                      \
        }                                                               \
    } while (0)

    /* Blocks 0..3 from lo */
    EXTRACT_BLOCK(lo, 0, 0);
    EXTRACT_BLOCK(lo, 1, 1);
    EXTRACT_BLOCK(lo, 2, 2);
    EXTRACT_BLOCK(lo, 3, 3);

    /* Blocks 4..7 from hi */
    EXTRACT_BLOCK(hi, 0, 4);
    EXTRACT_BLOCK(hi, 1, 5);
    EXTRACT_BLOCK(hi, 2, 6);
    EXTRACT_BLOCK(hi, 3, 7);

#undef EXTRACT_BLOCK
}

#undef LANE32

void cfx_chacha20_block8_impl(const uint8_t key[32],
                              uint32_t counter,
                              const uint8_t nonce[12],
                              uint8_t out[8][64])
{
    /* Expand key/nonce to u32 */
    uint32_t k[8];
    uint32_t n[3];

    for (size_t i = 0; i < 8; ++i) {
        k[i] = CFX_LOAD32_LE(key + 4*i);
    }
    for (size_t i = 0; i < 3; ++i) {
        n[i] = CFX_LOAD32_LE(nonce + 4*i);
    }

    /* lane offsets {0,1,2,3,4,5,6,7} */
    const __m256i lane = _mm256_setr_epi32(0,1,2,3,4,5,6,7);

    __m256i x[16], o[16];

    x[0]  = _mm256_set1_epi32((int32_t)_EXPA);
    x[1]  = _mm256_set1_epi32((int32_t)_ND_3);
    x[2]  = _mm256_set1_epi32((int32_t)_2_BY);
    x[3]  = _mm256_set1_epi32((int32_t)_TE_K);
    x[4]  = _mm256_set1_epi32((int32_t)k[0]);
    x[5]  = _mm256_set1_epi32((int32_t)k[1]);
    x[6]  = _mm256_set1_epi32((int32_t)k[2]);
    x[7]  = _mm256_set1_epi32((int32_t)k[3]);
    x[8]  = _mm256_set1_epi32((int32_t)k[4]);
    x[9]  = _mm256_set1_epi32((int32_t)k[5]);
    x[10] = _mm256_set1_epi32((int32_t)k[6]);
    x[11] = _mm256_set1_epi32((int32_t)k[7]);

    /* lane 0..7 has counter value: {counter, counter+1, ... counter+7} */
    x[12] = _mm256_add_epi32(_mm256_set1_epi32((int32_t)counter), lane);

    x[13] = _mm256_set1_epi32((int32_t)n[0]);
    x[14] = _mm256_set1_epi32((int32_t)n[1]);
    x[15] = _mm256_set1_epi32((int32_t)n[2]);

    /* save originals for final addition */
    for (size_t i = 0; i < 16; ++i) {
        o[i] = x[i];
    }

    /* 20 rounds (10 double-rounds) */
    for (size_t i = 0; i < 10; ++i) {
        MM256_QR(x[0], x[4], x[8],  x[12]);
        MM256_QR(x[1], x[5], x[9],  x[13]);
        MM256_QR(x[2], x[6], x[10], x[14]);
        MM256_QR(x[3], x[7], x[11], x[15]);
        MM256_QR(x[0], x[5], x[10], x[15]);
        MM256_QR(x[1], x[6], x[11], x[12]);
        MM256_QR(x[2], x[7], x[8],  x[13]);
        MM256_QR(x[3], x[4], x[9],  x[14]);
    }

    /* add original state */
    for (size_t i = 0; i < 16; ++i) {
        x[i] = _mm256_add_epi32(x[i], o[i]);
    }

    /* transpose and store */
    uint32_t (*out_words)[16] = (uint32_t (*)[16])out;
    transpose_16x8_to_blocks(x, out_words);

    CFX_MEMZERO_S(k, sizeof(k));
    CFX_MEMZERO_S(x, sizeof(x));
    CFX_MEMZERO_S(o, sizeof(o));
}

#undef MM256_QR

#endif /* CFX_TARGET_X86_64_AVX2 || CFX_TARGET_X86_64_AVX512 */
