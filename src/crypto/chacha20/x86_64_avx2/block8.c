/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "../chacha20_backend.h"
#include "cfx/memory.h"

#include <immintrin.h>

static inline __m256i rotl32(__m256i x, int n) {
    return _mm256_or_si256(_mm256_slli_epi32(x, n), _mm256_srli_epi32(x, 32 - n));
}

#define QR(a, b, c, d) do { \
    a = _mm256_add_epi32(a, b); d = _mm256_xor_si256(d, a); d = rotl32(d, 16); \
    c = _mm256_add_epi32(c, d); b = _mm256_xor_si256(b, c); b = rotl32(b, 12); \
    a = _mm256_add_epi32(a, b); d = _mm256_xor_si256(d, a); d = rotl32(d,  8); \
    c = _mm256_add_epi32(c, d); b = _mm256_xor_si256(b, c); b = rotl32(b,  7); \
} while (0)

#define LANE32(v, idx) ((uint32_t)_mm_extract_epi32(v, idx))

static inline void transpose(const __m256i x[16], uint32_t out[8][16]) {
    __m128i lo[16], hi[16];
    for (int w = 0; w < 16; ++w) {
        lo[w] = _mm256_extracti128_si256(x[w], 0);
        hi[w] = _mm256_extracti128_si256(x[w], 1);
    }
    for (int w = 0; w < 16; ++w) {
        out[0][w] = LANE32(lo[w], 0);
        out[1][w] = LANE32(lo[w], 1);
        out[2][w] = LANE32(lo[w], 2);
        out[3][w] = LANE32(lo[w], 3);
        out[4][w] = LANE32(hi[w], 0);
        out[5][w] = LANE32(hi[w], 1);
        out[6][w] = LANE32(hi[w], 2);
        out[7][w] = LANE32(hi[w], 3);
    }
}

void cfx_chacha20_block8_impl(const cfx_chacha20_state_t* ctx, uint32_t counter, uint8_t out[8][64]) {
    const __m256i lane = _mm256_setr_epi32(0, 1, 2, 3, 4, 5, 6, 7);
    __m256i x[16], s[16];

    /* Load pre-broadcast state from SoA layout - no per-call broadcast overhead */
    for (int i = 0; i < 16; ++i)
        s[i] = _mm256_loadu_si256((const __m256i*)ctx->s[i]);

    /* Set counters: counter+0, counter+1, ..., counter+7 */
    s[12] = _mm256_add_epi32(_mm256_set1_epi32((int32_t)counter), lane);

    for (int i = 0; i < 16; ++i) x[i] = s[i];

    for (int i = 0; i < 10; ++i) {
        QR(x[0], x[4], x[ 8], x[12]);
        QR(x[1], x[5], x[ 9], x[13]);
        QR(x[2], x[6], x[10], x[14]);
        QR(x[3], x[7], x[11], x[15]);
        QR(x[0], x[5], x[10], x[15]);
        QR(x[1], x[6], x[11], x[12]);
        QR(x[2], x[7], x[ 8], x[13]);
        QR(x[3], x[4], x[ 9], x[14]);
    }

    for (int i = 0; i < 16; ++i)
        x[i] = _mm256_add_epi32(x[i], s[i]);

    transpose(x, (uint32_t (*)[16])out);

    CFX_MEMZERO_S(x, sizeof(x));
    CFX_MEMZERO_S(s, sizeof(s));
}

#undef QR
#undef LANE32
