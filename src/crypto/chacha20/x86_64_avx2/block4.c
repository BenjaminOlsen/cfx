/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "../chacha20_backend.h"

#include <immintrin.h>

#define ROTL16(x) _mm_or_si128(_mm_slli_epi32(x, 16), _mm_srli_epi32(x, 16))
#define ROTL12(x) _mm_or_si128(_mm_slli_epi32(x, 12), _mm_srli_epi32(x, 20))
#define ROTL8(x)  _mm_or_si128(_mm_slli_epi32(x,  8), _mm_srli_epi32(x, 24))
#define ROTL7(x)  _mm_or_si128(_mm_slli_epi32(x,  7), _mm_srli_epi32(x, 25))

#define QR(a, b, c, d) do { \
            a = _mm_add_epi32(a, b); d = _mm_xor_si128(d, a); d = ROTL16(d); \
            c = _mm_add_epi32(c, d); b = _mm_xor_si128(b, c); b = ROTL12(b); \
            a = _mm_add_epi32(a, b); d = _mm_xor_si128(d, a); d = ROTL8(d);  \
            c = _mm_add_epi32(c, d); b = _mm_xor_si128(b, c); b = ROTL7(b);  \
} while (0)

void cfx_chacha20_block4_impl(const cfx_chacha20_state_t *ctx, uint32_t counter, uint8_t out[4][64]) {
    __m128i x[16], s[16];

    /* Load first 4 lanes from SoA layout */
    for (int i = 0; i < 16; ++i)
        s[i] = _mm_loadu_si128((const __m128i *)ctx->s[i]);

    /* Set counters: counter+0, counter+1, counter+2, counter+3 */
    s[12] = _mm_add_epi32(_mm_set1_epi32((int32_t)counter), _mm_set_epi32(3, 2, 1, 0));

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
        x[i] = _mm_add_epi32(x[i], s[i]);

    /* Transpose and store: convert from SoA to 4 separate 64-byte blocks */
    for (int blk = 0; blk < 4; ++blk) {
        for (int i = 0; i < 16; ++i) {
            uint32_t w = (uint32_t)_mm_extract_epi32(x[i], 0);
            x[i] = _mm_srli_si128(x[i], 4);
            out[blk][4*i + 0] = (uint8_t)(w);
            out[blk][4*i + 1] = (uint8_t)(w >> 8);
            out[blk][4*i + 2] = (uint8_t)(w >> 16);
            out[blk][4*i + 3] = (uint8_t)(w >> 24);
        }
    }
}
