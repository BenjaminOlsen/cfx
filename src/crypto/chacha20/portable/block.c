/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "../chacha20_backend.h"
#include "cfx/memory.h"

#define ROTL32(x, n) ((uint32_t)(((x) << (n)) | ((x) >> (32 - (n)))))

#define QR(a, b, c, d) \
        a += b; d ^= a; d = ROTL32(d, 16); \
        c += d; b ^= c; b = ROTL32(b, 12); \
        a += b; d ^= a; d = ROTL32(d,  8); \
        c += d; b ^= c; b = ROTL32(b,  7);

void cfx_chacha20_block_impl(const cfx_chacha20_state_t *ctx, uint32_t counter, uint8_t out[64]) {
    uint32_t w[16];
    for (int i = 0; i < 16; ++i) w[i] = ctx->s[i];
    w[12] = counter;

    for (int i = 0; i < 10; ++i) {
        QR(w[0], w[4], w[ 8], w[12])
        QR(w[1], w[5], w[ 9], w[13])
        QR(w[2], w[6], w[10], w[14])
        QR(w[3], w[7], w[11], w[15])
        QR(w[0], w[5], w[10], w[15])
        QR(w[1], w[6], w[11], w[12])
        QR(w[2], w[7], w[ 8], w[13])
        QR(w[3], w[4], w[ 9], w[14])
    }

    for (int i = 0; i < 12; ++i) w[i] += ctx->s[i];
    w[12] += counter;
    for (int i = 13; i < 16; ++i) w[i] += ctx->s[i];

    for (int i = 0; i < 16; ++i) CFX_STORE32_LE(out + 4*i, w[i]);
}

#undef QR
#undef ROTL32
