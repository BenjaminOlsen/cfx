/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/chacha20.h"
#include "cfx/memory.h"
#include "chacha20_backend.h"

#include <string.h>

#define _EXPA 0x61707865u
#define _ND_3 0x3320646eu
#define _2_BY 0x79622d32u
#define _TE_K 0x6b206574u

void cfx_chacha20_ctx_init(cfx_chacha20_ctx_t* ctx, const uint8_t key[32], const uint8_t nonce[12]) {
    cfx_chacha20_state_t* st = (cfx_chacha20_state_t*)ctx->opaque;
#if defined(CFX_CAP_AVX2)
    /* SoA layout: broadcast each word to all 8 lanes */
    for (int j = 0; j < 8; ++j) {
        st->s[0][j]  = _EXPA;
        st->s[1][j]  = _ND_3;
        st->s[2][j]  = _2_BY;
        st->s[3][j]  = _TE_K;
        st->s[4][j]  = CFX_LOAD32_LE(key + 0);
        st->s[5][j]  = CFX_LOAD32_LE(key + 4);
        st->s[6][j]  = CFX_LOAD32_LE(key + 8);
        st->s[7][j]  = CFX_LOAD32_LE(key + 12);
        st->s[8][j]  = CFX_LOAD32_LE(key + 16);
        st->s[9][j]  = CFX_LOAD32_LE(key + 20);
        st->s[10][j] = CFX_LOAD32_LE(key + 24);
        st->s[11][j] = CFX_LOAD32_LE(key + 28);
        st->s[12][j] = 0;
        st->s[13][j] = CFX_LOAD32_LE(nonce + 0);
        st->s[14][j] = CFX_LOAD32_LE(nonce + 4);
        st->s[15][j] = CFX_LOAD32_LE(nonce + 8);
    }
#else
    st->s[0]  = _EXPA;
    st->s[1]  = _ND_3;
    st->s[2]  = _2_BY;
    st->s[3]  = _TE_K;
    st->s[4]  = CFX_LOAD32_LE(key + 0);
    st->s[5]  = CFX_LOAD32_LE(key + 4);
    st->s[6]  = CFX_LOAD32_LE(key + 8);
    st->s[7]  = CFX_LOAD32_LE(key + 12);
    st->s[8]  = CFX_LOAD32_LE(key + 16);
    st->s[9]  = CFX_LOAD32_LE(key + 20);
    st->s[10] = CFX_LOAD32_LE(key + 24);
    st->s[11] = CFX_LOAD32_LE(key + 28);
    st->s[12] = 0;
    st->s[13] = CFX_LOAD32_LE(nonce + 0);
    st->s[14] = CFX_LOAD32_LE(nonce + 4);
    st->s[15] = CFX_LOAD32_LE(nonce + 8);
#endif
}

void cfx_chacha20_block(cfx_chacha20_ctx_t* ctx, uint32_t counter, uint8_t out[64]) {
    cfx_chacha20_block_impl((cfx_chacha20_state_t*)ctx->opaque, counter, out);
}

void cfx_chacha20_block4(cfx_chacha20_ctx_t* ctx, uint32_t counter, uint8_t out[4][64]) {
    cfx_chacha20_block4_impl((cfx_chacha20_state_t*)ctx->opaque, counter, out);
}

void cfx_chacha20_block8(cfx_chacha20_ctx_t* ctx, uint32_t counter, uint8_t out[8][64]) {
    cfx_chacha20_block8_impl((cfx_chacha20_state_t*)ctx->opaque, counter, out);
}

void cfx_chacha20_encrypt_ctx(cfx_chacha20_ctx_t* ctx, uint32_t* counter, const uint8_t* pt, size_t pt_len, uint8_t* ct) {
    uint8_t ks[64];
    while (pt_len) {
        cfx_chacha20_block(ctx, *counter, ks);
        ++*counter;
        size_t take = pt_len < 64 ? pt_len : 64;
        for (size_t i = 0; i < take; ++i) ct[i] = pt[i] ^ ks[i];
        pt += take;
        ct += take;
        pt_len -= take;
    }
    CFX_MEMZERO_S(ks, sizeof(ks));
}

void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t* pt, size_t pt_len, uint8_t* ct) {
    cfx_chacha20_ctx_t ctx;
    cfx_chacha20_ctx_init(&ctx, key, nonce);
    cfx_chacha20_encrypt_ctx(&ctx, &counter, pt, pt_len, ct);
    CFX_MEMZERO_S(&ctx, sizeof(ctx));
}

#define ROTL32(x, n) ((uint32_t)(((x) << (n)) | ((x) >> (32 - (n)))))
#define QR(a, b, c, d) \
    a += b; d ^= a; d = ROTL32(d, 16); \
    c += d; b ^= c; b = ROTL32(b, 12); \
    a += b; d ^= a; d = ROTL32(d,  8); \
    c += d; b ^= c; b = ROTL32(b,  7);

void cfx_hchacha20(uint8_t out[32], const uint8_t key[32], const uint8_t nonce[16]) {
    uint32_t w[16];
    w[0]  = _EXPA;
    w[1]  = _ND_3;
    w[2]  = _2_BY;
    w[3]  = _TE_K;
    w[4]  = CFX_LOAD32_LE(key + 0);
    w[5]  = CFX_LOAD32_LE(key + 4);
    w[6]  = CFX_LOAD32_LE(key + 8);
    w[7]  = CFX_LOAD32_LE(key + 12);
    w[8]  = CFX_LOAD32_LE(key + 16);
    w[9]  = CFX_LOAD32_LE(key + 20);
    w[10] = CFX_LOAD32_LE(key + 24);
    w[11] = CFX_LOAD32_LE(key + 28);
    w[12] = CFX_LOAD32_LE(nonce + 0);
    w[13] = CFX_LOAD32_LE(nonce + 4);
    w[14] = CFX_LOAD32_LE(nonce + 8);
    w[15] = CFX_LOAD32_LE(nonce + 12);

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

    CFX_STORE32_LE(out + 0,  w[0]);
    CFX_STORE32_LE(out + 4,  w[1]);
    CFX_STORE32_LE(out + 8,  w[2]);
    CFX_STORE32_LE(out + 12, w[3]);
    CFX_STORE32_LE(out + 16, w[12]);
    CFX_STORE32_LE(out + 20, w[13]);
    CFX_STORE32_LE(out + 24, w[14]);
    CFX_STORE32_LE(out + 28, w[15]);
    CFX_MEMZERO_S(w, sizeof(w));
}

#undef QR
#undef ROTL32

void cfx_xchacha20_ctx_init(cfx_chacha20_ctx_t* ctx, const uint8_t key[32], const uint8_t nonce[24]) {
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};
    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);
    cfx_chacha20_ctx_init(ctx, subkey, subnonce);
    CFX_MEMZERO_S(subkey, sizeof(subkey));
}

void cfx_xchacha20_encrypt_ctx(cfx_chacha20_ctx_t* ctx, uint32_t* counter, const uint8_t* pt, size_t pt_len, uint8_t* ct) {
    cfx_chacha20_encrypt_ctx(ctx, counter, pt, pt_len, ct);
}

void cfx_xchacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[24],
                           const uint8_t* pt, size_t pt_len, uint8_t* ct) {
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};
    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);
    cfx_chacha20_encrypt(subkey, counter, subnonce, pt, pt_len, ct);
    CFX_MEMZERO_S(subkey, sizeof(subkey));
}
