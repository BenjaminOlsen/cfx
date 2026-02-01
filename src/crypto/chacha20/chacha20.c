/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/chacha20.h"
#include "cfx/memory.h"
#include "chacha20_backend.h"

#include <string.h>
#if defined(CFX_CAP_AVX2)
#include <immintrin.h>
#endif

#define _EXPA 0x61707865u
#define _ND_3 0x3320646eu
#define _2_BY 0x79622d32u
#define _TE_K 0x6b206574u

void cfx_chacha20_ctx_init(cfx_chacha20_ctx_t *ctx, const uint8_t key[32], const uint8_t nonce[12]) {
    cfx_chacha20_state_t *st = (cfx_chacha20_state_t *)ctx->opaque;
#if defined(CFX_CAP_AVX2)
    /* SoA layout: broadcast each word to all 8 lanes (32-byte aligned) */
    __m256i *row = (__m256i *)st->s;
    row[0]  = _mm256_set1_epi32((int)_EXPA);
    row[1]  = _mm256_set1_epi32((int)_ND_3);
    row[2]  = _mm256_set1_epi32((int)_2_BY);
    row[3]  = _mm256_set1_epi32((int)_TE_K);
    row[4]  = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 0));
    row[5]  = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 4));
    row[6]  = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 8));
    row[7]  = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 12));
    row[8]  = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 16));
    row[9]  = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 20));
    row[10] = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 24));
    row[11] = _mm256_set1_epi32((int)CFX_LOAD32_LE(key + 28));
    row[12] = _mm256_setzero_si256();
    row[13] = _mm256_set1_epi32((int)CFX_LOAD32_LE(nonce + 0));
    row[14] = _mm256_set1_epi32((int)CFX_LOAD32_LE(nonce + 4));
    row[15] = _mm256_set1_epi32((int)CFX_LOAD32_LE(nonce + 8));
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

/* increment 96-bit nonce (s[13..15]) as little-endian integer */
void cfx_chacha20_ctx_inc_nonce(cfx_chacha20_ctx_t *ctx) {
    cfx_chacha20_state_t *st = (cfx_chacha20_state_t *)ctx->opaque;
#if defined(CFX_CAP_AVX2)
    for (int j = 0; j < 8; ++j) {
        if (++st->s[13][j] != 0) continue;
        if (++st->s[14][j] != 0) continue;
        ++st->s[15][j];
    }
#else
    if (++st->s[13] != 0) return;
    if (++st->s[14] != 0) return;
    ++st->s[15];
#endif
}

void cfx_chacha20_block(cfx_chacha20_ctx_t *ctx, uint32_t counter, uint8_t out[64]) {
    cfx_chacha20_block_impl((cfx_chacha20_state_t *)ctx->opaque, counter, out);
}

void cfx_chacha20_block4(cfx_chacha20_ctx_t *ctx, uint32_t counter, uint8_t out[4][64]) {
    cfx_chacha20_block4_impl((cfx_chacha20_state_t *)ctx->opaque, counter, out);
}

void cfx_chacha20_block8(cfx_chacha20_ctx_t *ctx, uint32_t counter, uint8_t out[8][64]) {
    cfx_chacha20_block8_impl((cfx_chacha20_state_t *)ctx->opaque, counter, out);
}

void cfx_chacha20_encrypt_ctx(cfx_chacha20_ctx_t *ctx, uint32_t *counter, const uint8_t *pt, size_t pt_len, uint8_t *ct) {
    /* 8 blocks (512 bytes) at a time */
    while (pt_len >= 512) {
        uint8_t ks[8][64];
        cfx_chacha20_block8(ctx, *counter, ks);
        *counter += 8;
#if defined(CFX_CAP_AVX2)
        for (int j = 0; j < 16; ++j) {
            __m256i k = _mm256_loadu_si256((const __m256i *)((const uint8_t *)ks + j * 32));
            __m256i p = _mm256_loadu_si256((const __m256i *)(pt + j * 32));
            _mm256_storeu_si256((__m256i *)(ct + j * 32), _mm256_xor_si256(p, k));
        }
#else
        for (int blk = 0; blk < 8; ++blk) {
            for (size_t i = 0; i < 64; ++i) {
                ct[blk * 64 + i] = pt[blk * 64 + i] ^ ks[blk][i];
            }
        }
#endif
        
        pt += 512;
        ct += 512;
        pt_len -= 512;
    }
    
    /* 4 blocks - 256 bytes at a time */
    while (pt_len >= 256) {
        uint8_t ks[4][64];
        cfx_chacha20_block4(ctx, *counter, ks);
        *counter += 4;
#if defined(CFX_CAP_AVX2)
        for (int j = 0; j < 8; ++j) {
            __m256i k = _mm256_loadu_si256((const __m256i *)((const uint8_t *)ks + j * 32));
            __m256i p = _mm256_loadu_si256((const __m256i *)(pt + j * 32));
            _mm256_storeu_si256((__m256i *)(ct + j * 32), _mm256_xor_si256(p, k));
        }
#else
        for (int blk = 0; blk < 4; ++blk) {
            for (size_t i = 0; i < 64; ++i) {
                ct[blk * 64 + i] = pt[blk * 64 + i] ^ ks[blk][i];
            }
        }
#endif
        pt += 256;
        ct += 256;
        pt_len -= 256;
    }

    /* tail: one block at a time */
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
    const uint8_t *pt, size_t pt_len, uint8_t *ct) {
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

void cfx_xchacha20_ctx_init(cfx_chacha20_ctx_t *ctx, const uint8_t key[32], const uint8_t nonce[24]) {
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};
    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);
    cfx_chacha20_ctx_init(ctx, subkey, subnonce);
    CFX_MEMZERO_S(subkey, sizeof(subkey));
}

void cfx_xchacha20_encrypt_ctx(cfx_chacha20_ctx_t *ctx, uint32_t *counter, const uint8_t *pt, size_t pt_len, uint8_t *ct) {
    cfx_chacha20_encrypt_ctx(ctx, counter, pt, pt_len, ct);
}

void cfx_xchacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[24],
    const uint8_t *pt, size_t pt_len, uint8_t *ct) {
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};
    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);
    cfx_chacha20_encrypt(subkey, counter, subnonce, pt, pt_len, ct);
    CFX_MEMZERO_S(subkey, sizeof(subkey));
}
