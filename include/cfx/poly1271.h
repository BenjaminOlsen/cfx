/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * poly1271.h - Polynomial MAC using Mersenne prime 2^127 - 1
 *
 * Similar to Poly1305 but uses the Mersenne prime M127 = 2^127 - 1
 * for faster modular reduction (just addition, no multiplication).
 *
 * Block size: 15 bytes (vs 16 for Poly1305)
 * Security: ~123 bits (vs ~103 for Poly1305)
 * Key: 32 bytes (r: 16 bytes, s: 16 bytes)
 * Tag: 16 bytes
 */

#ifndef CFX_POLY1271_H
#define CFX_POLY1271_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_POLY1271_BLOCK_SIZE 15
#define CFX_POLY1271_KEY_SIZE   32
#define CFX_POLY1271_TAG_SIZE   16

typedef struct {
    uint64_t acc[3];   /* accumulator (up to ~190 bits during computation) */
    uint64_t r[2];     /* r key (127 bits, clamped) */
    uint64_t r2[2];    /* r^2 mod M127 (precomputed) */
    uint64_t r3[2];    /* r^3 mod M127 (precomputed) */
    uint64_t r4[2];    /* r^4 mod M127 (precomputed for 4-block processing) */
    uint64_t s[2];     /* s key (128 bits) */
    uint8_t  buf[15];  /* partial block buffer */
    uint8_t  buflen;   /* bytes in buffer */
    uint8_t  blkcnt;   /* blocks since last full reduction (for lazy reduction) */
} cfx_poly1271_ctx_t;

/* initialize with 32-byte key (first 16 = r, second 16 = s) */
void cfx_poly1271_init(cfx_poly1271_ctx_t* ctx, const uint8_t key[32]);

/* process message data (streaming) */
void cfx_poly1271_update(cfx_poly1271_ctx_t* ctx, const uint8_t* msg, size_t len);

/* finalize and output 16-byte tag */
void cfx_poly1271_finish(cfx_poly1271_ctx_t* ctx, uint8_t tag[16]);

/* one-shot MAC */
void cfx_poly1271(uint8_t tag[16], const uint8_t* msg, size_t len, const uint8_t key[32]);

/* AVX2 implementation - radix-2^26, 4-way parallel */
#include "cfx/arch.h"

#if CFX_HAVE_AVX2
#include <immintrin.h>

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable: 4324) /* structure padded due to alignment */
#endif

typedef struct {
    uint64_t acc[5];     /* accumulator, radix-26 */
    uint64_t r[5];       /* r, radix-26 */
    uint64_t r2[5];      /* r^2 */
    uint64_t r3[5];      /* r^3 */
    uint64_t r4[5];      /* r^4 */
    __m256i rv[5];       /* rv[i] = {r[i], r2[i], r3[i], r4[i]} for combine */
    __m256i r4v[5];      /* r4v[i] = broadcast(r4[i]) for mul_r4 */
    uint64_t s[2];       /* s key */
    uint8_t  buf[60];    /* partial block buffer (4 blocks) */
    uint8_t  buflen;
} cfx_poly1271_avx2_ctx_t;

void cfx_poly1271_avx2_init(cfx_poly1271_avx2_ctx_t* ctx, const uint8_t key[32]);
void cfx_poly1271_avx2_update(cfx_poly1271_avx2_ctx_t* ctx, const uint8_t* msg, size_t len);
void cfx_poly1271_avx2_finish(cfx_poly1271_avx2_ctx_t* ctx, uint8_t tag[16]);
void cfx_poly1271_avx2(uint8_t tag[16], const uint8_t* msg, size_t len, const uint8_t key[32]);
int cfx_poly1271_avx2_verify(const uint8_t tag[16], const uint8_t* msg, size_t len, const uint8_t key[32]);

#ifdef _MSC_VER
#pragma warning(pop)
#endif

#endif /* CFX_HAVE_AVX2 */

#ifdef __cplusplus
}
#endif

#endif /* CFX_POLY1271_H */
