/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_CHACHA20_H
#define CFX_CHACHA20_H

#include "cfx/arch.h"
#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Context size depends on target:
 *   - AVX2: 512 bytes (pre-broadcast SoA layout for 8-lane SIMD, 32-byte aligned)
 *   - Others: 64 bytes (scalar layout)
 */
#if defined(CFX_CAP_AVX2)
#define CFX_CHACHA20_CTX_SIZE 512
#define CFX_CHACHA20_CTX_ALIGN 32
#else
#define CFX_CHACHA20_CTX_SIZE 64
#define CFX_CHACHA20_CTX_ALIGN 8
#endif

typedef union {
    CFX_ALIGNAS(CFX_CHACHA20_CTX_ALIGN) uint8_t opaque[CFX_CHACHA20_CTX_SIZE];
} cfx_chacha20_ctx_t;

void cfx_chacha20_ctx_init(cfx_chacha20_ctx_t *ctx, const uint8_t key[32], const uint8_t nonce[12]);
void cfx_chacha20_ctx_inc_nonce(cfx_chacha20_ctx_t *ctx);
void cfx_chacha20_block(cfx_chacha20_ctx_t *ctx, uint32_t counter, uint8_t out[64]);
void cfx_chacha20_block4(cfx_chacha20_ctx_t *ctx, uint32_t counter, uint8_t out[4][64]);
void cfx_chacha20_block8(cfx_chacha20_ctx_t *ctx, uint32_t counter, uint8_t out[8][64]);
void cfx_chacha20_encrypt_ctx(cfx_chacha20_ctx_t *ctx, uint32_t *counter, const uint8_t *pt, size_t pt_len, uint8_t *ct);

void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
    const uint8_t *pt, size_t pt_len, uint8_t *ct);

void cfx_hchacha20(uint8_t out[32], const uint8_t key[32], const uint8_t nonce[16]);

void cfx_xchacha20_ctx_init(cfx_chacha20_ctx_t *ctx, const uint8_t key[32], const uint8_t nonce[24]);
void cfx_xchacha20_encrypt_ctx(cfx_chacha20_ctx_t *ctx, uint32_t *counter, const uint8_t *pt, size_t pt_len, uint8_t *ct);
void cfx_xchacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[24],
    const uint8_t *pt, size_t pt_len, uint8_t *ct);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_CHACHA20_H */
