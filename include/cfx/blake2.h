/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * blake2.h - BLAKE2 cryptographic hash (RFC 7693)
 *
 * BLAKE2b: 64-bit optimized, up to 64-byte output, 128-byte blocks
 * BLAKE2s: 32-bit optimized, up to 32-byte output, 64-byte blocks
 *
 * Features:
 *   - Faster than MD5, SHA-1, SHA-2, SHA-3
 *   - No length extension attacks
 *   - Built-in keyed hashing (MAC mode)
 *   - Variable output length
 */

#ifndef CFX_BLAKE2_H
#define CFX_BLAKE2_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_BLAKE2B_BLOCKBYTES  128u
#define CFX_BLAKE2B_OUTBYTES    64u
#define CFX_BLAKE2B_KEYBYTES    64u
#define CFX_BLAKE2B_CTX_SIZE    256u

#define CFX_BLAKE2S_BLOCKBYTES  64u
#define CFX_BLAKE2S_OUTBYTES    32u
#define CFX_BLAKE2S_KEYBYTES    32u
#define CFX_BLAKE2S_CTX_SIZE    128u

typedef union {
    uint8_t opaque[CFX_BLAKE2B_CTX_SIZE];
    size_t aligner;
} cfx_blake2b_ctx_t;

typedef union {
    uint8_t opaque[CFX_BLAKE2S_CTX_SIZE];
    size_t aligner;
} cfx_blake2s_ctx_t;

/* BLAKE2b - 64-bit optimized */
int cfx_blake2b_init(cfx_blake2b_ctx_t *ctx, size_t outlen);
int cfx_blake2b_init_key(cfx_blake2b_ctx_t *ctx, size_t outlen,
    const void *key, size_t keylen);
int cfx_blake2b_update(cfx_blake2b_ctx_t *ctx, const void *in, size_t inlen);
int cfx_blake2b_final(cfx_blake2b_ctx_t *ctx, void *out);

/* one-shot */
int cfx_blake2b(void *out, size_t outlen,
    const void *in, size_t inlen,
    const void *key, size_t keylen);

/* BLAKE2s - 32-bit optimized */
int cfx_blake2s_init(cfx_blake2s_ctx_t *ctx, size_t outlen);
int cfx_blake2s_init_key(cfx_blake2s_ctx_t *ctx, size_t outlen,
    const void *key, size_t keylen);
int cfx_blake2s_update(cfx_blake2s_ctx_t *ctx, const void *in, size_t inlen);
int cfx_blake2s_final(cfx_blake2s_ctx_t *ctx, void *out);

/* one-shot */
int cfx_blake2s(void *out, size_t outlen,
    const void *in, size_t inlen,
    const void *key, size_t keylen);

#ifdef __cplusplus
}
#endif

#endif /* CFX_BLAKE2_H */
