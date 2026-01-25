/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * sha3.h - SHA-3 (Keccak) cryptographic hash (FIPS 202)
 *
 * Implements:
 *   SHA3-224, SHA3-256, SHA3-384, SHA3-512 (fixed output)
 *   SHAKE128, SHAKE256 (extendable output / XOF)
 *
 * Based on Keccak sponge construction with 1600-bit state.
 */

#ifndef CFX_SHA3_H
#define CFX_SHA3_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* state size: 1600 bits = 200 bytes = 25 x 64-bit lanes */
#define CFX_SHA3_STATE_SIZE  200u
#define CFX_SHA3_CTX_SIZE    400u

typedef union {
    uint8_t opaque[CFX_SHA3_CTX_SIZE];
    uint64_t aligner;
} cfx_sha3_ctx_t;

/* SHA3-224 (r=1152, c=448, output=28 bytes) */
void cfx_sha3_224_init(cfx_sha3_ctx_t *ctx);
void cfx_sha3_224_update(cfx_sha3_ctx_t *ctx, const void *data, size_t len);
void cfx_sha3_224_final(cfx_sha3_ctx_t *ctx, uint8_t out[28]);
void cfx_sha3_224(uint8_t out[28], const void *data, size_t len);

/* SHA3-256 (r=1088, c=512, output=32 bytes) */
void cfx_sha3_256_init(cfx_sha3_ctx_t *ctx);
void cfx_sha3_256_update(cfx_sha3_ctx_t *ctx, const void *data, size_t len);
void cfx_sha3_256_final(cfx_sha3_ctx_t *ctx, uint8_t out[32]);
void cfx_sha3_256(uint8_t out[32], const void *data, size_t len);

/* SHA3-384 (r=832, c=768, output=48 bytes) */
void cfx_sha3_384_init(cfx_sha3_ctx_t *ctx);
void cfx_sha3_384_update(cfx_sha3_ctx_t *ctx, const void *data, size_t len);
void cfx_sha3_384_final(cfx_sha3_ctx_t *ctx, uint8_t out[48]);
void cfx_sha3_384(uint8_t out[48], const void *data, size_t len);

/* SHA3-512 (r=576, c=1024, output=64 bytes) */
void cfx_sha3_512_init(cfx_sha3_ctx_t *ctx);
void cfx_sha3_512_update(cfx_sha3_ctx_t *ctx, const void *data, size_t len);
void cfx_sha3_512_final(cfx_sha3_ctx_t *ctx, uint8_t out[64]);
void cfx_sha3_512(uint8_t out[64], const void *data, size_t len);

/* SHAKE128 - extendable output function (r=1344, c=256) */
void cfx_shake128_init(cfx_sha3_ctx_t *ctx);
void cfx_shake128_update(cfx_sha3_ctx_t *ctx, const void *data, size_t len);
void cfx_shake128_final(cfx_sha3_ctx_t *ctx);  /* finalize absorbing */
void cfx_shake128_squeeze(cfx_sha3_ctx_t *ctx, void *out, size_t outlen);
void cfx_shake128(void *out, size_t outlen, const void *data, size_t len);

/* SHAKE256 - extendable output function (r=1088, c=512) */
void cfx_shake256_init(cfx_sha3_ctx_t *ctx);
void cfx_shake256_update(cfx_sha3_ctx_t *ctx, const void *data, size_t len);
void cfx_shake256_final(cfx_sha3_ctx_t *ctx);
void cfx_shake256_squeeze(cfx_sha3_ctx_t *ctx, void *out, size_t outlen);
void cfx_shake256(void *out, size_t outlen, const void *data, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_SHA3_H */
