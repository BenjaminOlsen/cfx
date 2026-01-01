/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * sha512.h - SHA-512 hash function
 *
 * Implementation of SHA-512 as specified in FIPS 180-4.
 * Required for Ed25519 signatures.
 */

#ifndef CFX_SHA512_H
#define CFX_SHA512_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_SHA512_DIGEST_SIZE 64
#define CFX_SHA512_BLOCK_SIZE  128
#define CFX_SHA512_CTX_SIZE    224u

typedef union {
    uint8_t  opaque[CFX_SHA512_CTX_SIZE];
    uint64_t aligner;
} cfx_sha512_ctx_t;

void cfx_sha512_init(cfx_sha512_ctx_t* ctx);
void cfx_sha512_update(cfx_sha512_ctx_t* ctx, const uint8_t* data, size_t len);
void cfx_sha512_final(cfx_sha512_ctx_t* ctx, uint8_t out[64]);
void cfx_sha512(uint8_t out[64], const uint8_t* data, size_t len);  /* one-shot */

#ifdef __cplusplus
}
#endif

#endif /* CFX_SHA512_H */
