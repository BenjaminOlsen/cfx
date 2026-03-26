/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_HMAC_H
#define CFX_HMAC_H

#include "cfx/arch.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* HMAC-SHA256 */

#define CFX_HMAC_SHA256_CTX_SIZE 192u

typedef union {
    CFX_ALIGNAS(CFX_CTX_ALIGN) uint8_t opaque[CFX_HMAC_SHA256_CTX_SIZE];
} cfx_hmac_sha256_ctx;

void cfx_hmac_sha256_init(cfx_hmac_sha256_ctx *ctx, const uint8_t *key, size_t key_len);
void cfx_hmac_sha256_update(cfx_hmac_sha256_ctx *ctx, const uint8_t *data, size_t len);
void cfx_hmac_sha256_final(cfx_hmac_sha256_ctx *ctx, uint8_t out[32]);

void cfx_hmac_sha256(uint8_t out[32], const uint8_t *key, size_t key_len, const uint8_t *data, size_t data_len);

/* HMAC-SHA512 */

#define CFX_HMAC_SHA512_CTX_SIZE 352u

typedef union {
    CFX_ALIGNAS(CFX_CTX_ALIGN) uint8_t opaque[CFX_HMAC_SHA512_CTX_SIZE];
} cfx_hmac_sha512_ctx;

void cfx_hmac_sha512_init(cfx_hmac_sha512_ctx *ctx, const uint8_t *key, size_t key_len);
void cfx_hmac_sha512_update(cfx_hmac_sha512_ctx *ctx, const uint8_t *data, size_t len);
void cfx_hmac_sha512_final(cfx_hmac_sha512_ctx *ctx, uint8_t out[64]);

void cfx_hmac_sha512(uint8_t out[64],
                      const uint8_t *key, size_t key_len,
                      const uint8_t *data, size_t data_len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_HMAC_H */
