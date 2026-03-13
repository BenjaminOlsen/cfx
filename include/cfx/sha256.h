/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_SHA256_H
#define CFX_SHA256_H

#include "cfx/arch.h"
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_SHA256_CTX_SIZE 128u

typedef union {
    CFX_ALIGNAS(CFX_CTX_ALIGN) uint8_t opaque[CFX_SHA256_CTX_SIZE];
} cfx_sha256_ctx;

void cfx_sha256_init(cfx_sha256_ctx *ctx);
void cfx_sha256_update(cfx_sha256_ctx *ctx, const uint8_t *data, size_t len);
void cfx_sha256_final(cfx_sha256_ctx *ctx, uint8_t out[32]);
void cfx_sha256(uint8_t out[32], const uint8_t *data, size_t len);

#ifdef __cplusplus
}
#endif

#endif
