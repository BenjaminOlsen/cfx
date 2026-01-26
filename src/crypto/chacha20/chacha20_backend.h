/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_CHACHA20_BACKEND_H
#define CFX_CHACHA20_BACKEND_H

#include "cfx/arch.h"
#include <stdint.h>

/*
 * Internal state layout depends on target:
 *   - AVX2: s[16][8] - SoA layout, pre-broadcast for 8-lane SIMD (32-byte aligned)
 *   - Others: s[16] - scalar layout
 */
#if defined(CFX_CAP_AVX2)
typedef struct {
    CFX_ALIGNAS(32) uint32_t s[16][8];
} cfx_chacha20_state_t;
#else
typedef struct {
    uint32_t s[16];
} cfx_chacha20_state_t;
#endif

void cfx_chacha20_block_impl(const cfx_chacha20_state_t *ctx, uint32_t counter, uint8_t out[64]);
void cfx_chacha20_block4_impl(const cfx_chacha20_state_t *ctx, uint32_t counter, uint8_t out[4][64]);
void cfx_chacha20_block8_impl(const cfx_chacha20_state_t *ctx, uint32_t counter, uint8_t out[8][64]);

#endif /* CFX_CHACHA20_BACKEND_H */
