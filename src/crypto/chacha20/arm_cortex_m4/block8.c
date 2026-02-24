/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "../chacha20_backend.h"

void cfx_chacha20_block8_impl(const cfx_chacha20_state_t *ctx, uint32_t counter, uint8_t out[8][64]) {
    for (int i = 0; i < 8; ++i) {
        cfx_chacha20_block_impl(ctx, counter + (uint32_t)i, out[i]);
    }
}
