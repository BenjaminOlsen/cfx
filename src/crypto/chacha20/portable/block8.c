/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Portable ChaCha20 block8 implementation
 *
 * Falls back to generating 8 blocks sequentially using the scalar
 * cfx_chacha20_block_rfc8439 function.
 */

#include "../chacha20_backend.h"
#include "cfx/chacha20.h"

void cfx_chacha20_block8_impl(const uint8_t key[32],
                              uint32_t counter,
                              const uint8_t nonce[12],
                              uint8_t out[8][64])
{
    for (int i = 0; i < 8; ++i) {
        cfx_chacha20_block_rfc8439(key, counter + (uint32_t)i, nonce, out[i]);
    }
}
