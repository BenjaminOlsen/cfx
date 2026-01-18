/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Portable ChaCha20 block8 implementation.
 *
 * Falls back to generating 8 blocks sequentially using the scalar
 * cfx_chacha20_block_rfc8439 function.
 *
 * Self-guarding: compiles when no optimized backend is selected.
 */

#if !defined(CFX_TARGET_X86_64_AVX2) && \
    !defined(CFX_TARGET_X86_64_AVX512) && \
    !defined(CFX_TARGET_ARM_CORTEX_M4) && \
    !defined(CFX_TARGET_AARCH64_NEON)

#include "chacha20_backend.h"
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

#endif /* !CFX_TARGET_* */
