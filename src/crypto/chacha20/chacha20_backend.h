/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ChaCha20 backend interface
 *
 * This header declares the target-specific functions that are implemented
 * in backend directories (base/, x86_avx2/, etc.). The facade in chacha20.c
 * calls these *_impl functions.
 */

#ifndef CFX_CHACHA20_BACKEND_H
#define CFX_CHACHA20_BACKEND_H

#include <stdint.h>

/*
 * Generate 8 ChaCha20 blocks in parallel.
 *
 * This is the primary function that benefits from SIMD optimization.
 * The base implementation falls back to calling cfx_chacha20_block_rfc8439
 * eight times, while AVX2 can process all 8 blocks simultaneously.
 *
 * Parameters:
 *   key     - 32-byte key
 *   counter - starting block counter (blocks use counter, counter+1, ..., counter+7)
 *   nonce   - 12-byte nonce
 *   out     - output buffer for 8 blocks (512 bytes total)
 */
void cfx_chacha20_block8_impl(const uint8_t key[32],
                              uint32_t counter,
                              const uint8_t nonce[12],
                              uint8_t out[8][64]);

#endif /* CFX_CHACHA20_BACKEND_H */
