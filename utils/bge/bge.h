/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_BGE_H
#define CFX_BGE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Encrypt plaintext into the binary BGE stream format.
 *
 * On success, returns 0 and transfers ownership of *output to the caller.
 * Release it with cfx_bge_free().
 * Empty plaintext and embedded '\0' bytes are ok.
 * The passphrase is treated as an arbitrary byte sequence.
 */
int cfx_bge_encrypt(const uint8_t *plaintext, size_t plaintext_len,
                    const uint8_t *passphrase, size_t passphrase_len,
                    uint8_t **output, size_t *output_len);

/*
 * Decrypt binary or armored BGE data.
 *
 * Returns 0 on success and transfers ownership of *plaintext to the caller.
 * Free the buffer with cfx_bge_free.
 * Returns -1 for invalid arguments/allocation failures, -2 for malformed
 * input, and -3 when authentication fails.
 */
int cfx_bge_decrypt(const uint8_t *input, size_t input_len,
                    const uint8_t *passphrase, size_t passphrase_len,
                    uint8_t **plaintext, size_t *plaintext_len);

/* Clear and release a buffer returned by this API. */
void cfx_bge_free(void *buffer, size_t buffer_len);

#ifdef __cplusplus
}
#endif

#endif
