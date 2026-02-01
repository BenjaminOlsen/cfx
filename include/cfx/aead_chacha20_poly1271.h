/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_AEAD_CHACHA20_POLY1271_H
#define CFX_AEAD_CHACHA20_POLY1271_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* chacha20-poly1271 AEAD (like RFC8439 but with poly1271, pads to 15-byte boundaries) */
int cfx_chacha20_poly1271_encrypt(
    uint8_t *ct,                    /* out: ciphertext, length = pt_len */
    uint8_t tag[16],                /* out */
    const uint8_t *pt,              /* in: plaintext, length = pt_len */
    size_t pt_len,
    const uint8_t *aad,             /* in, may be NULL iff aad_len == 0 */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[12]
    );

int cfx_chacha20_poly1271_decrypt(
    uint8_t *pt,                    /* out: plaintext, length = ct_len */
    const uint8_t *ct,              /* in: ciphertext, length = ct_len */
    size_t ct_len,
    const uint8_t *aad,             /* in, may be NULL iff aad_len == 0 */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[12],
    const uint8_t tag[16]
    );

/* xchacha20-poly1271: 24-byte nonce (safe for random nonces) */
int cfx_xchacha20_poly1271_encrypt(
    uint8_t *ct,
    uint8_t tag[16],
    const uint8_t *pt,
    size_t pt_len,
    const uint8_t *aad,
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[24]
    );

int cfx_xchacha20_poly1271_decrypt(
    uint8_t *pt,
    const uint8_t *ct,
    size_t ct_len,
    const uint8_t *aad,
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[24],
    const uint8_t tag[16]
    );

/*
 * Streaming AEAD (STREAM construction) over XChaCha20-Poly1271.
 *
 * Same STREAM design as the Poly1305 variant: per-chunk nonces derived from
 * base_nonce XOR'd with chunk counter + final flag, 64KB chunks, 16-byte tags.
 */
#define CFX_STREAM_1271_CHUNK_SIZE  65536u
#define CFX_STREAM_1271_TAG_SIZE    16u

int cfx_stream_xchacha20_poly1271_encrypt_chunk(
    uint8_t *ct, uint8_t tag[16],
    const uint8_t *pt, size_t pt_len,
    uint64_t chunk_counter, int is_final,
    const uint8_t key[32], const uint8_t base_nonce[24]);

int cfx_stream_xchacha20_poly1271_decrypt_chunk(
    uint8_t *pt,
    const uint8_t *ct, size_t ct_len,
    const uint8_t tag[16],
    uint64_t chunk_counter, int is_final,
    const uint8_t key[32], const uint8_t base_nonce[24]);

#ifdef __cplusplus
}
#endif

#endif /* CFX_AEAD_CHACHA20_POLY1271_H */
