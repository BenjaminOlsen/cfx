/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_AEAD_CHACHA20_POLY1305_H
#define CFX_AEAD_CHACHA20_POLY1305_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* implements the RFC8439 chacha20-poly1305 authenticated encryption with additional data algorithm */
int cfx_chacha20_poly1305_encrypt(
    uint8_t *ct,                    /* out: cyphertext, length = pt_len */
    uint8_t tag[16],                /* out */
    const uint8_t *pt,              /* in: plaintext,  length = pt_len */
    size_t pt_len,
    const uint8_t *aad,             /* in, may be NULL (iff aad_len == 0) */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[12]
    );


int cfx_chacha20_poly1305_decrypt(
    uint8_t *pt,                    /* plaintext  - out, length = ct_len */
    const uint8_t *ct,              /* cyphertext - in,  length = ct_len */
    size_t ct_len,
    const uint8_t *aad,             /* in, may be NULL iff aad_len == 0 */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[12],
    const uint8_t tag[16]
    );

/* XChaCha20-Poly1305: uses 24-byte nonce (safe for random nonces) */
int cfx_xchacha20_poly1305_encrypt(
    uint8_t *ct,                    /* out: ciphertext, length = pt_len */
    uint8_t tag[16],                /* out */
    const uint8_t *pt,              /* in: plaintext,  length = pt_len */
    size_t pt_len,
    const uint8_t *aad,             /* in, may be NULL iff aad_len == 0 */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[24]
    );

int cfx_xchacha20_poly1305_decrypt(
    uint8_t *pt,                    /* plaintext  - out, length = ct_len */
    const uint8_t *ct,              /* ciphertext - in,  length = ct_len */
    size_t ct_len,
    const uint8_t *aad,             /* in, may be NULL iff aad_len == 0 */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[24],
    const uint8_t tag[16]
    );

/*
 * Streaming AEAD (STREAM construction) over XChaCha20-Poly1305.
 *
 * encrypts/decrypts data in independent chunks, each authenticated with its
 * own Poly1305 tag. Per-chunk nonces are derived from a base nonce XOR'd with
 * the chunk counter and a final-chunk flag, preventing reordering and truncation.
 *
 * Max stream: 2^31 chunks * 64KB = 128 TB (top bit of byte 23 is the final flag).
 */
#define CFX_STREAM_CHUNK_SIZE  65536u
#define CFX_STREAM_TAG_SIZE    16u

/* Encrypt one chunk. ct must have room for pt_len bytes, tag for 16 bytes.
 * chunk_counter starts at 0, increments by 1 per chunk. is_final=1 for last chunk. */
int cfx_stream_xchacha20_poly1305_encrypt_chunk(
    uint8_t *ct, uint8_t tag[16],
    const uint8_t *pt, size_t pt_len,
    uint64_t chunk_counter, int is_final,
    const uint8_t key[32], const uint8_t base_nonce[24]);

/* Decrypt one chunk. Returns -1 on authentication failure.
 * chunk_counter and is_final must match what was used during encryption. */
int cfx_stream_xchacha20_poly1305_decrypt_chunk(
    uint8_t *pt,
    const uint8_t *ct, size_t ct_len,
    const uint8_t tag[16],
    uint64_t chunk_counter, int is_final,
    const uint8_t key[32], const uint8_t base_nonce[24]);

#ifdef __cplusplus
}
#endif

#endif /* CFX_AEAD_CHACHA20_POLY1305_H */
