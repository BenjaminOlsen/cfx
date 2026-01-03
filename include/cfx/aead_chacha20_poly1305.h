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
    uint8_t*        ct,             /* out: cyphertext, length = pt_len */
    uint8_t         tag[16],        /* out */
    const uint8_t*  pt,             /* in: plaintext,  length = pt_len */
    size_t          pt_len,
    const uint8_t*  aad,            /* in, may be NULL (iff aad_len == 0) */
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[12]
);


int cfx_chacha20_poly1305_decrypt(
    uint8_t*        pt,             /* plaintext  - out, length = ct_len */
    const uint8_t*  ct,             /* cyphertext - in,  length = ct_len */
    size_t          ct_len,
    const uint8_t*  aad,            /* in, may be NULL iff aad_len == 0 */
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[12],
    const uint8_t   tag[16]
);

/* XChaCha20-Poly1305: uses 24-byte nonce (safe for random nonces) */
int cfx_xchacha20_poly1305_encrypt(
    uint8_t*        ct,             /* out: ciphertext, length = pt_len */
    uint8_t         tag[16],        /* out */
    const uint8_t*  pt,             /* in: plaintext,  length = pt_len */
    size_t          pt_len,
    const uint8_t*  aad,            /* in, may be NULL iff aad_len == 0 */
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[24]
);

int cfx_xchacha20_poly1305_decrypt(
    uint8_t*        pt,             /* plaintext  - out, length = ct_len */
    const uint8_t*  ct,             /* ciphertext - in,  length = ct_len */
    size_t          ct_len,
    const uint8_t*  aad,            /* in, may be NULL iff aad_len == 0 */
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[24],
    const uint8_t   tag[16]
);

#ifdef __cplusplus
}
#endif

#endif /* CFX_AEAD_CHACHA20_POLY1305_H */
