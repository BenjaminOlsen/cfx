/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ed25519.h - Ed25519 digital signatures (RFC 8032)
 *
 * Ed25519 is a deterministic Schnorr signature scheme using the twisted
 * Edwards curve -x² + y² = 1 + d·x²·y² over GF(2^255-19).
 *
 * Key sizes:
 *   - Seed: 32 bytes (random)
 *   - Secret key: 64 bytes (expanded from seed)
 *   - Public key: 32 bytes
 *   - Signature: 64 bytes
 *
 * All operations are constant-time with respect to secret data.
 */

#ifndef CFX_ED25519_H
#define CFX_ED25519_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ED25519_SEED_SIZE       32
#define ED25519_PUBLIC_KEY_SIZE 32
#define ED25519_SECRET_KEY_SIZE 64
#define ED25519_SIGNATURE_SIZE  64

void cfx_ed25519_create_keypair(uint8_t pk[32], uint8_t sk[64], const uint8_t seed[32]);  /* derive keypair from seed */
void cfx_ed25519_sign(uint8_t sig[64], const uint8_t* msg, size_t msg_len, const uint8_t sk[64]);
int cfx_ed25519_verify(const uint8_t sig[64], const uint8_t* msg, size_t msg_len, const uint8_t pk[32]);  /* 0 = valid */
void cfx_ed25519_get_public_key(uint8_t pk[32], const uint8_t sk[64]);  /* extract pk from sk */

#ifdef __cplusplus
}
#endif

#endif /* CFX_ED25519_H */
