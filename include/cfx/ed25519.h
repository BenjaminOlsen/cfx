/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ed25519.h - Ed25519 digital signatures (RFC 8032)
 *
 * Ed25519 is a deterministic Schnorr signature scheme using the twisted
 * Edwards curve -x² + y² = 1 + d·x²·y² over GF(2^255-19).
 *
 * Key sizes:
 *   - Seed: 32 bytes (random)
 *   - Secret key: 64 bytes = seed || public_key (bundled to prevent misuse)
 *   - Public key: 32 bytes
 *   - Signature: 64 bytes
 *
 * The secret key bundles the seed with its public key to prevent a class of
 * attacks where signing the same message with mismatched keys leaks the
 * private scalar. The sign function extracts pk from sk internally.
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
#define ED25519_PREHASH_SIZE    64  /* SHA-512 output for Ed25519ph */

/* Ed25519 - sign/verify full message (RFC 8032 pure mode) */
void cfx_ed25519_create_keypair(uint8_t pk[32], uint8_t sk[64], const uint8_t seed[32]);  /* derive keypair from seed */
void cfx_ed25519_sign(uint8_t sig[64], const uint8_t *msg, size_t msg_len, const uint8_t sk[64]);
int cfx_ed25519_verify(const uint8_t sig[64], const uint8_t *msg, size_t msg_len, const uint8_t pk[32]);  /* 0 = valid */
void cfx_ed25519_get_public_key(uint8_t pk[32], const uint8_t sk[64]);  /* extract pk from sk */

/*
 * Ed25519ph - pre-hashed mode (RFC 8032 section 5.1)
 *
 * For signing large files: first compute PH = SHA-512(message),
 * then sign/verify the pre-hash with domain separation.
 *
 * The domain prefix prevents confusion between Ed25519 and Ed25519ph signatures:
 *   dom2(1, "") = "SigEd25519 no Ed25519 collisions" || 0x01 || 0x00
 *
 * Usage:
 *   uint8_t prehash[64];
 *   cfx_sha512(prehash, large_file, file_len);  // or stream with sha512_ctx
 *   cfx_ed25519ph_sign(sig, prehash, sk);
 */
void cfx_ed25519ph_sign(uint8_t sig[64], const uint8_t prehash[64], const uint8_t sk[64]);
int cfx_ed25519ph_verify(const uint8_t sig[64], const uint8_t prehash[64], const uint8_t pk[32]);  /* 0 = valid */

#ifdef __cplusplus
}
#endif

#endif /* CFX_ED25519_H */
