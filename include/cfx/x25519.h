/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * x25519.h - X25519 Diffie-Hellman key exchange (RFC 7748)
 *
 * Montgomery ladder on Curve25519, x-coordinate only.
 * All operations are constant-time.
 */

#ifndef CFX_X25519_H
#define CFX_X25519_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* key sizes */
#define X25519_KEY_SIZE     32
#define X25519_SCALAR_SIZE  32

extern const uint8_t cfx_x25519_basepoint[32];  /* basepoint x = 9 */

int cfx_x25519(uint8_t shared[32], const uint8_t scalar[32], const uint8_t point[32]);  /* full DH, returns 0 on success */
void cfx_x25519_base(uint8_t public_key[32], const uint8_t private_key[32]);  /* public = scalar * basepoint */
void cfx_x25519_scalarmult(uint8_t out[32], const uint8_t scalar[32], const uint8_t point[32]);  /* raw, no clamping */

#ifdef __cplusplus
}
#endif

#endif /* CFX_X25519_H */
