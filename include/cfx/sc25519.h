/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * sc25519.h - Scalar arithmetic modulo the group order L
 *
 * L = 2^252 + 27742317777372353535851937790883648493
 *   = 0x1000000000000000 00000000000000 14def9dea2f79cd6 5812631a5cf5d3ed
 *
 * Used for Ed25519 signature arithmetic: s = r + k*a (mod L)
 *
 * Scalars are stored as 32-byte little-endian arrays.
 */

#ifndef CFX_SC25519_H
#define CFX_SC25519_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void cfx_sc25519_reduce(uint8_t out[32], const uint8_t in[64]);  /* reduce 64-byte to 32-byte mod L */
void cfx_sc25519_muladd(uint8_t s[32], const uint8_t a[32], const uint8_t b[32], const uint8_t c[32]);  /* s = (a*b + c) mod L */

#ifdef __cplusplus
}
#endif

#endif /* CFX_SC25519_H */
