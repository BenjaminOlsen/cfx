/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * fe25519.h - Field element arithmetic in GF(2^255 - 19)
 *
 * Radix-2^51 representation: 5 limbs of 51 bits each, stored in uint64_t.
 * Each limb has 13 bits of headroom for lazy carry propagation.
 *
 * All operations are constant-time.
 */

#ifndef CFX_FE25519_H
#define CFX_FE25519_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Field element in GF(p) where p = 2^255 - 19.
 * Represented as 5 limbs in radix 2^51.
 */
typedef struct {
    uint64_t v[5];
} fe25519_t;


extern const fe25519_t cfx_fe25519_zero;    /* zero constant */
extern const fe25519_t cfx_fe25519_one;     /* one constant */

void cfx_fe25519_0(fe25519_t *h);   /* set to zero */
void cfx_fe25519_1(fe25519_t *h);   /* set to one */
void cfx_fe25519_copy(fe25519_t *h, const fe25519_t *f);
void cfx_fe25519_add(fe25519_t *h, const fe25519_t *f, const fe25519_t *g); /* h = f + g (no reduction needed) */
void cfx_fe25519_sub(fe25519_t *h, const fe25519_t *f, const fe25519_t *g); /* h = f - g (adds 2p to avoid underflow) */
void cfx_fe25519_neg(fe25519_t *h, const fe25519_t *f); /* h = -f */
void cfx_fe25519_mul(fe25519_t *h, const fe25519_t *f, const fe25519_t *g); /* h = f * g (with reduction mod p) */
void cfx_fe25519_sqr(fe25519_t *h, const fe25519_t *f);
void cfx_fe25519_sqr_n(fe25519_t *h, const fe25519_t *f, unsigned n);   /* h = f^(2^n) (repeated squaring) */
void cfx_fe25519_mul121666(fe25519_t *h, const fe25519_t *f);   /* h = f * 121666 */
void cfx_fe25519_inv(fe25519_t *h, const fe25519_t *f); /* h = 1/f */
void cfx_fe25519_pow22523(fe25519_t *h, const fe25519_t *f);    /* h = f^((p-5)/8)  */
void cfx_fe25519_carry(fe25519_t *h);
void cfx_fe25519_reduce(fe25519_t *h);      /* full reduction to canonical [0, p) range */
void cfx_fe25519_cswap(fe25519_t *f, fe25519_t *g, unsigned int b);
void cfx_fe25519_cmov(fe25519_t *h, const fe25519_t *f, unsigned int b);
int cfx_fe25519_iszero(const fe25519_t *f);
int cfx_fe25519_isnegative(const fe25519_t *f);
int cfx_fe25519_eq(const fe25519_t *f, const fe25519_t *g);
void cfx_fe25519_frombytes(fe25519_t *h, const uint8_t s[32]);  /* deserialize from 32 bytes (little-endian, reduced mod p) */
void cfx_fe25519_tobytes(uint8_t s[32], const fe25519_t *h);    /* serialize to 32 bytes (little-endian, fully reduced) */

#ifdef __cplusplus
}
#endif

#endif /* CFX_FE25519_H */
