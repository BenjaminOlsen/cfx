/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ge25519.h - Group element (point) operations on Ed25519 curve
 *
 * Ed25519 uses the twisted Edwards curve: -x² + y² = 1 + d·x²·y² (mod p)
 * where p = 2^255 - 19 and d = -121665/121666 mod p.
 *
 * Points are represented in extended coordinates: (X : Y : Z : T)
 * where x = X/Z, y = Y/Z, T = X·Y/Z
 *
 * All operations are constant-time.
 */

#ifndef CFX_GE25519_H
#define CFX_GE25519_H

#include "fe25519.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Group element in extended coordinates.
 * Represents a point (x, y) where x = X/Z, y = Y/Z, T = X·Y/Z
 */
typedef struct {
    fe25519_t X;
    fe25519_t Y;
    fe25519_t Z;
    fe25519_t T;
} ge25519_t;

/*
 * Precomputed point for faster addition.
 * Stores (y-x, y+x, 2*d*x*y) in affine form.
 */
typedef struct {
    fe25519_t yplusx;
    fe25519_t yminusx;
    fe25519_t xy2d;
} ge25519_precomp_t;

/*
 * Cached point for readdition.
 * Stores (Y-X, Y+X, Z, 2*d*T) in projective form.
 */
typedef struct {
    fe25519_t YplusX;
    fe25519_t YminusX;
    fe25519_t Z;
    fe25519_t T2d;
} ge25519_cached_t;

extern const fe25519_t cfx_ge25519_d;        /* curve constant d = -121665/121666 mod p */
extern const fe25519_t cfx_ge25519_2d;       /* 2*d for faster doubling */
extern const fe25519_t cfx_ge25519_sqrtm1;   /* sqrt(-1) mod p */
extern const ge25519_t cfx_ge25519_base;     /* basepoint B (generator) */

void cfx_ge25519_0(ge25519_t *r);            /* identity (neutral element): (0, 1, 1, 0) */
void cfx_ge25519_copy(ge25519_t *r, const ge25519_t *p);
void cfx_ge25519_add_cached(ge25519_t *r, const ge25519_t *p, const ge25519_cached_t *q);  /* r = p + q */
void cfx_ge25519_sub_cached(ge25519_t *r, const ge25519_t *p, const ge25519_cached_t *q);  /* r = p - q */
void cfx_ge25519_add_precomp(ge25519_t *r, const ge25519_t *p, const ge25519_precomp_t *q); /* r = p + q (affine q) */
void cfx_ge25519_sub_precomp(ge25519_t *r, const ge25519_t *p, const ge25519_precomp_t *q); /* r = p - q (affine q) */
void cfx_ge25519_double(ge25519_t *r, const ge25519_t *p);    /* r = 2*p */
void cfx_ge25519_to_cached(ge25519_cached_t *r, const ge25519_t *p);
void cfx_ge25519_scalarmult(ge25519_t *r, const uint8_t s[32], const ge25519_t *p);       /* r = [s]p */
void cfx_ge25519_scalarmult_base(ge25519_t *r, const uint8_t s[32]);                       /* r = [s]B */
void cfx_ge25519_pack(uint8_t out[32], const ge25519_t *p);   /* pack point to 32 bytes */
int cfx_ge25519_unpack(ge25519_t *p, const uint8_t in[32]);   /* unpack, returns 0 on success */
int cfx_ge25519_is_on_curve(const ge25519_t *p);
int cfx_ge25519_is_identity(const ge25519_t *p);
void cfx_ge25519_cmov(ge25519_t *r, const ge25519_t *p, unsigned int b);  /* conditional move */
void cfx_ge25519_cneg(ge25519_t *p, unsigned int b);          /* conditional negate */
void cfx_ge25519_neg(ge25519_t *r, const ge25519_t *p);       /* r = -p */

#ifdef __cplusplus
}
#endif

#endif /* CFX_GE25519_H */
