/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Poly1305 backend interface
 *
 * This header declares the target-specific functions for Poly1305's
 * inner field arithmetic. The main poly1305.c calls these *_impl functions.
 *
 * The core operation is computing:
 *   d[i] = sum of h[j] * r[k] products (5 terms each)
 *
 * On Cortex-M4, UMULL/UMLAL chains provide significant speedup.
 */

#ifndef CFX_POLY1305_BACKEND_H
#define CFX_POLY1305_BACKEND_H

#include <stdint.h>

/*
 * Poly1305 field multiply-accumulate for a single block.
 *
 * Computes h = (h + m) * r mod (2^130 - 5) where:
 *   - h is the accumulator in 5x26-bit limb form
 *   - m is the input block (decoded into 26-bit limbs)
 *   - r is the key in 26-bit limbs with precomputed s[i] = 5*r[i]
 *
 * Parameters:
 *   h0..h4   - accumulator limbs (in/out)
 *   t0..t4   - message block limbs (already decoded)
 *   r0..r4   - key limbs
 *   s1..s4   - precomputed 5*r[1..4] for reduction
 *
 * The 5-term MAC for each accumulator:
 *   d0 = h0*r0 + h1*s4 + h2*s3 + h3*s2 + h4*s1
 *   d1 = h0*r1 + h1*r0 + h2*s4 + h3*s3 + h4*s2
 *   d2 = h0*r2 + h1*r1 + h2*r0 + h3*s4 + h4*s3
 *   d3 = h0*r3 + h1*r2 + h2*r1 + h3*r0 + h4*s4
 *   d4 = h0*r4 + h1*r3 + h2*r2 + h3*r1 + h4*r0
 */
void cfx_poly1305_block_impl(
    uint32_t *h0, uint32_t *h1, uint32_t *h2, uint32_t *h3, uint32_t *h4,
    uint32_t t0, uint32_t t1, uint32_t t2, uint32_t t3, uint32_t t4,
    uint32_t r0, uint32_t r1, uint32_t r2, uint32_t r3, uint32_t r4,
    uint32_t s1, uint32_t s2, uint32_t s3, uint32_t s4);

#endif /* CFX_POLY1305_BACKEND_H */
