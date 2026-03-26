/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_BIG_BACKEND_H
#define CFX_BIG_BACKEND_H

/*
 * Big Integer Backend Interface
 *
 * This header defines the internal interface for backend implementations.
 * The facade (big.c) calls these *_impl() functions after handling edge
 * cases and validation.
 *
 * Implementations MUST be:
 *   - Bit-compatible with base/ for deterministic operations
 *   - Constant-time where noted (for crypto use)
 *
 * Backend selection is done at compile time via CMake. Each backend
 * implements these functions in its own directory (e.g., base/, x86_64_bmi2/).
 */

#include "cfx/big.h"
#include "cfx/arith.h"

#include <stddef.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/* MULTIPLICATION BACKENDS */

/*
 * Core multiplication: out = a * b
 *
 * Preconditions:
 *   - a, b are non-NULL and initialized
 *   - Neither a nor b is zero (caller handles this case)
 *   - out may alias a or b
 *
 * Postconditions:
 *   - out contains the product, properly trimmed
 *   - out->n <= a->n + b->n
 */
void cfx_big_mul_impl(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);

/*
 * Limb-level schoolbook multiplication (helper for backends)
 *
 * Computes: out[0..na+nb-1] = a[0..na-1] * b[0..nb-1]
 *
 * Preconditions:
 *   - out has space for na+nb limbs, pre-zeroed
 *   - a, b are non-NULL
 *   - na, nb > 0
 *
 * Note: This is the raw limb-level operation; doesn't trim or handle big_t.
 */
void cfx_big_mul_limbs_impl(cfx_limb_t *out,
    const cfx_limb_t *a, size_t na,
    const cfx_limb_t *b, size_t nb);

/* ADDITION/SUBTRACTION BACKENDS */

/*
 * Limb-level addition: dst += src
 *
 * Adds src[0..n-1] to dst[0..n-1], propagating carry.
 *
 * Preconditions:
 *   - dst has at least n limbs
 *   - src has at least n limbs
 *
 * Returns:
 *   - Final carry (0 or 1)
 */
cfx_limb_t cfx_big_add_limbs_impl(cfx_limb_t *dst,
    const cfx_limb_t *src,
    size_t n);

/*
 * Limb-level subtraction: dst -= src
 *
 * Subtracts src[0..n-1] from dst[0..n-1], propagating borrow.
 *
 * Preconditions:
 *   - dst has at least n limbs
 *   - src has at least n limbs
 *   - Mathematically, dst >= src (no underflow)
 *
 * Returns:
 *   - Final borrow (0 or 1; should be 0 if preconditions met)
 */
cfx_limb_t cfx_big_sub_limbs_impl(cfx_limb_t *dst,
    const cfx_limb_t *src,
    size_t n);

/* MONTGOMERY MULTIPLICATION BACKEND */

/*
 * Montgomery multiplication: T = a * b * R^{-1} mod n
 *
 * Uses CIOS (Coarsely Integrated Operand Scanning) algorithm.
 *
 * Parameters:
 *   T     - output limb array, pre-allocated with k+2 limbs (zeroed)
 *   a     - first operand limbs (at least k limbs, or zero-padded)
 *   b     - second operand limbs (at least k limbs, or zero-padded)
 *   n     - modulus limbs (k limbs)
 *   n0inv - -n[0]^{-1} mod 2^LIMB_BITS
 *   k     - number of limbs in modulus
 *   a_n   - actual number of limbs in a (for zero-padding)
 *   b_n   - actual number of limbs in b (for zero-padding)
 *
 * Preconditions:
 *   - n is odd (n[0] & 1 == 1)
 *   - All inputs properly allocated
 *   - T is pre-zeroed with k+2 limbs
 *
 * Postconditions:
 *   - T[0..k-1] contains result (may need final subtraction)
 *   - Result is in range [0, 2n), caller does final normalization
 *
 * Note: This is the core loop; caller handles big_t wrapper and normalization.
 */
void cfx_big_mont_mul_impl(cfx_limb_t *T,
    const cfx_limb_t *a, size_t a_n,
    const cfx_limb_t *b, size_t b_n,
    const cfx_limb_t *n,
    cfx_limb_t n0inv,
    size_t k);

/* BACKEND IDENTIFICATION */

/*
 * Backend name string for debugging/logging.
 * Defined in each backend's source file.
 */
#if defined(CFX_TARGET_X86_64_BMI2) || defined(CFX_TARGET_X86_64_AVX2) || defined(CFX_TARGET_X86_64_AVX512)
    #define CFX_BIG_BACKEND_NAME "x86_64_bmi2"
#elif defined(CFX_TARGET_ARM_CORTEX_M4)
    #define CFX_BIG_BACKEND_NAME "arm_cortex_m4"
#elif defined(CFX_TARGET_AARCH64_NEON)
    #define CFX_BIG_BACKEND_NAME "aarch64_neon"
#else
    #define CFX_BIG_BACKEND_NAME "base"
#endif

#ifdef __cplusplus
}
#endif

#endif /* CFX_BIG_BACKEND_H */
