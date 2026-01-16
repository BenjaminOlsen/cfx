/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Generic (portable) implementation of Montgomery multiplication.
 *
 * This implements the CIOS (Coarsely Integrated Operand Scanning)
 * algorithm for Montgomery multiplication.
 */

#include "../big_backend.h"

/*
 * Montgomery multiplication core: T = a * b * R^{-1} mod n
 *
 * CIOS algorithm:
 *   for i = 0..k-1:
 *     T += a * b[i]
 *     m = (T[0] * n0inv) mod 2^LIMB_BITS
 *     T += m * n
 *     T >>= LIMB_BITS   (drop least-significant limb)
 *
 * Parameters:
 *   T     - accumulator array with k+2 limbs (pre-zeroed)
 *   a     - first operand limbs
 *   a_n   - number of non-zero limbs in a
 *   b     - second operand limbs
 *   b_n   - number of non-zero limbs in b
 *   n     - modulus limbs (k limbs)
 *   n0inv - -n[0]^{-1} mod 2^LIMB_BITS
 *   k     - number of limbs in modulus
 *
 * Postconditions:
 *   - T[0..k] contains result (may be >= n, caller does final reduction)
 */
void cfx_big_mont_mul_impl(cfx_limb_t* T,
                           const cfx_limb_t* a, size_t a_n,
                           const cfx_limb_t* b, size_t b_n,
                           const cfx_limb_t* n,
                           cfx_limb_t n0inv,
                           size_t k)
{
    for (size_t i = 0; i < k; ++i) {
        /* Get b[i], zero-pad if beyond actual length */
        const cfx_limb_t bi = (i < b_n) ? b[i] : 0;

        /* T += a * b[i] */
        cfx_acc_t carry = 0;
        for (size_t j = 0; j < k; ++j) {
            const cfx_limb_t aj = (j < a_n) ? a[j] : 0;
            cfx_acc_t sum = (cfx_acc_t)T[j] + (cfx_acc_t)aj * bi + carry;
            T[j] = (cfx_limb_t)sum;
            carry = sum >> CFX_LIMB_BITS;
        }
        cfx_acc_t top = (cfx_acc_t)T[k] + carry;
        T[k] = (cfx_limb_t)top;
        T[k + 1] += (cfx_limb_t)(top >> CFX_LIMB_BITS);

        /* m = (T[0] * n0inv) mod 2^LIMB_BITS */
        const cfx_limb_t m = T[0] * n0inv;

        /* T += m * n */
        cfx_acc_t carry2 = 0;
        for (size_t j = 0; j < k; ++j) {
            cfx_acc_t sum = (cfx_acc_t)T[j] + (cfx_acc_t)m * n[j] + carry2;
            T[j] = (cfx_limb_t)sum;
            carry2 = sum >> CFX_LIMB_BITS;
        }
        cfx_acc_t top2 = (cfx_acc_t)T[k] + carry2;
        T[k] = (cfx_limb_t)top2;
        T[k + 1] += (cfx_limb_t)(top2 >> CFX_LIMB_BITS);

        /* T >>= LIMB_BITS: shift down by one limb (drop T[0]) */
        memmove(&T[0], &T[1], (k + 1) * sizeof(cfx_limb_t));
        T[k + 1] = 0;  /* clear scratch for next iteration */
    }
}
