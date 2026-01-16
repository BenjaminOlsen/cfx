/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ARM Cortex-M4 optimized implementation of Montgomery multiplication.
 *
 * Uses UMULL/UMLAL instructions for the CIOS algorithm inner loops,
 * providing single-cycle multiply-accumulate operations.
 *
 * Requires: -mcpu=cortex-m4 -mthumb
 */

#include "../big_backend.h"

/*
 * Montgomery multiplication core: T = a * b * R^{-1} mod n
 *
 * CIOS (Coarsely Integrated Operand Scanning) algorithm with ARM DSP
 * instructions for the multiply-accumulate chains.
 *
 * Parameters:
 *   T     - accumulator array with k+2 limbs (pre-zeroed)
 *   a     - first operand limbs
 *   a_n   - number of non-zero limbs in a
 *   b     - second operand limbs
 *   b_n   - number of non-zero limbs in b
 *   n     - modulus limbs (k limbs)
 *   n0inv - -n[0]^{-1} mod 2^32
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
        const uint32_t bi = (i < b_n) ? b[i] : 0;

        /*
         * Step 1: T += a * b[i]
         *
         * Inner loop uses UMLAL for fused multiply-accumulate:
         *   {carry_hi:T[j]} += a[j] * bi
         */
        uint32_t carry_lo = 0;
        uint32_t carry_hi = 0;

        for (size_t j = 0; j < k; ++j) {
            const uint32_t aj = (j < a_n) ? a[j] : 0;

            /* {carry_hi:carry_lo} = T[j] + a[j] * bi + carry */
            uint32_t lo = T[j];
            uint32_t hi = 0;

            /* Add previous carry */
            __asm__ volatile (
                "adds %[lo], %[lo], %[clo]\n\t"
                "adc  %[hi], %[hi], %[chi]"
                : [lo] "+r" (lo), [hi] "+r" (hi)
                : [clo] "r" (carry_lo), [chi] "r" (carry_hi)
                : "cc"
            );

            /* Multiply-accumulate: {hi:lo} += aj * bi */
            __asm__ volatile (
                "umlal %[lo], %[hi], %[aj], %[bi]"
                : [lo] "+r" (lo), [hi] "+r" (hi)
                : [aj] "r" (aj), [bi] "r" (bi)
            );

            T[j] = lo;
            carry_lo = hi;
            carry_hi = 0;  /* UMLAL into 64 bits doesn't overflow for single product */
        }

        /* Propagate carry to T[k] and T[k+1] */
        {
            uint32_t lo = T[k];
            uint32_t hi = T[k + 1];

            __asm__ volatile (
                "adds %[lo], %[lo], %[clo]\n\t"
                "adc  %[hi], %[hi], #0"
                : [lo] "+r" (lo), [hi] "+r" (hi)
                : [clo] "r" (carry_lo)
                : "cc"
            );

            T[k] = lo;
            T[k + 1] = hi;
        }

        /*
         * Step 2: m = (T[0] * n0inv) mod 2^32
         *
         * This is just a single 32-bit multiply (low bits only).
         */
        const uint32_t m = T[0] * n0inv;

        /*
         * Step 3: T += m * n
         *
         * Same UMLAL chain as step 1.
         */
        carry_lo = 0;
        carry_hi = 0;

        for (size_t j = 0; j < k; ++j) {
            uint32_t lo = T[j];
            uint32_t hi = 0;

            /* Add previous carry */
            __asm__ volatile (
                "adds %[lo], %[lo], %[clo]\n\t"
                "adc  %[hi], %[hi], %[chi]"
                : [lo] "+r" (lo), [hi] "+r" (hi)
                : [clo] "r" (carry_lo), [chi] "r" (carry_hi)
                : "cc"
            );

            /* Multiply-accumulate: {hi:lo} += m * n[j] */
            __asm__ volatile (
                "umlal %[lo], %[hi], %[m], %[nj]"
                : [lo] "+r" (lo), [hi] "+r" (hi)
                : [m] "r" (m), [nj] "r" (n[j])
            );

            T[j] = lo;
            carry_lo = hi;
            carry_hi = 0;
        }

        /* Propagate carry to T[k] and T[k+1] */
        {
            uint32_t lo = T[k];
            uint32_t hi = T[k + 1];

            __asm__ volatile (
                "adds %[lo], %[lo], %[clo]\n\t"
                "adc  %[hi], %[hi], #0"
                : [lo] "+r" (lo), [hi] "+r" (hi)
                : [clo] "r" (carry_lo)
                : "cc"
            );

            T[k] = lo;
            T[k + 1] = hi;
        }

        /*
         * Step 4: T >>= 32 (shift down by one limb, drop T[0])
         *
         * After adding m*n, T[0] is guaranteed to be 0 (that's the point
         * of Montgomery reduction). We shift the array down.
         */
        memmove(&T[0], &T[1], (k + 1) * sizeof(cfx_limb_t));
        T[k + 1] = 0;  /* Clear scratch for next iteration */
    }
}
