/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * x86-64 BMI2-optimized implementation of Montgomery multiplication.
 *
 * Uses _mulx_u64 for widening multiply and _addcarry_u64 for
 * addition with carry flag in the CIOS algorithm.
 *
 * Requires: -mbmi2 compiler flag
 */

#include "../big_backend.h"

#include <immintrin.h>

/*
 * Montgomery multiplication core: T = a * b * R^{-1} mod n
 *
 * CIOS algorithm with BMI2 intrinsics for better performance.
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
        const unsigned long long bi = (i < b_n) ? b[i] : 0;

        /* === T += a * b[i] === */
        unsigned long long carry = 0;
        for (size_t j = 0; j < k; ++j) {
            const unsigned long long aj = (j < a_n) ? a[j] : 0;
            unsigned long long hi;
            unsigned long long lo = _mulx_u64(aj, bi, &hi);

            /* Accumulate: T[j] += lo + carry */
            unsigned char c = 0;
            c = _addcarry_u64(c, lo, T[j], &T[j]);
            c = _addcarry_u64(c, hi, carry, &carry);
            carry += c;
        }
        /* Propagate carry to T[k] and T[k+1] */
        {
            unsigned char c = 0;
            c = _addcarry_u64(c, T[k], carry, &T[k]);
            (void)_addcarry_u64(c, T[k + 1], 0, &T[k + 1]);
        }

        /* === m = (T[0] * n0inv) mod 2^64 === */
        const unsigned long long m = T[0] * n0inv;

        /* === T += m * n === */
        carry = 0;
        for (size_t j = 0; j < k; ++j) {
            unsigned long long hi;
            unsigned long long lo = _mulx_u64(m, n[j], &hi);

            unsigned char c = 0;
            c = _addcarry_u64(c, lo, T[j], &T[j]);
            c = _addcarry_u64(c, hi, carry, &carry);
            carry += c;
        }
        /* Propagate carry to T[k] and T[k+1] */
        {
            unsigned char c = 0;
            c = _addcarry_u64(c, T[k], carry, &T[k]);
            (void)_addcarry_u64(c, T[k + 1], 0, &T[k + 1]);
        }

        /* === T >>= 64: shift down by one limb (drop T[0]) === */
        memmove(&T[0], &T[1], (k + 1) * sizeof(cfx_limb_t));
        T[k + 1] = 0;  /* clear scratch for next iteration */
    }
}
