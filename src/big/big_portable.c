/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Portable (base) backend for big integer operations.
 *
 * This is the reference implementation using standard C99. It serves as:
 *   - The fallback for platforms without optimized backends
 *   - The reference for correctness testing
 *   - The baseline for performance comparison
 *
 * Self-guarding: compiles to empty translation unit if another backend is active.
 */

#if !defined(CFX_TARGET_X86_64_BMI2) && \
    !defined(CFX_TARGET_X86_64_AVX2) && \
    !defined(CFX_TARGET_X86_64_AVX512) && \
    !defined(CFX_TARGET_ARM_CORTEX_M4) && \
    !defined(CFX_TARGET_AARCH64_NEON)

#include "big_backend.h"

/* ============================================================================
 * MULTIPLICATION
 * ============================================================================ */

/*
 * Limb-level schoolbook multiplication
 *
 * out[0..na+nb-1] = a[0..na-1] * b[0..nb-1]
 *
 * Classic O(n*m) schoolbook algorithm with accumulator for carries.
 */
void cfx_big_mul_limbs_impl(cfx_limb_t* out,
                            const cfx_limb_t* a, size_t na,
                            const cfx_limb_t* b, size_t nb)
{
    /* Zero the output buffer */
    memset(out, 0, (na + nb) * sizeof(cfx_limb_t));

    for (size_t i = 0; i < nb; ++i) {
        cfx_acc_t carry = 0;
        cfx_acc_t bi = (cfx_acc_t)b[i];

        size_t k = i;
        for (size_t j = 0; j < na; ++j, ++k) {
            cfx_acc_t s = (cfx_acc_t)out[k]
                        + bi * (cfx_acc_t)a[j]
                        + carry;
            out[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
        }

        /* Propagate remaining carry */
        while (carry) {
            cfx_acc_t s = (cfx_acc_t)out[k] + carry;
            out[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
            ++k;
        }
    }
}

/*
 * Big integer multiplication: out = a * b
 *
 * Handles memory allocation and trimming around the limb-level operation.
 */
void cfx_big_mul_impl(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b)
{
    size_t na = a->n;
    size_t nb = b->n;

    /* Handle aliasing: if out overlaps with a or b, use temporary */
    cfx_big_t tmp;
    cfx_big_t* result = out;

    if (out == a || out == b) {
        cfx_big_init(&tmp);
        result = &tmp;
    }

    cfx_big_reserve(result, na + nb);
    result->n = na + nb;

    cfx_big_mul_limbs_impl(result->limb, a->limb, na, b->limb, nb);

    /* Trim leading zeros */
    while (result->n && result->limb[result->n - 1] == 0) {
        --result->n;
    }

    /* Move result to out if we used a temporary */
    if (result != out) {
        cfx_big_move(out, &tmp);
    }
}

/* ============================================================================
 * ADDITION / SUBTRACTION
 * ============================================================================ */

/*
 * Limb-level addition: dst += src
 *
 * Adds n limbs from src to dst, propagating carry.
 * Returns the final carry (0 or 1).
 */
cfx_limb_t cfx_big_add_limbs_impl(cfx_limb_t* dst,
                                  const cfx_limb_t* src,
                                  size_t n)
{
    cfx_limb_t carry = 0;

    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t s = (cfx_acc_t)dst[i] + (cfx_acc_t)src[i] + carry;
        dst[i] = (cfx_limb_t)s;
        carry = (cfx_limb_t)(s >> CFX_LIMB_BITS);
    }

    return carry;
}

/*
 * Limb-level subtraction: dst -= src
 *
 * Subtracts n limbs of src from dst, propagating borrow.
 * Returns the final borrow (0 or 1).
 *
 * Caller must ensure dst >= src mathematically (no underflow).
 */
cfx_limb_t cfx_big_sub_limbs_impl(cfx_limb_t* dst,
                                  const cfx_limb_t* src,
                                  size_t n)
{
    cfx_limb_t borrow = 0;

    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t s = (cfx_acc_t)src[i] + borrow;
        cfx_limb_t di = dst[i];
        dst[i] = di - (cfx_limb_t)s;
        borrow = (di < s) ? 1 : 0;
    }

    return borrow;
}

/* ============================================================================
 * MONTGOMERY MULTIPLICATION
 * ============================================================================ */

/*
 * Montgomery multiplication core: T = a * b * R^{-1} mod n
 *
 * CIOS (Coarsely Integrated Operand Scanning) algorithm.
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

#endif /* !CFX_TARGET_* */
