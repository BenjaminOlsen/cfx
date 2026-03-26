/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * x86-64 BMI2-optimized backend for big integer operations.
 *
 * Uses BMI2 intrinsics (_mulx_u64, _addcarry_u64) for widening multiply
 * and carry-chain operations. Provides significant speedup over portable
 * code for multiplication and Montgomery reduction.
 *
 * Requires: -mbmi2 compiler flag
 *
 * Self-guarding: only compiles when CFX_TARGET_X86_64_BMI2 (or AVX2/AVX512) is defined.
 */

#if defined(CFX_TARGET_X86_64_BMI2) || \
    defined(CFX_TARGET_X86_64_AVX2) || \
    defined(CFX_TARGET_X86_64_AVX512)

#include "big_backend.h"

#include <immintrin.h>

/* MULTIPLICATION */

/*
 * Limb-level schoolbook multiplication with BMI2 intrinsics
 *
 * out[0..na+nb-1] = a[0..na-1] * b[0..nb-1]
 *
 * Uses _mulx_u64 for the widening multiply.
 */
void cfx_big_mul_limbs_impl(cfx_limb_t *out,
    const cfx_limb_t *a, size_t na,
    const cfx_limb_t *b, size_t nb){
    /* Zero the output buffer */
    memset(out, 0, (na + nb) * sizeof(cfx_limb_t));

    for (size_t i = 0; i < nb; ++i) {
        const cfx_limb_t bi = b[i];
        cfx_acc_t carry = 0;

        size_t k = i;
        for (size_t j = 0; j < na; ++j, ++k) {
            /*
             * _mulx_u64 computes 64x64->128 without affecting flags.
             * lo = low 64 bits, hi = high 64 bits
             */
            unsigned long long hi;
            unsigned long long lo = _mulx_u64(a[j], bi, &hi);

            /*
             * Accumulate: out[k] + (hi:lo) + carry
             * Using cfx_acc_t (128-bit) for the full accumulation
             */
            cfx_acc_t product = ((cfx_acc_t)hi << CFX_LIMB_BITS) | lo;
            cfx_acc_t s = (cfx_acc_t)out[k] + product + carry;
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
void cfx_big_mul_impl(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b){
    size_t na = a->n;
    size_t nb = b->n;

    /* Handle aliasing: if out overlaps with a or b, use temporary */
    cfx_big_t tmp;
    cfx_big_t *result = out;

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

/* ADDITION / SUBTRACTION
 *
 * Use portable implementation - BMI2 doesn't help much here.
 * The portable cfx_acc_t approach is already efficient for add/sub.
 */

cfx_limb_t cfx_big_add_limbs_impl(cfx_limb_t *dst,
    const cfx_limb_t *src,
    size_t n){
    cfx_limb_t carry = 0;

    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t s = (cfx_acc_t)dst[i] + (cfx_acc_t)src[i] + carry;
        dst[i] = (cfx_limb_t)s;
        carry = (cfx_limb_t)(s >> CFX_LIMB_BITS);
    }

    return carry;
}

cfx_limb_t cfx_big_sub_limbs_impl(cfx_limb_t *dst,
    const cfx_limb_t *src,
    size_t n){
    cfx_limb_t borrow = 0;

    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t s = (cfx_acc_t)src[i] + borrow;
        cfx_limb_t di = dst[i];
        dst[i] = di - (cfx_limb_t)s;
        borrow = (di < s) ? 1 : 0;
    }

    return borrow;
}

/* MONTGOMERY MULTIPLICATION */

/*
 * Montgomery multiplication core: T = a * b * R^{-1} mod n
 *
 * CIOS algorithm with BMI2 intrinsics for the multiply chains.
 */
void cfx_big_mont_mul_impl(cfx_limb_t *T,
    const cfx_limb_t *a, size_t a_n,
    const cfx_limb_t *b, size_t b_n,
    const cfx_limb_t *n,
    cfx_limb_t n0inv,
    size_t k){
    for (size_t i = 0; i < k; ++i) {
        /* Get b[i], zero-pad if beyond actual length */
        const unsigned long long bi = (i < b_n) ? b[i] : 0;

        /* T += a * b[i] */
        unsigned long long carry = 0;
        for (size_t j = 0; j < k; ++j) {
            const unsigned long long aj = (j < a_n) ? a[j] : 0;
            unsigned long long hi;
            unsigned long long lo = _mulx_u64(aj, bi, &hi);

            /*
             * Accumulate: T[j] += lo + carry_low, carry = hi + overflow
             *
             * We need to add lo AND the low part of carry to T[j].
             * Since carry can be > 2^64 (up to ~2^64 + k), we use 128-bit arithmetic.
             */
            cfx_acc_t sum = (cfx_acc_t)T[j] + lo + carry;
            T[j] = (cfx_limb_t)sum;
            carry = hi + (sum >> CFX_LIMB_BITS);
        }
        /* Propagate carry to T[k] and T[k+1] */
        {
            cfx_acc_t sum = (cfx_acc_t)T[k] + carry;
            T[k] = (cfx_limb_t)sum;
            T[k + 1] += (cfx_limb_t)(sum >> CFX_LIMB_BITS);
        }

        /* m = (T[0] * n0inv) mod 2^64 */
        const unsigned long long m = T[0] * n0inv;

        /* T += m * n */
        carry = 0;
        for (size_t j = 0; j < k; ++j) {
            unsigned long long hi;
            unsigned long long lo = _mulx_u64(m, n[j], &hi);

            cfx_acc_t sum = (cfx_acc_t)T[j] + lo + carry;
            T[j] = (cfx_limb_t)sum;
            carry = hi + (sum >> CFX_LIMB_BITS);
        }
        /* Propagate carry to T[k] and T[k+1] */
        {
            cfx_acc_t sum = (cfx_acc_t)T[k] + carry;
            T[k] = (cfx_limb_t)sum;
            T[k + 1] += (cfx_limb_t)(sum >> CFX_LIMB_BITS);
        }

        /* T >>= 64: shift down by one limb (drop T[0]) */
        memmove(&T[0], &T[1], (k + 1) * sizeof(cfx_limb_t));
        T[k + 1] = 0;  /* clear scratch for next iteration */
    }
}

#endif /* CFX_TARGET_X86_64_BMI2 || CFX_TARGET_X86_64_AVX2 || CFX_TARGET_X86_64_AVX512 */
