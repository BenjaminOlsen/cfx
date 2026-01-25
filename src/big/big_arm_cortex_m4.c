/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ARM Cortex-M4 optimized backend for big integer operations.
 *
 * Uses ARM DSP instructions (UMULL/UMLAL, ADDS/ADCS, SUBS/SBCS) for
 * efficient 32-bit limb operations. The Cortex-M4 has single-cycle
 * multiply-accumulate which makes this significantly faster than
 * emulating 64-bit operations.
 *
 * Requires: -mcpu=cortex-m4 -mthumb
 *
 * Self-guarding: only compiles when CFX_TARGET_ARM_CORTEX_M4 is defined.
 */

#ifdef CFX_TARGET_ARM_CORTEX_M4

#include "big_backend.h"

/* ============================================================================
 * MULTIPLICATION
 * ============================================================================ */

/*
 * Limb-level schoolbook multiplication with UMULL/UMLAL
 *
 * out[0..na+nb-1] = a[0..na-1] * b[0..nb-1]
 *
 * The inner loop accumulates: out[k] += a[j] * b[i] + carry
 * Using UMLAL, this becomes a single-cycle operation for the multiply-add.
 */
void cfx_big_mul_limbs_impl(cfx_limb_t *out,
    const cfx_limb_t *a, size_t na,
    const cfx_limb_t *b, size_t nb){
    /* Zero the output buffer */
    memset(out, 0, (na + nb) * sizeof(cfx_limb_t));

    for (size_t i = 0; i < nb; ++i) {
        const uint32_t bi = b[i];
        uint32_t carry = 0;

        size_t k = i;
        for (size_t j = 0; j < na; ++j, ++k) {
            /*
             * Compute: {carry:out[k]} = out[k] + a[j] * bi + carry
             *
             * ARM Cortex-M4 UMLAL does: Rd += Rm * Rn (64-bit accumulate)
             * We use two registers (lo, hi) as a 64-bit accumulator.
             */
            uint32_t lo = out[k];
            uint32_t hi = carry;

            /* UMLAL: {hi:lo} += a[j] * bi */
            __asm__ volatile (
                "umlal %[lo], %[hi], %[aj], %[bi]"
                : [lo] "+r" (lo), [hi] "+r" (hi)
                : [aj] "r" (a[j]), [bi] "r" (bi)
                );

            out[k] = lo;
            carry = hi;
        }

        /* Propagate remaining carry through higher limbs */
        while (carry) {
            uint32_t old = out[k];
            uint32_t sum = old + carry;
            out[k] = sum;
            carry = (sum < old) ? 1 : 0;  /* Carry if overflow */
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

/* ============================================================================
 * ADDITION / SUBTRACTION
 * ============================================================================ */

/*
 * Limb-level addition: dst += src
 *
 * Uses ARM ADDS/ADCS instructions to chain carries through the flags register.
 * Returns the final carry (0 or 1).
 */
cfx_limb_t cfx_big_add_limbs_impl(cfx_limb_t *dst,
    const cfx_limb_t *src,
    size_t n){
    if (n == 0) {
        return 0;
    }

    uint32_t carry = 0;

    /*
     * Process limbs using ADDS/ADC chain.
     * We maintain carry between iterations using explicit carry handling.
     */
    for (size_t i = 0; i < n; ++i) {
        uint32_t d = dst[i];
        uint32_t s = src[i];
        uint32_t result;

        /*
         * result = d + s + carry
         * carry_out = 1 if overflow
         */
        __asm__ volatile (
            "adds %[res], %[d], %[s]\n\t"      /* res = d + s, set C */
            "adc  %[cout], %[cin], #0\n\t"    /* cout = cin + C */
            "adds %[res], %[res], %[cin]\n\t" /* res += cin, set C */
            "adc  %[cout], %[cout], #0"       /* cout += C */
            : [res] "=&r" (result), [cout] "=&r" (carry)
            : [d] "r" (d), [s] "r" (s), [cin] "r" (carry)
            : "cc"
            );

        dst[i] = result;
    }

    return carry;
}

/*
 * Limb-level subtraction: dst -= src
 *
 * Uses ARM SUBS/SBCS instructions for borrow propagation.
 * Returns the final borrow (0 or 1).
 */
cfx_limb_t cfx_big_sub_limbs_impl(cfx_limb_t *dst,
    const cfx_limb_t *src,
    size_t n){
    if (n == 0) {
        return 0;
    }

    uint32_t borrow = 0;

    /*
     * Process limbs using SUBS/SBC chain.
     * We compute: dst[i] = dst[i] - src[i] - borrow
     */
    for (size_t i = 0; i < n; ++i) {
        uint32_t d = dst[i];
        uint32_t s = src[i];
        uint32_t result;

        __asm__ volatile (
            "subs %[res], %[d], %[s]\n\t"        /* res = d - s, C = no borrow */
            "sbc  %[bout], %[bin], %[bin]\n\t"  /* bout = bin - bin - !C = -!C = borrow */
            "subs %[res], %[res], %[bin]\n\t"   /* res -= incoming borrow */
            "sbc  %[bout], %[bout], #0\n\t"     /* bout -= !C (borrow) */
            "rsb  %[bout], %[bout], #0"         /* bout = -bout (invert) */
            : [res] "=&r" (result), [bout] "=&r" (borrow)
            : [d] "r" (d), [s] "r" (s), [bin] "r" (borrow)
            : "cc"
            );

        dst[i] = result;
    }

    return borrow;
}

/* ============================================================================
 * MONTGOMERY MULTIPLICATION
 * ============================================================================ */

/*
 * Montgomery multiplication core: T = a * b * R^{-1} mod n
 *
 * CIOS algorithm with ARM DSP instructions for the multiply-accumulate chains.
 */
void cfx_big_mont_mul_impl(cfx_limb_t *T,
    const cfx_limb_t *a, size_t a_n,
    const cfx_limb_t *b, size_t b_n,
    const cfx_limb_t *n,
    cfx_limb_t n0inv,
    size_t k){
    for (size_t i = 0; i < k; ++i) {
        /* Get b[i], zero-pad if beyond actual length */
        const uint32_t bi = (i < b_n) ? b[i] : 0;

        /*
         * Step 1: T += a * b[i]
         */
        uint32_t carry_lo = 0;
        uint32_t carry_hi = 0;

        for (size_t j = 0; j < k; ++j) {
            const uint32_t aj = (j < a_n) ? a[j] : 0;

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
         * Step 2: m = (T[0] * n0inv) mod 2^32
         */
        const uint32_t m = T[0] * n0inv;

        /*
         * Step 3: T += m * n
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
         */
        memmove(&T[0], &T[1], (k + 1) * sizeof(cfx_limb_t));
        T[k + 1] = 0;
    }
}

#endif /* CFX_TARGET_ARM_CORTEX_M4 */
