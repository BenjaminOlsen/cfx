/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ARM Cortex-M4 optimized implementation of big integer addition/subtraction.
 *
 * Uses ADDS/ADCS (add with carry) and SUBS/SBCS (subtract with borrow)
 * instructions for efficient carry propagation using the ARM flags register.
 *
 * Requires: -mcpu=cortex-m4 -mthumb
 */

#include "../big_backend.h"

/*
 * Limb-level addition: dst += src
 *
 * Uses ARM ADDS/ADCS instructions to chain carries through the flags register.
 * This is significantly faster than extracting carries from 64-bit results.
 *
 * Returns the final carry (0 or 1).
 */
cfx_limb_t cfx_big_add_limbs_impl(cfx_limb_t* dst,
                                  const cfx_limb_t* src,
                                  size_t n)
{
    if (n == 0) {
        return 0;
    }

    uint32_t carry = 0;

    /*
     * Process limbs using ADDS/ADCS chain.
     *
     * ADDS sets the carry flag if overflow occurs.
     * ADCS adds with the carry flag and updates it.
     *
     * We process in a loop, maintaining carry between iterations
     * by using ADDS with explicit carry handling.
     */
    for (size_t i = 0; i < n; ++i) {
        uint32_t d = dst[i];
        uint32_t s = src[i];
        uint32_t result;

        /*
         * result = d + s + carry
         * carry_out = 1 if overflow
         *
         * We use inline asm to access the flags directly.
         * The sequence: ADDS then ADC #0 captures the full carry chain.
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
 *
 * Returns the final borrow (0 or 1).
 * Caller must ensure dst >= src mathematically (no underflow).
 */
cfx_limb_t cfx_big_sub_limbs_impl(cfx_limb_t* dst,
                                  const cfx_limb_t* src,
                                  size_t n)
{
    if (n == 0) {
        return 0;
    }

    uint32_t borrow = 0;

    /*
     * Process limbs using SUBS/SBCS chain.
     *
     * ARM uses an inverted carry for subtraction:
     * - SUBS clears C if borrow needed
     * - SBCS subtracts (1 - C), so we track borrow explicitly
     *
     * We compute: dst[i] = dst[i] - src[i] - borrow
     */
    for (size_t i = 0; i < n; ++i) {
        uint32_t d = dst[i];
        uint32_t s = src[i];
        uint32_t result;

        /*
         * Compute d - s - borrow with borrow propagation.
         *
         * First compute d - s, check for borrow.
         * Then subtract incoming borrow, check again.
         */
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
