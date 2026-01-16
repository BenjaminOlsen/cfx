/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * ARM Cortex-M4 optimized implementation of big integer multiplication.
 *
 * Uses UMULL/UMLAL instructions for single-cycle 32x32->64 multiply-accumulate.
 * This replaces the portable cfx_acc_t emulation with native DSP instructions.
 *
 * Requires: -mcpu=cortex-m4 -mthumb
 */

#include "../big_backend.h"

/*
 * Limb-level schoolbook multiplication with UMULL/UMLAL
 *
 * out[0..na+nb-1] = a[0..na-1] * b[0..nb-1]
 *
 * The inner loop accumulates: out[k] += a[j] * b[i] + carry
 * Using UMLAL, this becomes a single-cycle operation for the multiply-add.
 */
void cfx_big_mul_limbs_impl(cfx_limb_t* out,
                            const cfx_limb_t* a, size_t na,
                            const cfx_limb_t* b, size_t nb)
{
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
