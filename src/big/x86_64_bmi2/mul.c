/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * x86-64 BMI2-optimized implementation of big integer multiplication.
 *
 * Uses _mulx_u64 for widening multiply. This avoids the portable
 * cfx_acc_t abstraction for slightly better codegen on x86-64.
 *
 * Requires: -mbmi2 compiler flag
 */

#include "../big_backend.h"

#include <immintrin.h>

/*
 * Limb-level schoolbook multiplication with BMI2 intrinsics
 *
 * out[0..na+nb-1] = a[0..na-1] * b[0..nb-1]
 *
 * Uses _mulx_u64 for the widening multiply, but keeps the same
 * accumulator logic as the generic implementation for correctness.
 */
void cfx_big_mul_limbs_impl(cfx_limb_t* out,
                            const cfx_limb_t* a, size_t na,
                            const cfx_limb_t* b, size_t nb)
{
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
