/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Generic (portable) implementation of big integer multiplication.
 *
 * This is the reference implementation using standard C with the
 * cfx accumulator abstraction for wide multiplies.
 */

#include "../big_backend.h"

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
