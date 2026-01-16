/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Generic (portable) implementation of big integer addition/subtraction.
 *
 * This is the reference implementation using standard C with the
 * cfx accumulator abstraction for carry/borrow propagation.
 */

#include "../big_backend.h"

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
