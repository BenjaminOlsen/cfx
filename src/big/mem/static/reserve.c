/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Static memory backend: capacity reservation
 *
 * If the big has no buffer (cap=0), grab one from the pool.
 * Otherwise, validate that the existing buffer has sufficient capacity.
 */

#include "../../mem_backend.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* External declarations - defined in init.c */
extern cfx_limb_t g_static_buffers[CFX_STATIC_POOL_SIZE][CFX_STATIC_LIMBS];
extern int g_static_in_use[CFX_STATIC_POOL_SIZE];

void cfx_big_reserve(cfx_big_t* b, size_t need)
{
    /* If big has no buffer, grab one from the pool (like dynamic mode's realloc(NULL)) */
    if (b->cap == 0) {
        for (int i = 0; i < CFX_STATIC_POOL_SIZE; ++i) {
            if (!g_static_in_use[i]) {
                g_static_in_use[i] = 1;
                b->limb = g_static_buffers[i];
                b->cap = CFX_STATIC_LIMBS;
                memset(b->limb, 0, CFX_STATIC_LIMBS * sizeof(cfx_limb_t));
                break;
            }
        }
        if (b->cap == 0) {
            fprintf(stderr, "cfx_big_reserve: static buffer pool exhausted "
                    "(max %d concurrent big integers)\n", CFX_STATIC_POOL_SIZE);
            abort();
        }
    }

    /* Validate that we have enough capacity */
    if (need > b->cap) {
        fprintf(stderr, "cfx_big_reserve: requested %zu limbs exceeds static capacity %zu "
                "(max %d limbs per big integer)\n", need, b->cap, CFX_STATIC_LIMBS);
        abort();
    }
}
