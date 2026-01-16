/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Static memory backend: initialization
 *
 * Allocates from a fixed pool of static buffers. No heap allocation.
 * Suitable for embedded systems and environments without malloc.
 *
 * Configuration:
 *   CFX_STATIC_LIMBS - limbs per buffer (default: 64)
 *   CFX_STATIC_POOL_SIZE - number of buffers in pool (default: 32)
 */

#include "../../mem_backend.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Static buffer pool - non-static so free.c can access them */
cfx_limb_t g_static_buffers[CFX_STATIC_POOL_SIZE][CFX_STATIC_LIMBS];
int g_static_in_use[CFX_STATIC_POOL_SIZE];

void cfx_big_init(cfx_big_t* b)
{
    /* Find a free buffer in the pool */
    for (int i = 0; i < CFX_STATIC_POOL_SIZE; ++i) {
        if (!g_static_in_use[i]) {
            g_static_in_use[i] = 1;
            b->limb = g_static_buffers[i];
            b->n = 0;
            b->cap = CFX_STATIC_LIMBS;
            memset(b->limb, 0, CFX_STATIC_LIMBS * sizeof(cfx_limb_t));
            return;
        }
    }

    /* Pool exhausted */
    fprintf(stderr, "cfx_big_init: static buffer pool exhausted "
            "(max %d concurrent big integers)\n", CFX_STATIC_POOL_SIZE);
    abort();
}
