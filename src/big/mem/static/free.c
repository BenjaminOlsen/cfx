/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Static memory backend: deallocation
 *
 * Returns buffer to the pool for reuse.
 */

#include "../../mem_backend.h"

#include <string.h>

/* External declarations - defined in init.c */
extern cfx_limb_t g_static_buffers[CFX_STATIC_POOL_SIZE][CFX_STATIC_LIMBS];
extern int g_static_in_use[CFX_STATIC_POOL_SIZE];

void cfx_big_free(cfx_big_t* b)
{
    if (!b->limb) return;

    /* Find which buffer this is and mark it free */
    for (int i = 0; i < CFX_STATIC_POOL_SIZE; ++i) {
        if (b->limb == g_static_buffers[i]) {
            g_static_in_use[i] = 0;
            break;
        }
    }

    /* Reset the structure */
    b->limb = NULL;
    b->n = 0;
    b->cap = 0;
}
