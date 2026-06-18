/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Static memory backend for big integers.
 *
 * Uses a fixed pool of static buffers instead of heap allocation.
 *
 * Configuration (via CMake or compile flags):
 *   CFX_STATIC_LIMBS     - limbs per buffer (default: 256)
 *   CFX_STATIC_POOL_SIZE - number of buffers in pool (default: 1024)
 *
 * Only compiles when CFX_MEMORY_STATIC is defined.
 */

#ifdef CFX_MEMORY_STATIC

#include "mem_backend.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* STATIC BUFFER POOL */

/*
 * The static buffer pool.
 *
 * Each slot is a fixed-size array of limbs. The in_use array tracks
 * which slots are currently allocated.
 */
static cfx_limb_t g_static_buffers[CFX_STATIC_POOL_SIZE][CFX_STATIC_LIMBS];
static int g_static_in_use[CFX_STATIC_POOL_SIZE];

/* INITIALIZATION */

/*
 * Initialize big integer by allocating a buffer from the pool.
 */
void cfx_big_init(cfx_big_t *b){
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

/* CAPACITY RESERVATION */

/*
 * Ensure capacity for at least 'need' limbs.
 *
 * In static mode, this validates that the existing buffer is large enough.
 * If the big has no buffer yet (cap=0), allocate one from the pool.
 */
void cfx_big_reserve(cfx_big_t *b, size_t need){
    /* If big has no buffer, grab one from the pool */
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

/* DEALLOCATION */

/*
 * Return buffer to the pool for reuse.
 */
void cfx_big_free(cfx_big_t *b){
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

#endif /* CFX_MEMORY_STATIC */
