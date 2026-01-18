/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Dynamic memory backend for big integers.
 *
 * Uses malloc/realloc/free for memory management. This is the default
 * mode for desktop and server applications where heap allocation is fine.
 *
 * Self-guarding: only compiles when CFX_MEMORY_STATIC is NOT defined.
 */

#ifndef CFX_MEMORY_STATIC

#include "mem_backend.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ============================================================================
 * INITIALIZATION
 * ============================================================================ */

/*
 * Initialize big integer to empty state with no allocated memory.
 */
void cfx_big_init(cfx_big_t* b)
{
    memset(b, 0, sizeof(*b));
}

/* ============================================================================
 * CAPACITY RESERVATION
 * ============================================================================ */

/*
 * Ensure capacity for at least 'need' limbs.
 *
 * Uses exponential growth strategy (doubling) to amortize allocation cost.
 * New memory is zeroed after allocation.
 */
void cfx_big_reserve(cfx_big_t* b, size_t need)
{
    if (need <= b->cap) return;

    size_t old_cap = b->cap;
    size_t new_cap = old_cap ? old_cap : 32;

    /* Exponential growth with overflow protection */
    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) {
            new_cap = need;
            break;
        }
        new_cap *= 2;
    }

    /* Check for allocation size overflow */
    if (new_cap > SIZE_MAX / sizeof(cfx_limb_t)) {
        fprintf(stderr, "cfx_big_reserve: allocation size overflow\n");
        abort();
    }

    size_t new_bytes = new_cap * sizeof(cfx_limb_t);

    void* tmp = realloc(b->limb, new_bytes);
    if (!tmp) {
        fprintf(stderr, "cfx_big_reserve: out of memory (requested %zu bytes)\n", new_bytes);
        abort();
    }

    b->limb = (cfx_limb_t*)tmp;

    /* Zero the newly added region */
    if (new_cap > old_cap) {
        size_t add = new_cap - old_cap;
        memset(b->limb + old_cap, 0, add * sizeof(cfx_limb_t));
    }

    b->cap = new_cap;
}

/* ============================================================================
 * DEALLOCATION
 * ============================================================================ */

/*
 * Free the limb array and reset the structure.
 */
void cfx_big_free(cfx_big_t* b)
{
    b->n = 0;
    b->cap = 0;
    if (b->limb) free(b->limb);
    b->limb = NULL;
}

#endif /* !CFX_MEMORY_STATIC */
