/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Dynamic memory backend: deallocation
 *
 * Frees the limb array and resets the structure.
 */

#include "../../mem_backend.h"

#include <stdlib.h>

void cfx_big_free(cfx_big_t* b)
{
    b->n = 0;
    b->cap = 0;
    if (b->limb) free(b->limb);
    b->limb = NULL;
}
