/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Dynamic memory backend: initialization
 *
 * Sets big integer to empty state with no allocated memory.
 */

#include "../../mem_backend.h"

#include <string.h>

void cfx_big_init(cfx_big_t* b)
{
    memset(b, 0, sizeof(*b));
}
