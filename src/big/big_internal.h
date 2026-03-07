/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_BIG_INTERNAL_H
#define CFX_BIG_INTERNAL_H

#include "cfx/big.h"
#include "cfx/arith.h"
#include "cfx/macros.h"

#include "big_backend.h"
#include "mem_backend.h"

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

static inline void cfx_big_trim(cfx_big_t *b) {
    while (b->n && b->limb[b->n-1] == 0) --b->n;
    if (b->n == 0 && b->cap){
        b->limb[0] = 0;
    }                                           /* dont trim last zero */
}

#endif /* CFX_BIG_INTERNAL_H */
