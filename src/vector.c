/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/vector.h"
#include "cfx/macros.h"

#include <stdio.h>
#include <stdlib.h>

void cfx_vec_init(cfx_vec_t *v) {
    v->data = NULL;
    v->size = 0;
    v->cap = 0;
}

int cfx_vec_reserve(cfx_vec_t* v, size_t need) {
    if (need < v->cap) return 0;
    
    const size_t max_size = SIZE_MAX / sizeof(cfx_limb_t);
    if (need > max_size) return -1;

    size_t new_cap = v->cap ? 2 * v->cap : CFX_VEC_MIN_CAP;
    
    while ((new_cap < need) && (new_cap < (max_size / 2))) {
        new_cap *= 2;
    }

    void *p = realloc(v->data, new_cap * sizeof(cfx_limb_t));
    if (!p) return -1;

    v->data = (cfx_limb_t*)p;
    v->cap = new_cap;
    // CFX_PRINT_DBG("reserved cap %zu\n", new_cap);
    return 0;
}

int cfx_vec_resize(cfx_vec_t* v, size_t sz, cfx_limb_t fill) {
    if (sz > v->cap) {
        if (cfx_vec_reserve(v, sz) != 0) return -1;
    }
    if (sz > v->size) {
        for (size_t k = v->size; k < sz; ++k) {
            v->data[k] = fill;
        }
    }
    v->size = sz;
    return 0;
}

int cfx_vec_push(cfx_vec_t* v, cfx_limb_t x) {
    if (v->size == v->cap) {
        if (cfx_vec_reserve(v, v->size + 1) != 0) return -1;
    }
    v->data[v->size] = x;
    ++v->size;
    return 0;
}

void cfx_vec_free(cfx_vec_t* v) {
    free(v->data);
    v->data = NULL;
    v->cap = 0;
    v->size = 0;
}

void cfx_vec_clear(cfx_vec_t* v) {
    v->size = 0;
}

void cfx_vec_print(const cfx_vec_t* v) {
    CFX_PRINT_DBG("v: %p, size: %zu\n", (void*)v, v->size);
    for (size_t i = 0; i < v->size; ++i) {
        CFX_PRINT_DBG(""CFX_PRIuLIMB" ", v->data[i]);
    }
    CFX_PRINT_DBG("\n");
}
