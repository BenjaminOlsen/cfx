/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/big_vec.h"

#include <stdlib.h>
#include <string.h>

static size_t cfx_big_vec_next_cap(size_t cur, size_t need) {
    size_t cap = (cur < CFX_BIG_VEC_MIN_CAP) ? CFX_BIG_VEC_MIN_CAP : cur;
    while (cap < need) {
        /* grow ~2x; avoid overflow */
        size_t next = cap * 2;
        if (next < cap) return need;   /* overflow fallback */
        cap = next;
    }
    return cap;
}

void cfx_big_vec_init(cfx_big_vec_t* v) {
    v->data = NULL;
    v->size = 0;
    v->cap  = 0;
}

void cfx_big_vec_clear(cfx_big_vec_t* v) {
    if (!v || !v->data) {
        if (v) v->size = 0;
        return;
    }
    for (size_t i = 0; i < v->size; i++) {
        cfx_big_free(&v->data[i]);
    }
    v->size = 0;
}

void cfx_big_vec_free(cfx_big_vec_t* v) {
    if (!v) return;
    cfx_big_vec_clear(v);
    free(v->data);
    v->data = NULL;
    v->cap  = 0;
}

int cfx_big_vec_reserve(cfx_big_vec_t* v, size_t need) {
    if (!v) return -1;
    if (need <= v->cap) return 0;

    size_t new_cap = cfx_big_vec_next_cap(v->cap, need);

    cfx_big_t* new_data = (cfx_big_t*)malloc(new_cap * sizeof(cfx_big_t));
    if (!new_data) return -1;

    /* move-construct existing elements into new buffer */
    for (size_t i = 0; i < v->size; i++) {
        cfx_big_init(&new_data[i]);
        cfx_big_move(&new_data[i], &v->data[i]);   /* steals limbs */
    }

    
    for (size_t i = v->size; i < new_cap; i++) {
        cfx_big_init(&new_data[i]);
    }

    /* old buffer elements are now empty; safe to free backing array */
    free(v->data);
    v->data = new_data;
    v->cap  = new_cap;
    return 0;
}

int cfx_big_vec_resize_zero(cfx_big_vec_t* v, size_t new_size) {
    if (!v) return -1;
    if (new_size > v->cap) {
        int rc = cfx_big_vec_reserve(v, new_size);
        if (rc) return rc;
    }

    if (new_size < v->size) {
        for (size_t i = new_size; i < v->size; i++) {
            cfx_big_free(&v->data[i]);
            cfx_big_init(&v->data[i]); /* keep slot in known state */
        }
    } else if (new_size > v->size) {
        for (size_t i = v->size; i < new_size; i++) {
            /* slot already init'd in reserve()  */
            /* cfx_big_init(&v->data[i]); */
            cfx_big_from_limb(&v->data[i], 0);
        }
    }

    v->size = new_size;
    return 0;
}

static int cfx_big_vec_ensure_push(cfx_big_vec_t* v) {
    if (v->size < v->cap) return 0;
    size_t need = (v->cap == 0) ? CFX_BIG_VEC_MIN_CAP : (v->cap + 1);
    return cfx_big_vec_reserve(v, need);
}

int cfx_big_vec_push_u64(cfx_big_vec_t* v, cfx_limb_t value) {
    if (!v) return -1;
    int rc = cfx_big_vec_ensure_push(v);
    if (rc) return rc;

    /* element slot is already init'd */
    if (cfx_big_from_limb(&v->data[v->size], value) != 0) return -1;
    v->size++;
    return 0;
}

int cfx_big_vec_push_copy(cfx_big_vec_t* v, const cfx_big_t* x) {
    if (!v || !x) return -1;
    int rc = cfx_big_vec_ensure_push(v);
    if (rc) return rc;

    if (cfx_big_copy(&v->data[v->size], x) != 0) return -1;    /* deep copy */
    v->size++;
    return 0;
}

int cfx_big_vec_push_move(cfx_big_vec_t* v, cfx_big_t* x) {
    if (!v || !x) return -1;
    int rc = cfx_big_vec_ensure_push(v);
    if (rc) return rc;

    /* slot already init'd; move steals x */
    cfx_big_move(&v->data[v->size], x);
    v->size++;
    return 0;
}

int cfx_big_vec_pop_move(cfx_big_vec_t* v, cfx_big_t* out) {
    if (!v || !out) return -1;
    if (v->size == 0) return -1;

    size_t idx = v->size - 1;

    /* out must be init'd by caller */
    cfx_big_move(out, &v->data[idx]);
    /* keep slot in known state */
    cfx_big_init(&v->data[idx]);

    v->size--;
    return 0;
}

void cfx_big_vec_pop_free(cfx_big_vec_t* v) {
    if (!v || v->size == 0) return;
    size_t idx = v->size - 1;
    cfx_big_free(&v->data[idx]);
    cfx_big_init(&v->data[idx]);
    v->size--;
}
