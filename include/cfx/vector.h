/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_VECTOR_H
#define CFX_VECTOR_H

#include "cfx/num.h"

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_VEC_MIN_CAP 16

typedef struct {
    cfx_limb_t* data;
    size_t size;
    size_t cap;
} cfx_vec_t;

void cfx_vec_init(cfx_vec_t *v);

int  cfx_vec_reserve(cfx_vec_t *v, size_t need);
int  cfx_vec_resize(cfx_vec_t *v, size_t sz, cfx_limb_t fill);
int  cfx_vec_push(cfx_vec_t *v, cfx_limb_t value);

void cfx_vec_clear(cfx_vec_t *v);
void cfx_vec_free(cfx_vec_t *v);

void cfx_vec_print(const cfx_vec_t *v);

#ifdef __cplusplus
}
#endif
#endif
