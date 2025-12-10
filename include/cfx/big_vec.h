/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
#ifndef CFX_BIG_VEC_H
#define CFX_BIG_VEC_H

#include "cfx/big.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_BIG_VEC_MIN_CAP 16

typedef struct {
    cfx_big_t* data;
    size_t     size;
    size_t     cap;
} cfx_big_vec_t;

void cfx_big_vec_init(cfx_big_vec_t* v);
void cfx_big_vec_clear(cfx_big_vec_t* v); /* frees elements, keeps buffer */
void cfx_big_vec_free(cfx_big_vec_t* v);  /* clear + free buffer */

int  cfx_big_vec_reserve(cfx_big_vec_t* v, size_t need);
int  cfx_big_vec_resize_zero(cfx_big_vec_t* v, size_t new_size); /* new elems init to 0 */

int  cfx_big_vec_push_u64 (cfx_big_vec_t* v, cfx_limb_t value);
int  cfx_big_vec_push_copy(cfx_big_vec_t* v, const cfx_big_t* x);
int  cfx_big_vec_push_move(cfx_big_vec_t* v, cfx_big_t* x); /* steals x */

int  cfx_big_vec_pop_move (cfx_big_vec_t* v, cfx_big_t* out); /* moves last into out */
void cfx_big_vec_pop_free (cfx_big_vec_t* v);                 /* frees last */

#ifdef __cplusplus
}
#endif

#endif /* CFX_BIG_VEC_H */
