#ifndef CFX_VECTOR_H
#define CFX_VECTOR_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_VEC_MIN_CAP 16

typedef struct {
    uint64_t* data;
    size_t size;
    size_t cap;
} cfx_vec_t;

void cfx_vec_init(cfx_vec_t *v);

int  cfx_vec_reserve(cfx_vec_t *v, size_t need);
int  cfx_vec_resize(cfx_vec_t *v, size_t sz, uint64_t fill);
int  cfx_vec_push(cfx_vec_t *v, uint64_t value);

void cfx_vec_clear(cfx_vec_t *v);
void cfx_vec_free(cfx_vec_t *v);

void cfx_vec_print(const cfx_vec_t *v);

#ifdef __cplusplus
}
#endif
#endif
