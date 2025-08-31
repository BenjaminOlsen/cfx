#ifndef CFX_VECTOR_H
#define CFX_VECTOR_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    uint32_t* data;
    size_t size;
} cfx_vec_u32;

void cfx_vec_u32_free(cfx_vec_u32* v);
void cfx_vec_u32_print(cfx_vec_u32* v);

#ifdef __cplusplus
}
#endif
#endif