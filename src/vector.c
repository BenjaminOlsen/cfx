#include "cfx/vector.h"
#include <stdio.h>
#include <stdlib.h>

void cfx_vec_u32_free(cfx_vec_u32* v) {
    free(v->data);
    v->size = 0;
}

void cfx_vec_u32_print(cfx_vec_u32* v) {
    printf("v: %p, size: %zu\n", v, v->size);
    for (size_t i = 0; i < v->size; ++i) {
        printf("%u ", v->data[i]);
    }
    printf("\n");
}