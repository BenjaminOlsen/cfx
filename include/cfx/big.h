#ifndef CFX_BIG_H
#define CFX_BIG_H

#include <stdint.h>
#include <stddef.h>

#include "factorization.h"

#ifdef __cplusplus
extern "C" {
#endif

/* minimal big int struct */
typedef struct {
    uint32_t* limb; /* the 'digits' of the number with base BIG_BASE*/
    size_t n;
    size_t cap;
} cfx_big_t;

/* 1e9 */
#define CFX_BIG_BASE 1000000000u



void cfx_big_init(cfx_big_t* b);
void cfx_big_free(cfx_big_t* b);
void cfx_big_reserve(cfx_big_t* b, size_t need);
void cfx_big_set_u32(cfx_big_t* b, uint32_t v);

void cfx_big_mul_u32(cfx_big_t* b, uint32_t m);
/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
void cfx_big_powmul_prime(cfx_big_t *b, uint32_t p, uint32_t e);
/* Materialize factorization into cfx_big_t */
void cfx_big_from_factorization(cfx_big_t *b, const cfx_factorization_t *f);
char* cfx_big_to_string(const cfx_big_t *b);

#ifdef __cplusplus
}
#endif
#endif