#ifndef CFX_BIG_H
#define CFX_BIG_H

#include "cfx/fac.h"
#include "cfx/types.h"

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* minimal big int struct */
typedef struct {
    uint64_t* limb; /* the 'digits' of the number with base BIG_BASE*/
    size_t n;
    size_t cap;
} cfx_big_t;



void cfx_big_init(cfx_big_t* b);
void cfx_big_free(cfx_big_t* b);
void cfx_big_reserve(cfx_big_t* b, size_t need);
void cfx_big_set_val(cfx_big_t* b, uint64_t v);
void cfx_big_add_sm(cfx_big_t* b, uint64_t n);
void cfx_big_sub_sm(cfx_big_t* b, uint64_t n);
void cfx_big_mul_sm(cfx_big_t* b, uint64_t m);

/* Multiply by p^e by repeated squaring using small chunks to avoid overflow */
void cfx_big_powmul_prime(cfx_big_t* b, uint64_t p, uint64_t e);

/* Materialize factorization into cfx_big_t */
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t* f);
void cfx_big_to_fac(cfx_fac_t* f, const cfx_big_t* b);

char* cfx_big_to_str(const cfx_big_t* b, size_t *sz_out);
int cfx_big_from_str(cfx_big_t* b, const char* str);

/**
 *  args: cfx_big_t * b : ptr to big
 *  fmt: format string
 * */
#define CFX_BIG_PRINTF(b, fmt, ...) \
do {\
    char *s = cfx_big_to_str(b, NULL); \
    CFX_PRINT_DBG(fmt "%s\n", ##__VA_ARGS__, s); \
    free(s); \
} while (0)

#ifdef __cplusplus
}
#endif
#endif
