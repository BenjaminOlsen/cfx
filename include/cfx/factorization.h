#ifndef CFX_FACTORIZATION_H
#define CFX_FACTORIZATION_H

#include <stdint.h>
#include <stddef.h>
#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

 /* a single prime factor: */
typedef struct {
    uint32_t p;  // prime
    uint32_t e;  // exponent
} cfx_pf_t;

/**
 * cfx_factorization_t: a struct to hold a prime factorization of a number
 */
typedef struct {
    cfx_pf_t* data;
    size_t len;
    size_t cap;
} cfx_factorization_t;



void cfx_fac_print(cfx_factorization_t* f);
void cfx_fac_init(cfx_factorization_t* f);
void cfx_fac_free(cfx_factorization_t* f) ;
void cfx_fac_reserve(cfx_factorization_t* f, size_t req_cap);
void cfx_fac_push(cfx_factorization_t* f, uint32_t p, uint32_t e);
/* Deep copy: dst becomes a copy of src, using cfx_fac_push for each element. */
void cfx_fac_copy(cfx_factorization_t *dst, const cfx_factorization_t *src);
// dst += src
void cfx_fac_add(cfx_factorization_t* dst, cfx_factorization_t* src);
/* dst -= src */
void cfx_fac_sub(cfx_factorization_t* dst, cfx_factorization_t* src) ;
/* calculate the factorization of n.
we pass in a list of primes to not have to calculate it on every call of this function */
cfx_factorization_t cfx_fac_factorial(uint32_t n, const cfx_vec_t *primes);
/* Factorization of C(n,k) = n! / (k! (n-k)!) */
cfx_factorization_t cfx_fac_binom(uint32_t n, uint32_t k);

#ifdef __cplusplus
}
#endif
#endif 