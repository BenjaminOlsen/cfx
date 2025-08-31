#ifndef CFX_PRIMES_H
#define CFX_PRIMES_H

#include <stdint.h>
#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

/* find all primes up to n using sieve of eratosthenes */
cfx_vec_u32 cfx_sieve_primes(uint32_t n);

/* uses legendre's formula to give th power of p that divides n! */
uint32_t cfx_legendre(uint32_t n, uint32_t p);

#ifdef __cplusplus
}
#endif
#endif