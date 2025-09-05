#ifndef CFX_PRIMES_H
#define CFX_PRIMES_H

#include <stdint.h>
#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

/* find all primes up to n using sieve of eratosthenes */
cfx_vec_t cfx_sieve_primes(uint64_t n);

/* uses legendre's formula to give th power of p that divides n! */
uint64_t cfx_legendre(uint64_t n, uint64_t p);

#ifdef __cplusplus
}
#endif
#endif
