#ifndef CFX_ALGO_H
#define CFX_ALGO_H

#include <stdint.h>
#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

/* find all primes up to n using sieve of eratosthenes */
cfx_vec_t cfx_sieve_primes(uint64_t n);

/* uses legendre's formula to give the power of p that divides n! */
uint64_t cfx_legendre(uint64_t n, uint64_t p);

int cfx_is_prime_u64(uint64_t n);

/* modular multiplication (a*b) mod m */
uint64_t cfx_mulmod_u64(uint64_t a, uint64_t b, uint64_t m);

/* modular exponentiation (a^e) mod m */
uint64_t cfx_powmod_u64(uint64_t a, uint64_t e, uint64_t m);

uint64_t cfx_gcd_u64(uint64_t a, uint64_t b);
uint64_t cfx_rho_brent(uint64_t n);
int cfx_factor_u64(cfx_vec_t* primes, cfx_vec_t* exps, uint64_t n);

#ifdef __cplusplus
}
#endif
#endif
