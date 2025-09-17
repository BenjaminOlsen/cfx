#ifndef CFX_ALGO_H
#define CFX_ALGO_H

#include <stdint.h>
#include "vector.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Accumulator holds a value as:  (acc_hi << 64) + acc_lo
   No cross-limb carry is propagated while accumulating. */
typedef struct {
    uint64_t lo;  // sum bits
    uint64_t hi;  // saved carries (conceptually << 1 bit each)
} csa128_t;

typedef struct {
    csa128_t* acc;     // nout entries
    uint64_t* spill;   // nout entries (small counters)
    size_t    cap;
} cfx_mul_scratch_t;

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

/* A * B -> R  (schoolbook) with carry-save accumulation.
   - A has na limbs, B has nb limbs
   - R must have space for na+nb limbs
*/
void cfx_mul_csa_portable(const uint64_t* A, size_t na,
                          const uint64_t* B, size_t nb,
                          uint64_t* R);

void cfx_mul_csa_portable_fast(const uint64_t* A, size_t na,
                               const uint64_t* B, size_t nb,
                               uint64_t* R,
                               cfx_mul_scratch_t* scratch);

void cfx_mul_scratch_alloc(cfx_mul_scratch_t* s, size_t needed);
void cfx_mul_scratch_free(cfx_mul_scratch_t* s);

// Zero only the first `nout` entries that the multiplication will touch
void cfx_mul_scratch_zero(cfx_mul_scratch_t* s, size_t nout);

/* portable CSA inner: accumulate rows i in [ia, ia+na_rows) into acc/spill */
void cfx_mul_csa_portable_fast_rows(const uint64_t* A, size_t ia, size_t na_rows,
                                    const uint64_t* B, size_t nb,
                                    csa128_t* acc, uint64_t* spill);

/* fold spills once and normalize: writes R[0..nout-1] */
void cfx_mul_csa_fold_and_normalize(csa128_t* acc, uint64_t* spill,
                                    size_t nout, uint64_t* R);

#ifdef __cplusplus
}
#endif
#endif
