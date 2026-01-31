/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* ntt.h - Number-Theoretic Transform for FFT multiplication */

#ifndef CFX_NTT_H
#define CFX_NTT_H

#include <stdint.h>
#include <stddef.h>
#include "cfx/arith.h"

/*
 * NTT Primes
 *
 * For radix-2 NTT we need primes p where (p-1) has a large power-of-2 factor.
 * Primes of form p = k * 2^m + 1 allow transforms up to size 2^m.
 *
 * We use three primes for CRT reconstruction of large products.
 */

/* p1 = 2^64 - 2^32 + 1 = 0xFFFFFFFF00000001 (allows 2^32 transform size)
 * This is the "Goldilocks prime" - well-known and tested */
#define CFX_NTT_P1 UINT64_C(0xFFFFFFFF00000001)

/* p2 = 71 * 2^57 + 1 = 0x8E00000000000001 (allows 2^57 transform)
 * This prime has p-1 = 2^57 * 71 */
#define CFX_NTT_P2 UINT64_C(0x8E00000000000001)

/* p3 = 75 * 2^57 + 1 = 0x9600000000000001 (allows 2^57 transform)
 * This prime has p-1 = 2^57 * 75 = 2^57 * 3 * 5^2 */
#define CFX_NTT_P3 UINT64_C(0x9600000000000001)

/* primitive roots for each prime (g such that g^((p-1)/2) = p-1 mod p) */
#define CFX_NTT_G1 UINT64_C(7)   /* primitive root of p1 */
#define CFX_NTT_G2 UINT64_C(3)   /* primitive root of p2 */
#define CFX_NTT_G3 UINT64_C(7)   /* primitive root of p3 */

/*
 * Modular Arithmetic (mod 64-bit prime)
 *
 * These operate on values in [0, p-1] and return results in [0, p-1].
 */

/* a + b mod p */
uint64_t cfx_ntt_mod_add(uint64_t a, uint64_t b, uint64_t p);

/* a - b mod p */
uint64_t cfx_ntt_mod_sub(uint64_t a, uint64_t b, uint64_t p);

/* a * b mod p (uses Montgomery or Barrett reduction internally) */
uint64_t cfx_ntt_mod_mul(uint64_t a, uint64_t b, uint64_t p);

/* base^exp mod p (binary exponentiation) */
uint64_t cfx_ntt_mod_pow(uint64_t base, uint64_t exp, uint64_t p);

/* a^(-1) mod p (modular inverse via Fermat's little theorem: a^(p-2) mod p) */
uint64_t cfx_ntt_mod_inv(uint64_t a, uint64_t p);

/*
 * Primitive Root Utilities
 */

/* compute g^((p-1)/n) mod p - the n-th root of unity */
uint64_t cfx_ntt_root_of_unity(uint64_t g, uint64_t p, size_t n);

/*
 * Twiddle Factor Table
 *
 * Precomputed powers of the primitive root for a given transform size.
 */
typedef struct {
    uint64_t *forward;   /* twiddles for forward NTT */
    uint64_t *inverse;   /* twiddles for inverse NTT */
    size_t n;            /* transform size (power of 2) */
    uint64_t p;          /* the prime modulus */
    uint64_t n_inv;      /* n^(-1) mod p for inverse scaling */
} cfx_ntt_twiddles_t;

/* allocate and precompute twiddle factors for size n with prime p and root g */
int cfx_ntt_twiddles_init(cfx_ntt_twiddles_t *tw, size_t n, uint64_t p, uint64_t g);

/* free twiddle table */
void cfx_ntt_twiddles_free(cfx_ntt_twiddles_t *tw);

/*
 * Bit-Reversal Permutation
 */

/* in-place bit-reversal permutation of array a of length n (n must be power of 2) */
void cfx_ntt_bit_reverse(uint64_t *a, size_t n);

/*
 * Core NTT Transforms
 */

/* forward NTT: a -> DFT(a) mod p */
void cfx_ntt_forward(uint64_t *a, size_t n, uint64_t p, const uint64_t *twiddles);

/* inverse NTT: a -> IDFT(a) mod p (includes 1/n scaling) */
void cfx_ntt_inverse(uint64_t *a, size_t n, uint64_t p,
    const uint64_t *twiddles, uint64_t n_inv);

/*
 * Polynomial Convolution via NTT
 *
 * Computes c = a * b (polynomial multiplication) mod p using NTT.
 * Result c has length n (must be >= len_a + len_b - 1).
 *
 * Parameters:
 *   c     - output array of length n (will be overwritten)
 *   a     - first input polynomial coefficients
 *   len_a - length of a
 *   b     - second input polynomial coefficients
 *   len_b - length of b
 *   n     - transform size (must be power of 2, >= len_a + len_b - 1)
 *   p     - prime modulus
 *   g     - primitive root of p
 *
 * Returns 0 on success, -1 on error (allocation failure or invalid n).
 */
int cfx_ntt_convolve(uint64_t *c, const uint64_t *a, size_t len_a,
    const uint64_t *b, size_t len_b,
    size_t n, uint64_t p, uint64_t g);

/*
 * Big Integer Multiplication via NTT
 */

/* chunk size for splitting limbs - 16 bits provides ample headroom */
#define CFX_NTT_CHUNK_BITS 16
#define CFX_NTT_CHUNK_MASK 0xFFFFu

/*
 * Convert limb array to NTT chunks (16-bit pieces).
 * Returns number of chunks written.
 *
 * Each limb is split into CFX_LIMB_BITS/16 chunks (2 for 32-bit, 4 for 64-bit).
 * Output chunks are stored little-endian (low chunks first).
 */
size_t cfx_ntt_limbs_to_chunks(uint64_t *chunks, size_t max_chunks,
    const cfx_limb_t *limbs, size_t n_limbs);

/*
 * Convert NTT output chunks back to limbs with carry propagation.
 * Returns number of limbs written (excluding trailing zeros).
 *
 * Input chunks may have values > 16 bits (from convolution sums).
 * Carry propagation reduces each to 16 bits and accumulates into limbs.
 */
size_t cfx_ntt_chunks_to_limbs(cfx_limb_t *limbs, size_t max_limbs,
    const uint64_t *chunks, size_t n_chunks);

/*
 * Multiply two limb arrays using NTT.
 * out must have space for (n_a + n_b) limbs.
 * Returns number of result limbs (excluding trailing zeros).
 */
size_t cfx_ntt_mul_limbs(cfx_limb_t *out, size_t out_cap,
    const cfx_limb_t *a, size_t n_a,
    const cfx_limb_t *b, size_t n_b);

#endif /* CFX_NTT_H */
