/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_BIG_H
#define CFX_BIG_H

#include "cfx/fac.h"
#include "cfx/algo.h"
#include "cfx/arith.h"

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_BIG_PRINT(bg) printf("b.n %zu: b.cap: %zu, b %s\n", bg.n, bg.cap, cfx_big_dec_alloc(&bg, NULL))

#define CFX_BIG_PRINT_LIMBS(bg) \
        do { \
            for (size_t i = 0; i < bg.n; ++i) { \
                printf("limb[%zu] = " CFX_PRIuLIMB "\n", i, bg.limb[i]); \
            } \
        } while (0)

/** cfx_big_t - arbitrary precision positive integer type, represented in base 2^64 limbs
 * limb: ptr to array of limbs
 * n: current number of limbs
 * cap: capacity of underlying limb array
 *
 * >>>>> Always call cfx_big_init after declaring a cfx_big_t !
 **/
struct cfx_big {
    cfx_limb_t *limb; /* the 'digits' of the number with base BIG_BASE */
    size_t n;
    size_t cap;
};
typedef struct cfx_big cfx_big_t;

void cfx_big_init(cfx_big_t *b);
void cfx_big_clear(cfx_big_t *b);
void cfx_big_free(cfx_big_t *b);
void cfx_big_reserve(cfx_big_t *b, size_t need);
void cfx_big_assign(cfx_big_t *dst, const cfx_big_t *src);
void cfx_big_assign_sm(cfx_big_t *dst, const cfx_limb_t src);
void cfx_big_assign_zero(cfx_big_t *b);
void cfx_big_assign_one(cfx_big_t *b);
int cfx_big_endswith_u64(const cfx_big_t *x, uint64_t value);
int cfx_big_copy(cfx_big_t *dst, const cfx_big_t *src);
void cfx_big_move(cfx_big_t *dst, cfx_big_t *src);
int cfx_big_is_zero(const cfx_big_t *b);

int cfx_big_is_one(const cfx_big_t *b);
int cfx_big_is_even(const cfx_big_t *b);
int cfx_big_eq_u64(const cfx_big_t *b, cfx_limb_t n);
int cfx_big_eq(const cfx_big_t *b1, const cfx_big_t *b2);
int cfx_big_cmp(const cfx_big_t *a, const cfx_big_t *b);
int cfx_big_cmp_sm(const cfx_big_t *a, cfx_limb_t n);
void cfx_big_swap(cfx_big_t *a, cfx_big_t *b);
void cfx_big_cswap(cfx_big_t *a, cfx_big_t *b, int condition); /* (maybe) constant-time conditional swap */

/* Bitwise operations: out = a OP b */
void cfx_big_and(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);
void cfx_big_or(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);
void cfx_big_xor(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);

/* In-place bitwise: a OP= b */
void cfx_big_and_eq(cfx_big_t *a, const cfx_big_t *b);
void cfx_big_or_eq(cfx_big_t *a, const cfx_big_t *b);
void cfx_big_xor_eq(cfx_big_t *a, const cfx_big_t *b);

void cfx_big_rotl(cfx_big_t *out, const cfx_big_t *b, unsigned r);
void cfx_big_rotr(cfx_big_t *out, const cfx_big_t *b, unsigned r);
void cfx_big_rotl_w(cfx_big_t *out, const cfx_big_t *b, unsigned r, unsigned w);
void cfx_big_rotr_w(cfx_big_t *out, const cfx_big_t *b, unsigned r, unsigned w);

void cfx_big_mask_bits(cfx_big_t *a, unsigned nbits);

/* Single-bit operations */
int cfx_big_bit_is_set(const cfx_big_t *x, size_t bit);
void cfx_big_bit_set(cfx_big_t *x, size_t bit);
void cfx_big_bit_clear(cfx_big_t *x, size_t bit);
void cfx_big_bit_flip(cfx_big_t *x, size_t bit);
size_t cfx_big_popcount(const cfx_big_t *x);

size_t cfx_big_bitlen(const cfx_big_t *b);      /* assumes b->limb[b->n - 1] != 0 */
int cfx_big_from_u64(cfx_big_t *b, uint64_t v);
int cfx_big_from_limb(cfx_big_t *b, cfx_limb_t v);
int cfx_big_to_bytes_be(uint8_t *out, size_t *out_len, const cfx_big_t *b);
int cfx_big_from_bytes_be(cfx_big_t *out, const uint8_t *be, size_t len);
void cfx_big_mul(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);

void cfx_big_mul_eq(cfx_big_t *b, const cfx_big_t *m);
void cfx_big_mul_eq_fft(cfx_big_t *b, const cfx_big_t *m);
void cfx_big_mul_eq_csa(cfx_big_t *b, const cfx_big_t *m);
void cfx_big_mul_fft(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);
/* assumes scratch is allocated with the appropriate size b->n + m->n already. */
void cfx_big_mul_csa_scratch(cfx_big_t *b, const cfx_big_t *m, cfx_mul_scratch_t *scratch);
/* if CFX_HAS_PTHREAD */
void cfx_big_mul_rows_pthreads(cfx_big_t *b, const cfx_big_t *m, int threads);

/* chooses the fastest in-place multiplication for the size of the multiplicands */
void cfx_big_mul_auto(cfx_big_t *b, const cfx_big_t *m);

void cfx_big_sq_eq(cfx_big_t *b);

/* Non-in-place arithmetic */
void cfx_big_add(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);
void cfx_big_sub(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b); /* assumes a >= b */

/* In-place arithmetic */
void cfx_big_add_eq(cfx_big_t *b, const cfx_big_t *a);
void cfx_big_add_sm_eq(cfx_big_t *b, cfx_limb_t n);
void cfx_big_sub_eq(cfx_big_t *a, const cfx_big_t *b);
void cfx_big_sub_sm_eq(cfx_big_t *b, cfx_limb_t n);
void cfx_big_mul_sm_eq(cfx_big_t *b, cfx_limb_t m);

/* Knuth's long division from TAOCP Vol. 2, 4.3.1 */
int cfx_big_divrem(cfx_big_t *q, cfx_big_t *r, const cfx_big_t *u, const cfx_big_t *v);

/* q = floor(a/b) */
int cfx_big_div(cfx_big_t *q, const cfx_big_t *a, const cfx_big_t *b);

/* In-place: a /= b; optional remainder r. */
int cfx_big_divrem_eq(cfx_big_t *a, const cfx_big_t *b, cfx_big_t *r /*nullable*/);

cfx_limb_t cfx_big_div_sm_eq(cfx_big_t *b, cfx_limb_t d);
uint32_t cfx_big_div_sm_u32_eq(cfx_big_t *b, uint32_t d);

/* out = n % m */
int cfx_big_mod(cfx_big_t *out, const cfx_big_t *n, const cfx_big_t *m);

/* returns b % m */
cfx_limb_t cfx_big_mod_sm(const cfx_big_t *b, cfx_limb_t m);

/* TODO: factor big into fac */
int cfx_fac_from_big(cfx_fac_t *fac, const cfx_big_t *in);

/* Bitshift operations */
/* b >>= s */
void cfx_big_shr_bits_eq(cfx_big_t *b, unsigned s);

/* b <<= s */
void cfx_big_shl_bits_eq(cfx_big_t *b, unsigned s);

/* out = x << s (0..63)  */
void cfx_big_shl_bits(cfx_big_t *out, const cfx_big_t *x, unsigned s);

/* out = x >> s (0..63)  */
void cfx_big_shr_bits(cfx_big_t *out, const cfx_big_t *x, unsigned s);

/* Multiply by p^e by repeated squaring using small chunks to avoid overflow */
void cfx_big_expmul_prime(cfx_big_t *b, cfx_limb_t p, cfx_limb_t e);

/* out = p^e */
void cfx_big_pow_sm(cfx_big_t *out, cfx_limb_t p, cfx_limb_t e);

/* out = n ^ p */
void cfx_big_exp(cfx_big_t *out, const cfx_big_t *n, const cfx_big_t *p);
void cfx_big_exp_u64(cfx_big_t *out, const cfx_big_t *n, cfx_limb_t p);

/* out = (n^p) mod m */
void cfx_big_mod_exp(cfx_big_t *out, const cfx_big_t *n, const cfx_big_t *p, const cfx_big_t *m);

/* out = (a*b) mod m*/
int cfx_big_mulmod(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b, const cfx_big_t *m);

int cfx_big_is_prime(const cfx_big_t *b);

/* Binary GCD algorithm for big integers */
void cfx_big_gcd(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b);

/* Extended GCD for big integers.
 * Computes g = gcd(a, b) and signed coefficients x, y such that ax + by = g.
 * x and y are cfx_sbig_t* (signed big integers) since they can be negative.
 * Any of g, x, y can be NULL if that output is not needed.
 * Forward declaration - cfx_sbig_t is defined in cfx/sbig.h */
typedef struct cfx_sbig cfx_sbig_t;
void cfx_big_xgcd(cfx_big_t *g, cfx_sbig_t *x, cfx_sbig_t *y,
    const cfx_big_t *a, const cfx_big_t *b);

/* Pollard-Rho factorization using Montgomery multiplication.
 * Returns a non-trivial factor, or a copy of n if n is prime/unfactorable. */
void cfx_big_pollard_rho(cfx_big_t *factor, const cfx_big_t *n);

void cfx_big_from_limbs(cfx_big_t *b, const cfx_limb_t *limbs, size_t n);
void cfx_big_from_fac(cfx_big_t *b, const cfx_fac_t *f);
void cfx_big_from_fac_fast(cfx_big_t *out, const cfx_fac_t *f);
void cfx_big_from_fac_faster(cfx_big_t *out, const cfx_fac_t *f);
/*
 * Factorize b into prime factors, storing results in f.
 * Small primes (≤64 bits) go in f->data, large primes (>64 bits) go in f->big_primes.
 * If remainder is non-NULL and factorization is incomplete, the unfactored
 * composite is written there.
 *
 * Returns:
 *   0  - complete factorization
 *   1  - incomplete (remainder holds unfactored composite)
 *  -1  - error
 */
int cfx_big_to_fac(cfx_fac_t *f, const cfx_big_t *b, cfx_big_t *remainder);

/* allocating string conversions - caller must free returned pointer */
char * cfx_big_dec_alloc(const cfx_big_t *b, size_t *len_out);
char * cfx_big_hex_alloc(const cfx_big_t *b, size_t *len_out);
char * cfx_big_bin_alloc(const cfx_big_t *b, size_t *len_out);
char * cfx_big_b64_alloc(const cfx_big_t *b, size_t *len_out);
int cfx_big_to_sci(const cfx_big_t *x, unsigned base, int sig_digits, char *out, size_t outsz);

/* snprintf-style: write to out[0..outlen-1], always NUL-terminate,
   return number of chars that would be written (excluding NUL). if out==NULL, just return size. */
size_t cfx_big_snprint_dec(const cfx_big_t *b, char *out, size_t outlen);
size_t cfx_big_snprint_hex(const cfx_big_t *b, char *out, size_t outlen);
size_t cfx_big_snprint_bin(const cfx_big_t *b, char *out, size_t outlen);

int cfx_big_from_str(cfx_big_t *out, const char *s);

int cfx_big_from_bin(cfx_big_t *out, const char *s);
int cfx_big_from_dec(cfx_big_t *out, const char *str);
int cfx_big_from_hex(cfx_big_t *out, const char *s);
int cfx_big_from_oct(cfx_big_t *out, const char *s);
int cfx_big_from_file(cfx_big_t *out, FILE *fp, int base);

int cfx_big_scan_num_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed);
int cfx_big_scan_hex_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed);
int cfx_big_scan_dec_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed);
int cfx_big_scan_bin_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed);
int cfx_big_scan_b64_n(cfx_big_t *out, const uint8_t *s, size_t len, size_t *consumed);
int cfx_big_from_oct_n(cfx_big_t *out, const uint8_t *in, size_t in_len, size_t *consumed);


/* ================== Montgomery ================== */
/**
 * Montgomery arithmetic context for a fixed odd modulus n.
 *
 * Let b = 2^64 and k = number of 64-bit limbs of n. Define R = b^k.
 * This context precomputes the constants needed for CIOS Montgomery
 * reduction and the Montgomory “to/from” transforms:
 *
 *   n      : modulus (normalized, odd)
 *   k      : limb count of n; also defines R = b^k
 *   n0inv  : (-n[0])^{-1} mod b, i.e. n[0] * n0inv ≡ -1 (mod b)
 *   rr     : R^2 mod n  (used to convert x -> xR mod n via MontMul)
 *   Additional caches:
 *   R1  : R mod n    (the Montgomery representation of 1)
 *   mu  : floor(b^(2k)/n) for Barrett reduction of large inputs
 *
 * With this context:
 *   MontMul(a,b) = a * b * R^{-1} mod n   (for 0 ≤ a,b < n)
 *   MontTo(x)    = MontMul(x, rr)  = xR mod n
 *   MontFrom(x)  = MontMul(x, 1)   = xR^{-1} mod n
 *
 * Requirements/notes:
 * - n must be odd.
 * - Operands are assumed reduced (0 ≤ x < n) unless you explicitly
 *   reduce them before MontTo.
 * - The context is immutable after init and may be shared across threads.
 */

typedef struct {
    cfx_big_t n;        /* modulus (odd), normalized */
    size_t k;           /* limb count of n (R = 2^(64*k)) */
    cfx_limb_t n0inv;     /* -n^{-1} mod 2^64 (requires n.limb[0] odd) */
    cfx_big_t rr;       /* R^2 mod n */
    cfx_big_t R1;       /* R mod n -- (Montgomery "1") */
} cfx_big_mont_ctx_t;

int cfx_big_mont_ctx_init(cfx_big_mont_ctx_t *ctx, const cfx_big_t *n);
void cfx_big_mont_ctx_free(cfx_big_mont_ctx_t *ctx);

/* Conversions between normal domain and Montgomery domain */
int cfx_big_mont_to   (cfx_big_t *out, const cfx_big_t *a,  const cfx_big_mont_ctx_t *ctx); /* out = a*R mod n */
int cfx_big_mont_from (cfx_big_t *out, const cfx_big_t *aR, const cfx_big_mont_ctx_t *ctx); /* out = aR*R^-1 */

/* Core Montgomery ops (operands/results in Montgomery domain) */
/**
 * Montgomery product (CIOS): out = a * b * R^{-1} mod n
 *
 * Base b = 2^64, R = b^k where k = ctx->k (number of 64-bit limbs of n).
 * Uses the Coarsely Integrated Operand Scanning (CIOS) algorithm:
 *   for i = 0..k-1:
 *     T += a * b[i]
 *     m = (T[0] * n0inv) mod b         n0inv ≡ -n[0]^{-1} (mod b)
 *     T += m * n
 *     T >>= 64                         drop least-significant limb
 * Final step does a conditional subtraction so that 0 ≤ out < n.
 *
 * Preconditions:
 *   - ctx is initialized for an odd modulus n (ctx->n), with ctx->k, ctx->n0inv, ctx->rr set.
 *   - Inputs are reduced: 0 ≤ a,b < n (i.e., at most k limbs). Higher limbs are NOT consumed.
 *
 * Postconditions:
 *   - out holds (a·b·R^{-1}) mod n in canonical form (no leading zero limbs), 0 ≤ out < n.
 *   - alias-safe: out may alias a or b.
 *
 * Complexity / resources:
 *   - Time: O(k^2) word multiplies/adds.
 *   - Space: O(k) temporaries (internal accumulator T of k+2 limbs).
 *
 * Side-channel note:
 *   - Loop bounds and memory access are fixed; the final conditional subtraction is value-dependent.
 *     Replace with a masked subtraction if strict constant-time behavior is required.
 *
 * Returns 1 on success, 0 on invalid arguments or internal errors.
 */

int cfx_big_mont_mul(cfx_big_t *out, const cfx_big_t *aR, const cfx_big_t *bR, const cfx_big_mont_ctx_t *ctx);
int cfx_big_modexp_binary(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *e, const cfx_big_mont_ctx_t *ctx);

static inline int cfx_big_mont_sqr(cfx_big_t *out, const cfx_big_t *aR, const cfx_big_mont_ctx_t *ctx) {
    return cfx_big_mont_mul(out, aR, aR, ctx);
}


/* one liners that use montgomery internally */
int cfx_big_mul_mod(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b, const cfx_big_t *n);
int cfx_big_sqr_mod(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *n);
int cfx_big_modexp(cfx_big_t *out, const cfx_big_t *base, const cfx_big_t *exp, const cfx_big_t *n);

double cfx_big_log(const cfx_big_t *b, double base);

/**
 *  args: cfx_big_t * b : ptr to big
 *  fmt: format string
 * */
#define CFX_BIG_PRINTF(b, ...) \
        do { \
            char *_cfx_big_printf_var = cfx_big_dec_alloc(b, NULL); \
            CFX_PRINT_DBG("%s\n", __VA_ARGS__, _cfx_big_printf_var); \
            free(_cfx_big_printf_var); \
        } while (0)

#ifdef __cplusplus
}
#endif
#endif
