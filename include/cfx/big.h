#ifndef CFX_BIG_H
#define CFX_BIG_H

#include "cfx/fac.h"
#include "cfx/fmt.h"
#include "cfx/algo.h"
#include "cfx/types.h"

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_BIG_PRINT(bg) printf("b.n %zu: b.cap: %zu, b %s\n", bg.n, bg.cap, cfx_big_to_str(&bg, NULL))

#define CFX_BIG_PRINT_LIMBS(bg) \
    do { \
        for (size_t i = 0; i < bg.n; ++i) { \
            printf("limb[%zu] = "U64F"\n", i, bg.limb[i]); \
        } \
    } while (0)

typedef enum {
    CFX_FAC_NONE,
    CFX_FAC_PARTIAL,
    CFX_FAC_FULL
} cfx_fac_state_e;

typedef struct cfx_fac_cache cfx_fac_cache_t;

/** cfx_big_t - arbitrary precision positive integer type, represented in base 2^64 limbs
 * limb: ptr to array of limbs
 * n: current number of limbs
 * cap: capacity of underlying limb array
 * 
 * >>>>> Always call cfx_big_init after declaring a cfx_big_t ! 
 **/
typedef struct {
    cfx_u64_t* limb; /* the 'digits' of the number with base BIG_BASE*/
    size_t n;
    size_t cap;
    cfx_fac_cache_t* cache;  /* == NULL if prime factor caching disabled */
} cfx_big_t;

/**  cache for prime factors of a cfx_big_t.
 * 
 * Invariants: (assume value N of cfx_big_t:)
 *  - N == (∏ p_i^e_i) * cofactor.
 *  - cofactor ≥ 1.
 *  - cofactor has no prime factor ≤ bnd.
 *  - state == CFX_FAC_FULL iff cofactor == 1.
 * 
 * for N = 1: primes empty, cofactor = 1, state = CFX_FAC_FULL.
*/
struct cfx_fac_cache {
    cfx_fac_t primes;      /* your (p, exp) list; can be sorted & coalesced */
    cfx_big_t cofactor;    /* remaining unfactored part (>= 1) */
    cfx_fac_state_e state;
    cfx_u64_t bnd;           /* proven “B-smoothness” bound:
                             * "all primes <= bnd have been factored out of cofactor already"
                             */
};

void cfx_big_init(cfx_big_t* b);
void cfx_big_clear(cfx_big_t* b);
void cfx_big_free(cfx_big_t* b);
int cfx_big_reserve(cfx_big_t* b, size_t need);
void cfx_big_assign(cfx_big_t* dst, const cfx_big_t* src);
int cfx_big_copy(cfx_big_t* dst, const cfx_big_t* src);
void cfx_big_move(cfx_big_t* dst, cfx_big_t* src);
int cfx_big_is_zero(const cfx_big_t* b);
int cfx_big_eq_u64(const cfx_big_t* b, cfx_u64_t n);
int cfx_big_eq(const cfx_big_t* b1, const cfx_big_t* b2);
int cfx_big_cmp(const cfx_big_t* a, const cfx_big_t* b);
void cfx_big_swap(cfx_big_t* a, cfx_big_t* b);

int cfx_big_from_u64(cfx_big_t* b, cfx_u64_t v);
void cfx_big_mul(cfx_big_t* b, const cfx_big_t* m);
void cfx_big_mul_fft(cfx_big_t* b, const cfx_big_t* m); /* todo */
void cfx_big_mul_csa(cfx_big_t* b, const cfx_big_t* m);
/* assumes scratch is allocated with the appropriate size b->n + m->n already. */
void cfx_big_mul_csa_scratch(cfx_big_t* b, const cfx_big_t* m, cfx_mul_scratch_t* scratch);
void cfx_big_mul_rows_pthreads(cfx_big_t* b, const cfx_big_t* m, int threads);
/* chooses the fastest multiplication for the size of the multiplicands */
void cfx_big_mul_auto(cfx_big_t* b, const cfx_big_t* m);

void cfx_big_sq(cfx_big_t* b);
int cfx_big_div(cfx_big_t* b, const cfx_big_t* d, cfx_big_t* r); /* b /= d; r remainder. */
void cfx_big_add(cfx_big_t* b, const cfx_big_t* a);
void cfx_big_add_sm(cfx_big_t* b, cfx_u64_t n);
void cfx_big_sub(cfx_big_t* a, const cfx_big_t* b);
void cfx_big_sub_sm(cfx_big_t* b, cfx_u64_t n);
void cfx_big_mul_sm(cfx_big_t* b, cfx_u64_t m);

/* Knuth's long division from TAOCP Vol. 2, 4.3.1 */
int cfx_big_divrem(cfx_big_t* q, cfx_big_t* r, const cfx_big_t* u, const cfx_big_t* v);

/* q = floor(n/d)*/
int cfx_big_div_out(cfx_big_t* q, const cfx_big_t* n, const cfx_big_t* d);

/* In-place: b := floor(b/d); optional remainder r. Alias-safe for any combination. */
int cfx_big_div_eq(cfx_big_t* b, const cfx_big_t* d, cfx_big_t* r /*nullable*/);

cfx_u64_t cfx_big_div_sm(cfx_big_t* b, cfx_u64_t d);
uint32_t cfx_big_div_sm_u32(cfx_big_t* b, uint32_t d);

/* out = n % m */
int cfx_big_mod(cfx_big_t* out, const cfx_big_t* n, const cfx_big_t* m);

/* returns b % m */
cfx_u64_t cfx_big_mod_sm(cfx_big_t* b, cfx_u64_t m);


/* Bitshift operations */
/* b >>= s */
void cfx_big_shr_bits_eq(cfx_big_t* b, unsigned s);
/* b <<= s */
void cfx_big_shl_bits_eq(cfx_big_t* b, unsigned s);
/* out = x << s (0..63)  */
void cfx_big_shl_bits(cfx_big_t* out, const cfx_big_t* x, unsigned s);
/* out = x >> s (0..63)  */
void cfx_big_shr_bits(cfx_big_t* out, const cfx_big_t* x, unsigned s);

void cfx_big_enable_cache(cfx_big_t* b);
void cfx_big_disable_cache(cfx_big_t* b);

/* Multiply by p^e by repeated squaring using small chunks to avoid overflow */
void cfx_big_expmul_prime(cfx_big_t* b, cfx_u64_t p, cfx_u64_t e);
/* out = n ^ p*/
void cfx_big_exp(cfx_big_t* out, const cfx_big_t* n, const cfx_big_t* p);
void cfx_big_exp_u64(cfx_big_t* out, const cfx_big_t* n, cfx_u64_t p);
/* out = (n^p) mod m */
void cfx_big_exp_mod(cfx_big_t* out, const cfx_big_t* n, const cfx_big_t* p, const cfx_big_t* m); 

void cfx_big_from_limbs(cfx_big_t* b, const cfx_u64_t* limbs, size_t n);
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t* f);
void cfx_big_to_fac(cfx_fac_t* f, const cfx_big_t* b);

char* cfx_big_to_str(const cfx_big_t* b, size_t *sz_out);
char* cfx_big_to_hex(const cfx_big_t* src, size_t *sz_out);
char* cfx_big_to_bin(const cfx_big_t* b, size_t *sz_out);
int cfx_big_from_str(cfx_big_t* b, const char* str);
int cfx_big_from_hex(cfx_big_t* out, const char* s);
int cfx_big_from_file(cfx_big_t* out, FILE* fp, int base);

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
 *   rr     : R^2 mod n  (used to convert x → xR mod n via MontMul)
 *   // Additional caches:
 *   // R1  : R mod n    (the Montgomery representation of 1)
 *   // mu  : floor(b^(2k)/n) for Barrett reduction of large inputs
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
    size_t    k;        /* limb count of n (R = 2^(64*k)) */
    cfx_u64_t  n0inv;    /* -n^{-1} mod 2^64 (requires n.limb[0] odd) */
    cfx_big_t rr;       /* R^2 mod n */
    cfx_big_t R1;       /* R mod n -- (Montgomery "1") */
} cfx_big_mont_ctx_t;

int cfx_big_mont_ctx_init(cfx_big_mont_ctx_t* ctx, const cfx_big_t* n);
void cfx_big_mont_ctx_free(cfx_big_mont_ctx_t* ctx);

/* Conversions between normal domain and Montgomery domain */
int cfx_big_mont_to   (cfx_big_t* out, const cfx_big_t* a,  const cfx_big_mont_ctx_t* ctx); /* out = a*R mod n */
int cfx_big_mont_from (cfx_big_t* out, const cfx_big_t* aR, const cfx_big_mont_ctx_t* ctx); /* out = aR*R^-1 */

/* Core Montgomery ops (operands/results in Montgomery domain) */
/**
 * Montgomery product (CIOS): out = a * b * R^{-1} mod n
 *
 * Base b = 2^64, R = b^k where k = ctx->k (number of 64-bit limbs of n).
 * Uses the Coarsely Integrated Operand Scanning (CIOS) algorithm:
 *   for i = 0..k-1:
 *     T += a * b[i]
 *     m = (T[0] * n0inv) mod b        // n0inv ≡ -n[0]^{-1} (mod b)
 *     T += m * n
 *     T >>= 64                        // drop least-significant limb
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

int cfx_big_mont_mul(cfx_big_t* out, const cfx_big_t* aR, const cfx_big_t* bR, const cfx_big_mont_ctx_t* ctx);
int cfx_big_modexp_binary(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* e, const cfx_big_mont_ctx_t* ctx);

static inline int cfx_big_mont_sqr(cfx_big_t* out, const cfx_big_t* aR, const cfx_big_mont_ctx_t* ctx) {
    return cfx_big_mont_mul(out, aR, aR, ctx);
}



/* Ergonomic one liners that use montgomery internally */
int cfx_big_sqr_mod (cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* n);
int cfx_big_modexp  (cfx_big_t* out, const cfx_big_t* base, const cfx_big_t* exp, const cfx_big_t* n);


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
