#ifndef CFX_BIG_H
#define CFX_BIG_H

#include "cfx/fac.h"
#include "cfx/fmt.h"
#include "cfx/algo.h"
#include "cfx/types.h"

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>

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
    uint64_t* limb; /* the 'digits' of the number with base BIG_BASE*/
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
    uint64_t bnd;           /* proven “B-smoothness” bound:
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
int cfx_big_eq_sm(const cfx_big_t* b, uint64_t n);
int cfx_big_eq(const cfx_big_t* b1, const cfx_big_t* b2);
int cfx_big_cmp(const cfx_big_t* a, const cfx_big_t* b);
void cfx_big_swap(cfx_big_t* a, cfx_big_t* b);

void cfx_big_from_u64(cfx_big_t* b, uint64_t v);
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
void cfx_big_add_sm(cfx_big_t* b, uint64_t n);
void cfx_big_sub(cfx_big_t* a, const cfx_big_t* b);
void cfx_big_sub_sm(cfx_big_t* b, uint64_t n);
void cfx_big_mul_sm(cfx_big_t* b, uint64_t m);

/* Knuth's long division from TAOCP Vol. 2, 4.3.1 */
int cfx_big_divrem(cfx_big_t* q, cfx_big_t* r, const cfx_big_t* u, const cfx_big_t* v);

int cfx_big_div_out(cfx_big_t* q, const cfx_big_t* n, const cfx_big_t* d);
int cfx_big_mod_out(cfx_big_t* r, const cfx_big_t* n, const cfx_big_t* d);
/* In-place: b := floor(b/d); optional remainder r. Alias-safe for any combination. */
int cfx_big_div_eq(cfx_big_t* b, const cfx_big_t* d, cfx_big_t* r /*nullable*/);

uint64_t cfx_big_div_sm(cfx_big_t* b, uint64_t d);
uint32_t cfx_big_div_sm_u32(cfx_big_t* b, uint32_t d);
uint64_t cfx_big_mod_sm(cfx_big_t* b, uint64_t m);

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
void cfx_big_powmul_prime(cfx_big_t* b, uint64_t p, uint64_t e);
/* out = n ^ p*/
void cfx_big_pow(cfx_big_t* out, const cfx_big_t* n, const cfx_big_t* p);

void cfx_big_from_limbs(cfx_big_t* b, const uint64_t* limbs, size_t n);
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t* f);
void cfx_big_to_fac(cfx_fac_t* f, const cfx_big_t* b);

char* cfx_big_to_str(const cfx_big_t* b, size_t *sz_out);
char* cfx_big_to_hex(const cfx_big_t* src, size_t *sz_out);
char* cfx_big_to_bin(const cfx_big_t* b, size_t *sz_out);
int cfx_big_from_str(cfx_big_t* b, const char* str);
int cfx_big_from_hex(cfx_big_t* out, const char* s);
int cfx_big_from_file(cfx_big_t* out, FILE* fp, int base);

/* ================== Montgomery ================== */

typedef struct {
    cfx_big_t n;        /* modulus (odd), normalized */
    size_t    k;        /* limb count of n (R = 2^(64*k)) */
    uint64_t  n0inv;    /* -n^{-1} mod 2^64 (requires n.limb[0] odd) */
    cfx_big_t rr;       /* R^2 mod n */
} cfx_big_mont_ctx_t;

int cfx_big_mont_ctx_init(cfx_big_mont_ctx_t* ctx, const cfx_big_t* n);
void cfx_big_mont_ctx_free(cfx_big_mont_ctx_t* ctx);

/* Conversions between normal domain and Montgomery domain */
int cfx_big_mont_to   (cfx_big_t* out, const cfx_big_t* a,  const cfx_big_mont_ctx_t* ctx); /* out = a*R mod n */
int cfx_big_mont_from (cfx_big_t* out, const cfx_big_t* aR, const cfx_big_mont_ctx_t* ctx); /* out = aR*R^-1 */

/* Core Montgomery ops (operands/results in Montgomery domain) */
int cfx_big_mont_mul1(cfx_big_t* out, const cfx_big_t* aR, const cfx_big_t* bR, const cfx_big_mont_ctx_t* ctx);
int cfx_big_mont_mul2(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, const cfx_big_mont_ctx_t* ctx);
int cfx_big_mont_mul(cfx_big_t* out, const cfx_big_t* aR, const cfx_big_t* bR, const cfx_big_mont_ctx_t* ctx);
static inline int cfx_big_mont_sqr(cfx_big_t* out, const cfx_big_t* aR, const cfx_big_mont_ctx_t* ctx) {
    return cfx_big_mont_mul(out, aR, aR, ctx);
}

/* Ergonomic one-liners that hide the context internally */
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
