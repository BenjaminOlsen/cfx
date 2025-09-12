#ifndef CFX_BIG_H
#define CFX_BIG_H

#include "cfx/fac.h"
#include "cfx/types.h"

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CFX_BIG_PRINT(bg) printf("b.n %zu: b.cap: %zu, b %s\n", bg.n, bg.cap, cfx_big_to_str(&bg, NULL))

#define CFX_BIG_PRINT_LIMBS(bg) \
    do { \
        for (size_t i = 0; i < bg.n; ++i) { \
            printf("limb[%zu] = %llu\n", i, bg.limb[i]); \
        } \
    } while (0)

typedef enum {
    CFX_FAC_NONE,
    CFX_FAC_PARTIAL,
    CFX_FAC_FULL
} cfx_fac_state_e;

typedef struct cfx_fac_cache cfx_fac_cache_t;

/* minimal big int struct */
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
void cfx_big_reserve(cfx_big_t* b, size_t need);
void cfx_big_copy(cfx_big_t* dst, const cfx_big_t* src);

void cfx_big_set_val(cfx_big_t* b, uint64_t v);
cfx_big_t cfx_big_mul(const cfx_big_t* b, const cfx_big_t* m);
void cfx_big_mul_eq(cfx_big_t* b, const cfx_big_t* m);

cfx_big_t cfx_big_sq(const cfx_big_t* b);
void cfx_big_sq_eq(cfx_big_t* b);
int cfx_big_div_eq(cfx_big_t* b, const cfx_big_t* d, cfx_big_t* r); /* b /= d; r remainder. */

void cfx_big_add_sm(cfx_big_t* b, uint64_t n);
void cfx_big_sub_sm(cfx_big_t* b, uint64_t n);
void cfx_big_mul_sm(cfx_big_t* b, uint64_t m);

uint64_t cfx_big_div_eq_sm(cfx_big_t* b, uint64_t d);
cfx_big_t cfx_big_div_sm(const cfx_big_t* b, uint64_t d, uint64_t* r);
uint32_t cfx_big_div_eq_sm_u32(cfx_big_t* b, uint32_t d);
cfx_big_t cfx_big_div_sm_u32(const cfx_big_t* b, uint32_t d);
uint64_t cfx_big_mod_eq_sm(cfx_big_t* b, uint64_t m);
cfx_big_t cfx_big_mod_sm(const cfx_big_t* b, uint64_t m);

int cfx_big_is_zero(const cfx_big_t* b);
int cfx_big_eq_sm(const cfx_big_t* b, uint64_t n);
int cfx_big_eq(const cfx_big_t* b1, const cfx_big_t* b2);
void cfx_big_swap(cfx_big_t* a, cfx_big_t* b);

void cfx_big_enable_fac(cfx_big_t* b);
void cfx_big_disable_fac(cfx_big_t* b);

/* Multiply by p^e by repeated squaring using small chunks to avoid overflow */
void cfx_big_powmul_prime(cfx_big_t* b, uint64_t p, uint64_t e);

void cfx_big_from_limbs(cfx_big_t* b, const uint64_t* limbs, size_t n);
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
