/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * cfx/sbig.h - Signed big integer wrapper
 *
 * Provides signed arithmetic on top of the unsigned cfx_big_t.
 *
 * Design: magnitude + sign representation
 *   - mag: the absolute value (always non-negative)
 *   - sign: -1 (negative), 0 (zero), or +1 (positive)
 *   - Zero always has sign = 0, regardless of operations
 */

#ifndef CFX_SBIG_H
#define CFX_SBIG_H

#include "cfx/big.h"
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Signed big integer type.
 * Note: typedef'd as cfx_sbig_t in cfx/big.h for forward declaration.
 */
struct cfx_sbig {
    cfx_big_t mag;   /* magnitude (absolute value) */
    int8_t sign;     /* -1 = negative, 0 = zero, +1 = positive */
};

/*
 * Lightweight signed reference to an existing unsigned bigint.
 * Useful when you don't want to copy the magnitude.
 * The referenced cfx_big_t must outlive this reference.
 */
typedef struct {
    const cfx_big_t* mag;  /* pointer to magnitude (borrowed, not owned) */
    int8_t sign;           /* -1, 0, or +1 */
} cfx_sref_t;


void cfx_sbig_init(cfx_sbig_t* s);
void cfx_sbig_free(cfx_sbig_t* s);
void cfx_sbig_clear(cfx_sbig_t* s);  /* set to zero without freeing */
void cfx_sbig_assign(cfx_sbig_t* dst, const cfx_sbig_t* src);
void cfx_sbig_assign_big(cfx_sbig_t* dst, const cfx_big_t* mag, int8_t sign);
void cfx_sbig_from_i64(cfx_sbig_t* s, int64_t val);
void cfx_sbig_from_u64(cfx_sbig_t* s, uint64_t val);

static inline int8_t cfx_sbig_sign(const cfx_sbig_t* s) { return s->sign; }
static inline int cfx_sbig_is_negative(const cfx_sbig_t* s) { return s->sign < 0; }
static inline int cfx_sbig_is_zero(const cfx_sbig_t* s) { return s->sign == 0; }
static inline int cfx_sbig_is_positive(const cfx_sbig_t* s) { return s->sign > 0; }
static inline const cfx_big_t* cfx_sbig_mag(const cfx_sbig_t* s) { return &s->mag; }

void cfx_sbig_neg_eq(cfx_sbig_t* s);
void cfx_sbig_abs_eq(cfx_sbig_t* s);
void cfx_sbig_neg(cfx_sbig_t* out, const cfx_sbig_t* a);
void cfx_sbig_abs(cfx_sbig_t* out, const cfx_sbig_t* a);

int cfx_sbig_cmp(const cfx_sbig_t* a, const cfx_sbig_t* b);
static inline int cfx_sbig_eq(const cfx_sbig_t* a, const cfx_sbig_t* b) {
    return cfx_sbig_cmp(a, b) == 0;
}

void cfx_sbig_add(cfx_sbig_t* out, const cfx_sbig_t* a, const cfx_sbig_t* b);
void cfx_sbig_sub(cfx_sbig_t* out, const cfx_sbig_t* a, const cfx_sbig_t* b);
void cfx_sbig_mul(cfx_sbig_t* out, const cfx_sbig_t* a, const cfx_sbig_t* b);
void cfx_sbig_divrem(cfx_sbig_t* q, cfx_sbig_t* rem,
                     const cfx_sbig_t* a, const cfx_sbig_t* b);

void cfx_sbig_add_eq(cfx_sbig_t* a, const cfx_sbig_t* b);
void cfx_sbig_sub_eq(cfx_sbig_t* a, const cfx_sbig_t* b);
void cfx_sbig_mul_eq(cfx_sbig_t* a, const cfx_sbig_t* b);

/* caller must free returned string */
char* cfx_sbig_dec_alloc(const cfx_sbig_t* s, size_t* len_out);
int cfx_sbig_from_dec(cfx_sbig_t* out, const char* str);
int cfx_sbig_from_hex(cfx_sbig_t* out, const char* str);

void cfx_sbig_normalize(cfx_sbig_t* s);
void cfx_sbig_swap(cfx_sbig_t* a, cfx_sbig_t* b);

/* Extended GCD for signed big integers.
 * Computes g = gcd(|a|, |b|) and coefficients x, y such that a*x + b*y = g.
 * g is always non-negative. Any of g, x, y can be NULL if not needed. */
void cfx_sbig_xgcd(cfx_sbig_t* g, cfx_sbig_t* x, cfx_sbig_t* y,
                   const cfx_sbig_t* a, const cfx_sbig_t* b);

#ifdef __cplusplus
}
#endif

#endif /* CFX_SBIG_H */
