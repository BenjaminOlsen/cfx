/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * Compile with these (optional) flags to force limb sizes:
 *   -DCFX_FORCE_LIMB_32    or  -DCFX_FORCE_LIMB_64
 *   -DCFX_FORCE_ACC_STRUCT (use {lo,hi} accumulator struct even if native wide exists)
 *   -DCFX_FORCE_NO_UINT128 (pretend no __uint128_t)
 */

#ifndef CFX_ARITH_H
#define CFX_ARITH_H

#include "cfx/arch.h"

#include <stdint.h>
#include <inttypes.h>

/* -------- limb selection: 32 or 64 -------- */
#if defined(CFX_FORCE_LIMB_32)

typedef uint32_t cfx_limb_t;
#define CFX_LIMB_BITS 32

#elif defined(CFX_FORCE_LIMB_64)

typedef uint64_t cfx_limb_t;
#define CFX_LIMB_BITS 64

#else
/* Auto: prefer 64-bit limbs when native 128 exists; otherwise 32-bit. */
#  if CFX_HAS_UINT128
typedef uint64_t cfx_limb_t;
#    define CFX_LIMB_BITS 64
#  else
typedef uint32_t cfx_limb_t;
#    define CFX_LIMB_BITS 32
#  endif
#endif

#if (CFX_LIMB_BITS == 64)
#  define CFX_LIMB_DIGITS_DEC 19u /* log_10(2^64) = 19.26.. -> 10^19 fits in cfx_limb_t */
#  define CFX_LIMB_DIGITS_HEX 15u /* log_16(2^64) = 16. 16^15 = 2^60 fits in cfx_limb_t */
#  define CFX_LIMB_MAX UINT64_MAX
#elif (CFX_LIMB_BITS == 32)
#  define CFX_LIMB_DIGITS_DEC 9u /* log_10(2^32) = 9.63*/
#  define CFX_LIMB_DIGITS_HEX 8u /* log_16(2^32) = 8. */
#  define CFX_LIMB_MAX UINT32_MAX
#endif

/* -------- print formats -------- */
#if (CFX_LIMB_BITS == 64)

#  define CFX_SQRT_ACC_MAX 0xFFFFFFFFFFFFFFFFllu   /* floor(sqrt(UINT128_MAX+1)) - 1 */
#  define CFX_PRIxLIMB "%" PRIx64
#  define CFX_PRIuLIMB "%" PRIu64
#  define CFX_PRI0xLIMB "%016" PRIx64
#  define CFX_PRI0uLIMB "%016" PRIu64

#elif (CFX_LIMB_BITS == 32)

#  define CFX_SQRT_ACC_MAX 0xFFFFFFFFllu           /* floor(sqrt(UINT64_MAX+1)) - 1 */
#  define CFX_PRIxLIMB "%" PRIx32
#  define CFX_PRIuLIMB "%" PRIu32
#  define CFX_PRI0xLIMB "%08" PRIx32
#  define CFX_PRI0uLIMB "%08" PRIu32

#else
/* todo - choose defaults or error  */
#endif

/* -------- accumulator selection: 2x limb width --------
 * If limb=64 and native 128 exists and not forced to struct -> use __uint128_t.
 * Else use a portable {lo,hi} struct with limb-size halves.
 */
#if (CFX_LIMB_BITS == 64)

#  if CFX_HAS_UINT128 && !defined(CFX_FORCE_ACC_STRUCT)
    /* 64-bit limbs, native 128-bit accumulator */
    typedef __uint128_t cfx_acc_t;
#    define CFX_ACC_NATIVE 1
#  else
    /* 64-bit limbs, no __uint128_t -> use struct {lo,hi} */
    typedef struct { cfx_limb_t lo, hi; } cfx_acc_t;
#    define CFX_ACC_STRUCT 1
#  endif

#elif (CFX_LIMB_BITS == 32)

typedef uint64_t cfx_acc_t;
#define CFX_ACC_NATIVE 1

#else
#  error "Unsupported CFX_LIMB_BITS (expected 32 or 64)"
#endif

#ifdef _POSIX_THREADS
#  define CFX_HAS_PTHREAD 1
#endif

CFX_INLINE cfx_acc_t cfx_acc_zero(void) {
#if defined(CFX_ACC_NATIVE)
    return (cfx_acc_t)0;
#else
    cfx_acc_t a;
    a.lo = 0;
    a.hi = 0;
    return a;
#endif
}

CFX_INLINE cfx_acc_t cfx_acc_from_lo(cfx_limb_t lo) {
#if defined(CFX_ACC_NATIVE)
    return (cfx_acc_t)lo;
#else
    cfx_acc_t a;
    a.lo = lo;
    a.hi = 0;
    return a;
#endif
}

/* Low/High extractors (limb-sized) */
CFX_INLINE cfx_limb_t cfx_acc_lo(cfx_acc_t a) {
#if defined(CFX_ACC_NATIVE)
    return (cfx_limb_t)(a);
#else
    return a.lo;
#endif
}

CFX_INLINE cfx_limb_t cfx_acc_hi(cfx_acc_t a) {
#if defined(CFX_ACC_NATIVE)
    return (cfx_limb_t)(a >> CFX_LIMB_BITS);
#else
    return a.hi;
#endif
}

CFX_INLINE void cfx_acc_mul(cfx_acc_t* out, cfx_limb_t x, cfx_limb_t y) {
#if (CFX_LIMB_BITS == 64) && CFX_HAS_UINT128

    *out = (__uint128_t)x * (__uint128_t)y;

#elif (CFX_LIMB_BITS == 64)

    const uint64_t x0 = (uint32_t)x;
    const uint64_t x1 = (uint64_t)x >> 32;
    const uint64_t y0 = (uint32_t)y;
    const uint64_t y1 = (uint64_t)y >> 32;

    uint64_t p00 = x0 * y0;
    uint64_t p01 = x0 * y1;
    uint64_t p10 = x1 * y0;
    uint64_t p11 = x1 * y1;

    uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
    out->lo = (cfx_limb_t)((p00 & 0xffffffffu) | (mid << 32));
    out->hi = (cfx_limb_t)((p11 + (p01 >> 32)));

#else /* CFX_LIMB_BITS == 32 */

    *out = (uint64_t)x * (uint64_t)y;

#endif
}

/* -------- primitive: limb*limb -> {hi,lo} (each limb-sized) -------- */
CFX_INLINE void cfx_mul_wide(cfx_limb_t x, cfx_limb_t y,
                             cfx_limb_t* hi, cfx_limb_t* lo) {
#if (CFX_ACC_NATIVE)
    /* 64*64->128: use native when present (even if acc is struct) */
    cfx_acc_t p = ( (cfx_acc_t)x * (cfx_acc_t)y );
    *lo = (cfx_limb_t)p;
    *hi = (cfx_limb_t)(p >> CFX_LIMB_BITS);
#elif (CFX_LIMB_BITS == 64)
    /* Portable 64*64->128 via 32-bit halves */
    const uint64_t x0 = (uint32_t)x;
    const uint64_t x1 = (uint64_t)x >> 32;
    const uint64_t y0 = (uint32_t)y;
    const uint64_t y1 = (uint64_t)y >> 32;

    uint64_t p00 = x0 * y0;
    uint64_t p01 = x0 * y1;
    uint64_t p10 = x1 * y0;
    uint64_t p11 = x1 * y1;

    uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
    *lo = (cfx_limb_t)((p00 & 0xffffffffu) | (mid << 32));
    *hi = (cfx_limb_t)(p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32));
#else /* CFX_LIMB_BITS == 32 */
    /* 32*32->64 fits in uint64_t portably */
    uint64_t p = (uint64_t)x * (uint64_t)y;
    *lo = (cfx_limb_t)p;
    *hi = (cfx_limb_t)(p >> 32);
#endif
}

CFX_INLINE void cfx_acc_mul_eq(cfx_acc_t* a, const cfx_acc_t* b) {
#if defined(CFX_ACC_NATIVE)

    *a *= *b;

#else

    const cfx_limb_t a_lo = a->lo;
    const cfx_limb_t a_hi = a->hi;
    const cfx_limb_t b_lo = b->lo;
    const cfx_limb_t b_hi = b->hi;

    /* partial products */
    cfx_limb_t p0_hi, p0_lo;
    cfx_limb_t tmp_hi, p1_lo;
    cfx_limb_t tmp2_hi, p2_lo;

    /* p0 = a_lo * b_lo = p0_hi:B + p0_lo */
    cfx_mul_wide(a_lo, b_lo, &p0_hi, &p0_lo);

    /* p1 = a_hi * b_lo, we only need low limb of p1 */
    cfx_mul_wide(a_hi, b_lo, &tmp_hi, &p1_lo);

    /* p2 = a_lo * b_hi, we only need low limb of p2 */
    cfx_mul_wide(a_lo, b_hi, &tmp2_hi, &p2_lo);

    /* Full product:
         (a_hi B + a_lo)(b_hi B + b_lo) ≡ (p1_lo + p2_lo + p0_hi)·B + p0_lo  (mod B^2)
       High part of p1, p2 and a_hi*b_hi live in B^2 and above, discarded.
    */

    cfx_limb_t lo = p0_lo;
    cfx_limb_t hi = p0_hi;
    hi += p1_lo;
    hi += p2_lo;   /* natural limb wrap gives mod B */

    a->hi = hi;
    a->lo = lo;

#endif
}

/* -------- primitive: acc += low limb -------- */
CFX_INLINE void cfx_acc_add_lo(cfx_acc_t* acc, cfx_limb_t add) {
#if defined(CFX_ACC_NATIVE)
    *acc += (cfx_acc_t)add;
#else
    cfx_limb_t old = acc->lo;
    acc->lo = old + add;
    acc->hi += (acc->lo < old); /* carry */
#endif
}

/* -------- primitive: acc += high limb (i.e., add at +LIMB_BITS) -------- */
CFX_INLINE void cfx_acc_add_hi(cfx_acc_t* acc, cfx_limb_t add_hi) {
#if defined(CFX_ACC_NATIVE)
    *acc += ( (cfx_acc_t)add_hi << CFX_LIMB_BITS );
#else
    /* adding into the high half directly */
    cfx_limb_t old = acc->hi;
    acc->hi = old + add_hi;
    (void)old;
#endif
}

/* -------- primitive: acc += x*y -------- */
CFX_INLINE void cfx_acc_mac(cfx_acc_t* acc, cfx_limb_t x, cfx_limb_t y) {
#if defined(CFX_ACC_NATIVE)
    *acc += ( (cfx_acc_t)x * (cfx_acc_t)y );
#else
    cfx_limb_t hi, lo;
    cfx_mul_wide(x, y, &hi, &lo);
    cfx_acc_add_lo(acc, lo);
    cfx_acc_add_hi(acc, hi);
#endif
}

CFX_INLINE int cfx_acc_lt(const cfx_acc_t* u, const cfx_acc_t* v) {
#if defined(CFX_ACC_NATIVE)
    return *u < *v;
#else
    return u->hi == v->hi ? u->lo < v->lo : u->hi < v->hi;
#endif
}

CFX_INLINE int cfx_acc_le(const cfx_acc_t* u, const cfx_acc_t* v) {
#if defined(CFX_ACC_NATIVE)
    return *u <= *v;
#else
    return u->hi == v->hi ? u->lo <= v->lo : u->hi <= v->hi;
#endif
}

CFX_INLINE int cfx_acc_eq(const cfx_acc_t* u, const cfx_acc_t* v) {
#if defined(CFX_ACC_NATIVE)
    return *u == *v;
#else
    return u->hi == v->hi && u->lo == v->lo;
#endif
}

CFX_INLINE void cfx_acc_divrem(cfx_acc_t* q, cfx_acc_t* r,
                               const cfx_acc_t* u, const cfx_acc_t* v) {
#if defined(CFX_ACC_NATIVE)

    *q = *u / *v;
    *r = *u % *v;

#else
    /* Portable 2-limb / 2-limb unsigned division:
       q = u / v, r = u % v, where each is (hi, lo) in base 2^CFX_LIMB_BITS.
       Works for CFX_LIMB_BITS == 32 or 64.
    */

    cfx_acc_t U = *u;
    cfx_acc_t V = *v;

    /* division by zero */
    if (V.hi == 0 && V.lo == 0) {
        /* Here: q = 0, r = u */
        q->hi = 0;
        q->lo = 0;
        *r = U;
        return;
    }

    /* If U < V then quotient = 0, remainder = U */
    if (U.hi < V.hi || (U.hi == V.hi && U.lo < V.lo)) {
        q->hi = 0;
        q->lo = 0;
        *r = U;
        return;
    }

    cfx_acc_t Q = cfx_acc_zero();
    cfx_acc_t R = cfx_acc_zero();

    const int TOTAL_BITS = 2 * CFX_LIMB_BITS;

    for (int i = TOTAL_BITS - 1; i >= 0; --i) {
        /* R <<= 1 */
        cfx_limb_t carry = R.lo >> (CFX_LIMB_BITS - 1);
        R.lo <<= 1;
        R.hi = (R.hi << 1) | carry;

        /* Bring in next bit from U (MSB-first) */
        cfx_limb_t bit;
        if (i >= CFX_LIMB_BITS) {
            bit = (U.hi >> (i - CFX_LIMB_BITS)) & (cfx_limb_t)1;
        } else {
            bit = (U.lo >> i) & (cfx_limb_t)1;
        }
        R.lo |= bit;

        /* If R >= V, subtract and set quotient bit */
        if (R.hi > V.hi || (R.hi == V.hi && R.lo >= V.lo)) {
            /* R -= V (with borrow) */
            cfx_limb_t old_lo = R.lo;
            R.lo -= V.lo;
            R.hi -= V.hi + (R.lo > old_lo ? 1 : 0);

            /* Set bit i in Q */
            if (i >= CFX_LIMB_BITS) {
                Q.hi |= (cfx_limb_t)1 << (i - CFX_LIMB_BITS);
            } else {
                Q.lo |= (cfx_limb_t)1 << i;
            }
        }
    }

    *q = Q;
    *r = R;

#endif
}

#endif /* CFX_ARITH_H */
