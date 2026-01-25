/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_MONT_H
#define CFX_MONT_H

/*
 * Montgomery space for 64 bit moduli.
 */

#include "cfx/arch.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Montgomery context for repeated modular multiplications with same modulus.
 */
typedef struct {
    uint64_t n;      /* the modulus (must be odd) */
    uint64_t n_inv;  /* -n^(-1) mod 2^64 */
    uint64_t r2;     /* R^2 mod n, where R = 2^64 */
    uint64_t r;      /* R mod n (for fallback path) */
    uint64_t r_inv;  /* R^(-1) mod n (for fallback path) */
} cfx_mont64_t;

/* Compute modular inverse: a^(-1) mod m using extended Euclidean algorithm
 * Uses careful unsigned arithmetic to handle full 64-bit range.
 * Returns 0 if inverse doesn't exist (gcd(a, m) != 1) */
CFX_INLINE uint64_t cfx_modinv64(uint64_t a, uint64_t m) {
    if (m == 1 || a == 0) return 0;
    a = a % m;
    if (a == 0) return 0;
    if (a == 1) return 1;

    /* Extended Euclidean algorithm with unsigned arithmetic.
     * We track coefficients as (value, negative_flag) pairs to avoid signed overflow.
     * x1 = a^(-1) mod m when the algorithm terminates with gcd = 1.
     */
    uint64_t u = a, v = m;
    uint64_t x1 = 1, x2 = 0;
    int x1_neg = 0, x2_neg = 0;

    while (u != 1 && v != 1) {
        while ((u & 1) == 0) {
            u >>= 1;
            if ((x1 & 1) == 0) {
                x1 >>= 1;
            } else {
                /* x1 = (x1 + m) / 2, but we need to handle potential overflow */
                uint64_t sum = x1 + m;
                int overflow = (sum < x1);  /* Detect wrap-around */
                x1 = sum >> 1;
                if (overflow) x1 |= (1ULL << 63);  /* Carry the overflow bit */
            }
        }
        while ((v & 1) == 0) {
            v >>= 1;
            if ((x2 & 1) == 0) {
                x2 >>= 1;
            } else {
                uint64_t sum = x2 + m;
                int overflow = (sum < x2);
                x2 = sum >> 1;
                if (overflow) x2 |= (1ULL << 63);
            }
        }
        if (u >= v) {
            u -= v;
            /* x1 = x1 - x2 (mod m), handling sign */
            if (x1_neg == x2_neg) {
                if (x1 >= x2) {
                    x1 = x1 - x2;
                } else {
                    x1 = x2 - x1;
                    x1_neg = !x1_neg;
                }
            } else {
                x1 = x1 + x2;
                if (x1 >= m) x1 -= m;
            }
        } else {
            v -= u;
            /* x2 = x2 - x1 (mod m), handling sign */
            if (x2_neg == x1_neg) {
                if (x2 >= x1) {
                    x2 = x2 - x1;
                } else {
                    x2 = x1 - x2;
                    x2_neg = !x2_neg;
                }
            } else {
                x2 = x2 + x1;
                if (x2 >= m) x2 -= m;
            }
        }
    }

    uint64_t result = (u == 1) ? x1 : x2;
    int result_neg = (u == 1) ? x1_neg : x2_neg;

    if (result_neg) {
        result = m - result;
    }
    return result % m;
}

/* Compute -n^(-1) mod 2^64 using Newton's method */
CFX_INLINE uint64_t cfx_mont64_inv(uint64_t n) {
    /* n must be odd; find x such that n*x ≡ -1 (mod 2^64) */
    /* Newton iteration: x = x * (2 + n*x) converges for odd n */
    uint64_t x = n;  /* initial guess: n*n ≡ 1 (mod 8) for odd n */
    x *= 2 - n * x;  /* 4 bits */
    x *= 2 - n * x;  /* 8 bits */
    x *= 2 - n * x;  /* 16 bits */
    x *= 2 - n * x;  /* 32 bits */
    x *= 2 - n * x;  /* 64 bits */
    return (uint64_t)(-(int64_t)x); /* return -x = -n^(-1) mod 2^64 */
}

/* Initialize Montgomery context for modulus n (must be odd) */
CFX_INLINE void cfx_mont64_init(cfx_mont64_t *ctx, uint64_t n) {
    ctx->n = n;
    ctx->n_inv = cfx_mont64_inv(n);
    /* Compute R mod n where R = 2^64 */
    /* -n as uint64_t gives 2^64 - n, then mod n gives R mod n */
    ctx->r = (uint64_t)(-((int64_t)n)) % n;
    /* Compute R^2 mod n = (R mod n)^2 mod n */
    ctx->r2 = cfx_mulmod64(ctx->r, ctx->r, n);
    /* Compute R^(-1) mod n for fallback path (n is odd so gcd(R,n)=1) */
    ctx->r_inv = cfx_modinv64(ctx->r, n);
}

/* Convert a to Montgomery form: aR mod n */
/* Uses MonPro(a, R²) = a * R² * R⁻¹ mod n = a*R mod n */
CFX_INLINE uint64_t cfx_mont64_to(const cfx_mont64_t *ctx, uint64_t a) {
    a = a % ctx->n;
#if defined(CFX_128BIT_NATIVE)
    __uint128_t t = (__uint128_t)a * ctx->r2;
    uint64_t m = (uint64_t)t * ctx->n_inv;
    __uint128_t mn = (__uint128_t)m * ctx->n;
    t += mn;
    int overflow = (t < mn);  /* detect 128-bit overflow */
    uint64_t u = (uint64_t)(t >> 64);
    return (overflow || u >= ctx->n) ? (u - ctx->n) : u;
#elif defined(CFX_128BIT_MSVC)
    uint64_t t_hi, t_lo;
    t_lo = _umul128(a, ctx->r2, &t_hi);
    uint64_t m = t_lo * ctx->n_inv;
    uint64_t mn_hi, mn_lo;
    mn_lo = _umul128(m, ctx->n, &mn_hi);
    /* 128-bit add: t += m*n, tracking overflow into bit 128 */
    uint64_t carry = 0;
    t_lo += mn_lo;
    if (t_lo < mn_lo) carry = 1;
    uint64_t old_t_hi = t_hi;
    t_hi += mn_hi + carry;
    /* overflow if t_hi wrapped around */
    int overflow = (t_hi < old_t_hi) || (t_hi < mn_hi);
    /* If overflow, result >= 2^64 > n, so subtract n */
    /* Also subtract if t_hi >= n (could be up to 2n-1) */
    return (overflow || t_hi >= ctx->n) ? (t_hi - ctx->n) : t_hi;
#else
    /* Fallback without 128-bit: compute a*R mod n using precomputed R mod n */
    return cfx_mulmod64(a, ctx->r, ctx->n);
#endif
}

/* Convert from Montgomery form: a/R mod n */
CFX_INLINE uint64_t cfx_mont64_from(const cfx_mont64_t *ctx, uint64_t aR) {
    /* MonPro(aR, 1) = a*R*1/R mod n = a mod n */
#if defined(CFX_128BIT_NATIVE)
    __uint128_t t = aR;
    uint64_t m = (uint64_t)t * ctx->n_inv;
    __uint128_t mn = (__uint128_t)m * ctx->n;
    t += mn;
    int overflow = (t < mn);
    uint64_t u = (uint64_t)(t >> 64);
    return (overflow || u >= ctx->n) ? (u - ctx->n) : u;
#elif defined(CFX_128BIT_MSVC)
    uint64_t t_lo = aR, t_hi = 0;
    uint64_t m = t_lo * ctx->n_inv;
    uint64_t mn_hi, mn_lo;
    mn_lo = _umul128(m, ctx->n, &mn_hi);
    /* t += m*n, tracking overflow */
    uint64_t carry = 0;
    t_lo += mn_lo;
    if (t_lo < mn_lo) carry = 1;
    t_hi += mn_hi + carry;
    /* For from(), t_hi starts at 0, so overflow only if mn_hi + carry overflows */
    int overflow = (t_hi < mn_hi);
    return (overflow || t_hi >= ctx->n) ? (t_hi - ctx->n) : t_hi;
#else
    /* Fallback: multiply by R^(-1) to convert back from Montgomery form */
    return cfx_mulmod64(aR, ctx->r_inv, ctx->n);
#endif
}

/* Montgomery multiplication: MonPro(aR, bR) = (a*b*R) mod n */
CFX_INLINE uint64_t cfx_mont64_mul(const cfx_mont64_t *ctx, uint64_t aR, uint64_t bR) {
#if defined(CFX_128BIT_NATIVE)
    __uint128_t t = (__uint128_t)aR * bR;
    uint64_t m = (uint64_t)t * ctx->n_inv;
    __uint128_t mn = (__uint128_t)m * ctx->n;
    t += mn;
    int overflow = (t < mn);
    uint64_t u = (uint64_t)(t >> 64);
    return (overflow || u >= ctx->n) ? (u - ctx->n) : u;
#elif defined(CFX_128BIT_MSVC)
    uint64_t t_hi, t_lo;
    t_lo = _umul128(aR, bR, &t_hi);
    uint64_t m = t_lo * ctx->n_inv;
    uint64_t mn_hi, mn_lo;
    mn_lo = _umul128(m, ctx->n, &mn_hi);
    /* t += m*n, tracking overflow */
    uint64_t carry = 0;
    t_lo += mn_lo;
    if (t_lo < mn_lo) carry = 1;
    uint64_t old_t_hi = t_hi;
    t_hi += mn_hi + carry;
    int overflow = (t_hi < old_t_hi) || (t_hi < mn_hi);
    return (overflow || t_hi >= ctx->n) ? (t_hi - ctx->n) : t_hi;
#else
    /* Fallback: regular mulmod then convert */
    uint64_t a = cfx_mont64_from(ctx, aR);
    uint64_t b = cfx_mont64_from(ctx, bR);
    return cfx_mont64_to(ctx, cfx_mulmod64(a, b, ctx->n));
#endif
}

/* Montgomery squaring (slightly optimized) */
CFX_INLINE uint64_t cfx_mont64_sqr(const cfx_mont64_t *ctx, uint64_t aR) {
    return cfx_mont64_mul(ctx, aR, aR);
}

/* Montgomery add: (aR + bR) mod n (no division needed) */
CFX_INLINE uint64_t cfx_mont64_add(const cfx_mont64_t *ctx, uint64_t aR, uint64_t bR) {
    uint64_t sum = aR + bR;
    /* Handle overflow and reduction */
    if (sum < aR || sum >= ctx->n) {
        sum -= ctx->n;
    }
    return sum;
}

/* Montgomery subtract: (aR - bR) mod n */
CFX_INLINE uint64_t cfx_mont64_sub(const cfx_mont64_t *ctx, uint64_t aR, uint64_t bR) {
    return (aR >= bR) ? (aR - bR) : (ctx->n - (bR - aR));
}

/* Montgomery modular exponentiation */
CFX_INLINE uint64_t cfx_mont64_pow(const cfx_mont64_t *ctx, uint64_t aR, uint64_t exp) {
    uint64_t result = cfx_mont64_to(ctx, 1);  /* 1 in Montgomery form = R mod n */
    while (exp > 0) {
        if (exp & 1) result = cfx_mont64_mul(ctx, result, aR);
        aR = cfx_mont64_sqr(ctx, aR);
        exp >>= 1;
    }
    return result;
}

#ifdef __cplusplus
}
#endif

#endif /* CFX_MONT_H */
