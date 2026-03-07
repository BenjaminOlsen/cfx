/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include "cfx/algo.h"

/* Platform-specific includes */
#if defined(_WIN32) || defined(_WIN64)
    #include <intrin.h>  /* MSVC intrinsics for _mulx_u64, etc. */
#endif


static inline void _mul_sm_fast(cfx_big_t *b, cfx_limb_t m) {
    size_t n = b->n;
    cfx_big_reserve(b, n + 1);
    cfx_limb_t *p = b->limb;

#if (CFX_USE_X86_INTRINSICS == 1) && defined(__BMI2__) && (CFX_LIMB_BITS == 64)
    cfx_limb_t carry = 0;
    for (size_t i = 0; i < n; ++i) {
        unsigned long long lo, hi;

        /* lo = (p[i] * m) low 64, hi = high 64 — flags NOT clobbered */
        lo = _mulx_u64(p[i], m, &hi);

        /* add previous carry into lo, get carry-out in c1 */
        unsigned char c1 = _addcarry_u64(0, lo, carry, &lo);
        hi += c1;
        p[i] = (cfx_limb_t)lo;
        carry = (cfx_limb_t)hi;
    }

    p[n] = carry;
    b->n = n + (carry != 0);
#else /* Portable fallback (cfx_acc_t) */
    cfx_limb_t carry = 0;
    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t t;
        cfx_acc_mul(&t, p[i], m);
        cfx_acc_add_lo(&t, carry);
        p[i] = cfx_acc_lo(t);
        carry = cfx_acc_hi(t);
    }
    p[n] = carry;
    b->n = n + (carry != 0);
#endif
}

/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
void cfx_big_expmul_prime(cfx_big_t *b, cfx_limb_t p, cfx_limb_t e) {

    /* Find largest t so that p^t fits in 32 bits -> p^2t fits in 64 */
    cfx_limb_t t = 1;
    cfx_acc_t acc = cfx_acc_from_lo(p);
    const cfx_acc_t lim = cfx_acc_from_lo(CFX_SQRT_ACC_MAX);

    cfx_acc_t quot, rem;
    cfx_acc_divrem(&quot, &rem, &lim, &acc);
    while (cfx_acc_le(&acc, &quot)) {
        cfx_acc_mul_eq(&acc, &acc);
        t *= 2u;
        cfx_acc_divrem(&quot, &rem, &lim, &acc);
    }
    /* now t is the largest s.t. p^t <= lim */
    CFX_PRINT_DBG("multiplying by " CFX_PRIuLIMB "^" CFX_PRIuLIMB " by breaking it into (" CFX_PRIuLIMB "^" CFX_PRIuLIMB ")^(" CFX_PRIuLIMB "/" CFX_PRIuLIMB ") \n", p, e, p, t, e, t);

    /* Now multiply by (p^t)^(e/t) and then the remainder */
    /* Compute p^t */
    /* compute power safely */
    cfx_limb_t pow_t = p;

    /* fast power to get pow_t = p^t by using binary expansion of t:
       p^t = p^(b_i*2^i) * p^(b_(i-1)*2^(i-1) * ... */
    pow_t = 1;
    acc = p;
    cfx_limb_t tt = t;
    while (tt) {
        if (tt & 1u) pow_t = (cfx_limb_t)(pow_t * acc);
        tt >>= 1u;
        acc = acc*acc;
    }
    cfx_limb_t q = e / t;
    cfx_limb_t r = e % t;

    for (cfx_limb_t i = 0; i < q; i++) _mul_sm_fast(b, pow_t);

    /* multiply remainder by binary exponentiation (still fits u32) */
    cfx_limb_t rempow = 1;
    acc = p;
    cfx_limb_t rr = r;

    while (rr) {
        if (rr & 1u) rempow = (cfx_limb_t)(rempow * acc);
        rr >>= 1u;
        acc = acc*acc;
    }
    if (rempow != 1) _mul_sm_fast(b, rempow);
}

static inline void cfx_big_mul_small_inplace(cfx_big_t *out, cfx_limb_t m) {
    if (m == 0) {
        cfx_big_from_limb(out, 0); return;
    }
    if (m == 1) return;
    _mul_sm_fast(out, m);
}

/* out = p^e (p is a limb prime) */
void cfx_big_pow_sm(cfx_big_t *out, cfx_limb_t p, cfx_limb_t e) {
    cfx_big_from_limb(out, 1);
    if (e == 0) return;


    if (p == 2) {
        /* out = 1 << e */
        cfx_big_t one;
        cfx_big_init(&one);
        cfx_big_from_limb(&one, 1);
        cfx_big_shl_bits(out, &one, (unsigned)e);
        cfx_big_free(&one);
        return;
    }

    if (e <= 3) {
        for (cfx_limb_t i = 0; i < e; i++) cfx_big_mul_small_inplace(out, p);
        return;
    }

    /* Binary exponentiation */
    cfx_big_t base;
    cfx_big_init(&base);
    cfx_big_from_limb(&base, p);

    cfx_limb_t exp = e;
    while (exp) {
        if (exp & 1u) {
            cfx_big_mul_auto(out, &base); /* out *= base */
        }
        exp >>= 1u;
        if (exp) {
            /* base = base^2 */
            cfx_big_mul_auto(&base, &base);
        }
    }

    cfx_big_free(&base);
}

/* out = a + b */
void cfx_big_add(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
    /* handle aliasing: if out aliases a or b, work on a temp */
    if (out == a) {
        cfx_big_add_eq(out, b);
        return;
    }
    if (out == b) {
        cfx_big_add_eq(out, a);
        return;
    }

    /* out is distinct from a and b */
    cfx_big_assign(out, a);
    cfx_big_add_eq(out, b);
}

/* b += a */
void cfx_big_add_eq(cfx_big_t *b, const cfx_big_t *a) {
    cfx_limb_t carry = 0;
    size_t i = 0;

    while (i < a->n || carry) {
        cfx_big_reserve(b, i + 1);

        cfx_acc_t bi = (i < b->n) ? b->limb[i] : 0;
        cfx_acc_t ai = (i < a->n) ? a->limb[i] : 0;
        cfx_acc_t s  = bi + ai + carry;

        if (i >= b->n) b->n = i + 1;      /* we're extending b */
        b->limb[i] = (cfx_limb_t)s;
        carry      = (cfx_limb_t)(s >> CFX_LIMB_BITS);
        ++i;
    }

    /* If a->n > b->n and carry ended at 0, we may still need to bump b->n */
    if (i > b->n) b->n = i;
}


void cfx_big_add_sm_eq(cfx_big_t *b, cfx_limb_t n) {
    if (n == 0) return;
    if (b->n == 0) {
        cfx_big_from_limb(b, n);
        return;
    }

    size_t i = 0;
    while (i < b->n) {
        cfx_acc_t s = (cfx_acc_t)b->limb[i] + n;
        b->limb[i] = (cfx_limb_t)s;
        n = (cfx_limb_t)(s >> CFX_LIMB_BITS);
        if (n == 0) {
            return;
        }
        ++i;
    }

    /* carry left; append one limb */
    cfx_big_reserve(b, b->n + 1);
    b->limb[b->n++] = n;      /* n < 2^64 here */
}

/* out = a - b (assumes a >= b) */
void cfx_big_sub(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
    /* handle aliasing */
    if (out == a) {
        cfx_big_sub_eq(out, b);
        return;
    }
    if (out == b) {
        /* out = a - out, need temp */
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_assign(&tmp, a);
        cfx_big_sub_eq(&tmp, b);
        cfx_big_swap(out, &tmp);
        cfx_big_free(&tmp);
        return;
    }

    /* out is distinct from a and b */
    cfx_big_assign(out, a);
    cfx_big_sub_eq(out, b);
}

void cfx_big_sub_eq(cfx_big_t *a, const cfx_big_t *b) {
    /* assumes a >= b; subtract b from a */
    cfx_acc_t borrow = 0;
    size_t i = 0, n = b->n;
    for (; i < n; ++i) {
        cfx_acc_t ai = a->limb[i];
        cfx_acc_t s  = (cfx_acc_t)b->limb[i] + borrow;
        cfx_acc_t r  = ai - s;
        a->limb[i] = (cfx_limb_t)r;
        borrow = (ai < s);
    }
    while (borrow && i < a->n) {
        cfx_limb_t ai = a->limb[i];
        cfx_limb_t r  = ai - 1u;
        a->limb[i++] = r;
        borrow = (ai == 0u);
    }
    cfx_big_trim(a);
}

void cfx_big_sub_sm_eq(cfx_big_t *b, cfx_limb_t n) {
    if (n == 0 || b->n == 0) return;

    /* Subtract n from limb[0], propagate borrow if needed */
    cfx_limb_t borrow = (b->limb[0] < n) ? 1 : 0;
    b->limb[0] -= n;

    /* Propagate borrow through remaining limbs */
    size_t i = 1;
    while (borrow && i < b->n) {
        cfx_limb_t old = b->limb[i];
        b->limb[i] -= 1;
        borrow = (old == 0);
        ++i;
    }

    /* If borrow remains, we had underflow (b < n). This is a caller bug.
     * The result is now garbage (wrapped), but we don't crash. */
    assert(!borrow && "cfx_big_sub_sm_eq: underflow (b < n)");

    cfx_big_trim(b);
}

/* subtract uint64_t in place - handles 32-bit limb builds */
void cfx_big_sub_u64_eq(cfx_big_t *b, uint64_t n) {
#if CFX_LIMB_BITS == 64
    cfx_big_sub_sm_eq(b, n);
#else
    /* 32-bit limbs: subtract low part first, then high */
    uint32_t lo = (uint32_t)n;
    uint32_t hi = (uint32_t)(n >> 32);

    if (lo == 0 && hi == 0) return;
    if (b->n == 0) {
        assert(0 && "cfx_big_sub_u64_eq: underflow");
        return;
    }

    /* subtract lo from limb[0] */
    cfx_limb_t borrow = (b->limb[0] < lo) ? 1 : 0;
    b->limb[0] -= lo;

    /* subtract hi + borrow from limb[1] if needed */
    if (hi || borrow) {
        if (b->n < 2) {
            assert(0 && "cfx_big_sub_u64_eq: underflow");
            return;
        }
        uint64_t sub = (uint64_t)hi + borrow;
        borrow = (b->limb[1] < sub) ? 1 : 0;
        b->limb[1] -= (cfx_limb_t)sub;
    }

    /* propagate remaining borrow */
    size_t i = 2;
    while (borrow && i < b->n) {
        cfx_limb_t old = b->limb[i];
        b->limb[i] -= 1;
        borrow = (old == 0);
        ++i;
    }

    assert(!borrow && "cfx_big_sub_u64_eq: underflow");
    cfx_big_trim(b);
#endif
}

void cfx_big_mul_sm_eq(cfx_big_t *b, cfx_limb_t m) {
    if (m == 1) return;
    if (m == 0 || b->n == 0) {
        cfx_big_init(b);
        cfx_big_from_limb(b, 0);
        return;
    }
    _mul_sm_fast(b, m);
}
