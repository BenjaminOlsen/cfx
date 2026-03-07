/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include "cfx/algo.h"

#include <inttypes.h>

cfx_limb_t cfx_big_div_sm_eq(cfx_big_t *b, cfx_limb_t d) {
    cfx_acc_t rem = 0;
    for (size_t i = b->n; i--;) {
        cfx_acc_t cur = ((cfx_acc_t)rem << CFX_LIMB_BITS) | b->limb[i];
        cfx_limb_t q = (cfx_limb_t)(cur / d);
        rem = (cfx_limb_t)(cur % d);
        b->limb[i] = q;
    }
    return rem;
}

/* Divides x (base 2^64) by uint32_t d, returns remainder. */
uint32_t cfx_big_div_sm_u32_eq(cfx_big_t *b, uint32_t d) {
    cfx_limb_t rem = 0;
    for (size_t i = b->n; i--;) {
        cfx_acc_t cur = ((cfx_acc_t)rem << CFX_LIMB_BITS) | b->limb[i];
        cfx_limb_t q = (cfx_limb_t)(cur / d);   /* 128/32 -> 128/64 promoted; OK */
        rem = (cfx_limb_t)(cur % d);
        b->limb[i] = q;
    }
    cfx_big_trim(b);
    return (uint32_t)rem;
}

/* uses Horner's rule to evaluate the polynomial with x = B=2^64 */
/* P(x)=(...((an​x+an−1​)x+an−2​)x+...+a1​)x+a0 as a running sum:
   acc0 = (an*B + an-1)
   acc1 = acc0*B + an-2 = (an*B + an-1)B + an-2 = anB^2 + an-1B + an-2
   acc2 = acc1*B + an-3 .. etc

   and each step, take % m because the remainder of the sum div m
   is the sum of the remainders modulo m
 */
cfx_limb_t cfx_big_mod_sm(const cfx_big_t *b, cfx_limb_t m) {
    if (b->n == 0) return 0;
    cfx_acc_t acc = 0;
    for (size_t i = b->n; i--;) {
        /* acc << 64 is acc * 2^64, which is (a_n * B); limb[i] is a_(n-1)*/
        acc = ((acc << CFX_LIMB_BITS) + b->limb[i]) % m;
    }
    return (cfx_limb_t)acc;
}

int cfx_big_mulmod(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b, const cfx_big_t *m) {
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_mul(&tmp, a, b);
    int ret = cfx_big_mod(out, &tmp, m);
    cfx_big_free(&tmp);
    return ret;
}

/* Core: q = u / v; r = u % v; any of q or r may be NULL. Returns 0, or -1 if v==0. */
int cfx_big_divrem(cfx_big_t *q, cfx_big_t *r,
    const cfx_big_t *u, const cfx_big_t *v) {
    #if CFX_LIMB_BITS==64
    const cfx_limb_t B_hi_bit = 1ull << 63;
    #elif CFX_LIMB_BITS==32
    const cfx_limb_t B_hi_bit = 1ul << 31;
    #endif

    /* ---- Trivial cases */
    if (cfx_big_is_zero(v)) return -1;

    if (cfx_big_is_zero(u)) {
        if (q) cfx_big_from_limb(q, 0);
        if (r) cfx_big_from_limb(r, 0);
        return 0;
    }

    if (cfx_big_cmp(u, v) < 0) {
        if (q) cfx_big_from_limb(q, 0);
        if (r) cfx_big_copy(r, u);
        return 0;
    }

    /* ---- single-limb fast path */
    if (v->n == 1) {
        /* PRINT_DBG("[FAST] single-limb divisor\n"); */
        cfx_limb_t div = v->limb[0];
        cfx_acc_t rem = 0;
        if (q) {
            cfx_big_reserve(q, u->n); q->n = u->n;
        }
        for (size_t i = u->n; i--;) {
            cfx_acc_t cur = (rem << CFX_LIMB_BITS) | u->limb[i];
            cfx_limb_t qi = (cfx_limb_t)(cur / div);
            rem = cur % div;
            if (q) q->limb[i] = qi;
            /* PRINT_DBG("  [FAST] i=%zd  cur=(" CFX_PRI0xLIMB "|" CFX_PRI0xLIMB ")  qi=" CFX_PRI0xLIMB "  rem=" CFX_PRI0xLIMB "\n", */
            /*        i, (cfx_limb_t)(cur >> CFX_LIMB_BITS), (cfx_limb_t)cur, */
            /*        qi, (cfx_limb_t)rem); */
        }
        if (q) cfx_big_trim(q);
        if (r) cfx_big_from_limb(r, (cfx_limb_t)rem);
        return 0;
    }

    /* ---- Algorithm D (Knuth) */
    /* D1. Normalize so MSB of V's top limb is set. */
    unsigned s = cfx_clz(v->limb[v->n - 1]);
    cfx_big_t V, U;
    cfx_big_init(&V);
    cfx_big_init(&U);
    cfx_big_shl_bits(&V, v, s);    /* V = v << s */
    cfx_big_shl_bits(&U, u, s);    /* U = u << s */

    PRINT_DBG("---- Normalize: s = %u (shift left by s bits)\n", s);
#ifdef CFX_PRINT_DEBUG
    PRINT_BIG("u (orig)     ", u);
    PRINT_BIG("v (orig)     ", v);
    PRINT_BIG("U = u << s   ", &U);
    PRINT_BIG("V = v << s   ", &V);
#endif
    /* Ensure the invariant V[n-1] >= B/2 */
    if (!(V.limb[V.n - 1] & B_hi_bit)) {
        PRINT_DBG("  [WARN] V not normalized: V[n-1] MSB not set!\n");
    }

    /* Append one extra zero limb to U (Knuth D1 convenience) */
    cfx_big_reserve(&U, U.n + 1);
    U.limb[U.n] = 0;
    U.n += 1;

    size_t n = V.n;          /* divisor length after normalization */
    size_t m_plus_n = U.n;   /* dividend length incl the extra top limb */
    size_t m = (m_plus_n >= n) ? (m_plus_n - n) : 0;  /* number of quotient digits (qlen) */

    /* PRINT_DBG("---- Sizes: n(divisor limbs)=%zu, U.n=%zu, qlen=%zu\n", n, U.n, m); */

    cfx_big_t Q;
    cfx_big_init(&Q);
    cfx_big_reserve(&Q, m);
    Q.n = m;
    for (size_t i = 0; i < m; ++i) Q.limb[i] = 0;

    const cfx_limb_t Vn1 = V.limb[n - 1];
    const cfx_limb_t Vn2 = V.limb[n - 2]; /* safe since n>=2 here */

    /* D2..D7 main loop: j = m-1 down to 0 */
    for (size_t j = m; j--;) {
        /* D3: compute qhat, rhat from top two limbs of U and top limb of V */
        cfx_limb_t Ujm  = U.limb[j + n];       /* U_{j+n} */
        cfx_limb_t Ujm1 = U.limb[j + n - 1];   /* U_{j+n-1} */

        cfx_limb_t qhat;
        cfx_acc_t rhat;  /* Use wider type to detect overflow */
        int skip_d3_refine = 0;

        if (Ujm == Vn1) {
            /* Ujm == Vn1 means the trial quotient would be >= b, so clamp to b-1.
             * According to Knuth, rhat = (Ujm*b + Ujm1) - qhat*Vn1
             *                         = (Vn1*b + Ujm1) - (b-1)*Vn1
             *                         = Ujm1 + Vn1
             * If this overflows (>= b), then D3 refinement is skipped entirely
             * since rhat >= b means the test qhat*Vn2 > rhat*b + Ujm2 is always false. */
            qhat = CFX_LIMB_MAX;               /* b-1 */
            rhat = (cfx_acc_t)Ujm1 + Vn1;      /* may be >= b */
            if (rhat >= ((cfx_acc_t)1 << CFX_LIMB_BITS)) {
                skip_d3_refine = 1;            /* rhat >= b, skip refinement */
            }
            /* PRINT_DBG("j=%zd  qhat CLAMP (Ujm==Vn1)  qhat=" CFX_PRI0xLIMB " rhat=0x%" PRIx64 " skip=%d\n", */
            /*        j, qhat, (uint64_t)rhat, skip_d3_refine); */
        } else {
            cfx_acc_t top = ((cfx_acc_t)Ujm << CFX_LIMB_BITS) | Ujm1;
            qhat = (cfx_limb_t)(top / Vn1);
            rhat = top % Vn1;
            /* PRINT_DBG("j=%zd  top=[" CFX_PRI0xLIMB "|" CFX_PRI0xLIMB "]/" CFX_PRI0xLIMB " -> qhat=" CFX_PRI0xLIMB " rhat=" CFX_PRI0xLIMB "\n", */
            /*        j, Ujm, Ujm1, */
            /*        Vn1, qhat, (cfx_limb_t)rhat); */
        }

        /* D3 refine: while qhat*V[n-2] > rhat*b + U[j+n-2], decrement qhat.
         * Only do this if rhat < b (i.e., fits in a limb). */
        if (!skip_d3_refine) {
            cfx_limb_t Ujm2 = U.limb[j + n - 2];
            cfx_acc_t lhs = (cfx_acc_t)qhat * Vn2;
            cfx_acc_t rhs = (rhat << CFX_LIMB_BITS) | Ujm2;

            /* First possible decrement */
            if (lhs > rhs) {
                /* PRINT_DBG("  D3: adjust qhat (lhs>rhs): qhat--, rhat+=Vn1\n"); */
                qhat--;
                rhat += Vn1; /* may now be >= b */
                /* Second possible decrement only if rhat < b */
                if (rhat < ((cfx_acc_t)1 << CFX_LIMB_BITS)) {
                    lhs = (cfx_acc_t)qhat * Vn2;
                    rhs = (rhat << CFX_LIMB_BITS) | Ujm2;
                    if (lhs > rhs) {
                        /* PRINT_DBG("  D3: second adjust qhat: qhat--, rhat+=Vn1\n"); */
                        qhat--;
                    }
                }
            }
        }

        /* D4: Multiply-subtract qhat * V from U window [j .. j+n] */
        cfx_acc_t carry = 0;
        for (size_t i = 0; i < n; ++i) {
            cfx_acc_t t = (cfx_acc_t)qhat * V.limb[i] + carry;
            cfx_limb_t t_lo = (cfx_limb_t)t;
            cfx_limb_t t_hi = (cfx_limb_t)(t >> CFX_LIMB_BITS);

            cfx_limb_t uji = U.limb[j + i];
            cfx_limb_t uji_new = uji - t_lo;
            cfx_limb_t borrow = (uji_new > uji);   /* 1 if underflow */

            U.limb[j + i] = uji_new;
            carry = (cfx_acc_t)t_hi + borrow;

            /* PRINT_DBG("  D4: j=%zd i=%zu  q*V=" CFX_PRI0xLIMB "|" CFX_PRI0xLIMB "  U=" CFX_PRI0xLIMB " -> " CFX_PRI0xLIMB "  carry=" CFX_PRI0xLIMB "(+%u)\n", */
            /*        j, i, */
            /*        t_hi, t_lo, */
            /*        uji, uji_new, */
            /*        (cfx_limb_t)carry, (unsigned)((cfx_limb_t)(carry >> CFX_LIMB_BITS))); */
        }

        /* D6: subtract final carry from U[j+n] */
        cfx_acc_t ujmw = (cfx_acc_t)U.limb[j + n];
        cfx_acc_t diff = ujmw - carry;
        cfx_limb_t ujm_new = (cfx_limb_t)diff;
        cfx_limb_t borrow_out = (ujmw < carry);  /* true iff underflow */

        /* PRINT_DBG("  D6: j=%zd  U[j+n]=" CFX_PRI0xLIMB " - carry=" CFX_PRI0xLIMB "(+%u) -> " CFX_PRI0xLIMB "  borrow_out=%u\n", */
        /*        j, */
        /*        (cfx_limb_t)ujmw, */
        /*        (cfx_limb_t)carry, */
        /*        (unsigned)((cfx_limb_t)(carry >> CFX_LIMB_BITS)), */
        /*        ujm_new, */
        /*        (unsigned)borrow_out); */

        U.limb[j + n] = ujm_new;

        /* D7: if we subtracted too much, add V back and decrement qhat */
        if (borrow_out) {
            /* PRINT_DBG("  D7: add-back (qhat too large): qhat=" CFX_PRI0xLIMB " -> " CFX_PRI0xLIMB "\n", */
            /*        qhat, (qhat - 1)); */
            qhat--;

            cfx_limb_t c = 0;
            for (size_t i = 0; i < n; ++i) {
                cfx_acc_t ssum = (cfx_acc_t)U.limb[j + i] + V.limb[i] + c;
                U.limb[j + i] = (cfx_limb_t)ssum;
                c = (cfx_limb_t)(ssum >> CFX_LIMB_BITS);
                /* PRINT_DBG("    add-back: i=%zu  U+=V+c -> U=" CFX_PRI0xLIMB "  c=" CFX_PRI0xLIMB "\n", */
                /*        i, U.limb[j + i], c); */
            }
            U.limb[j + n] += c;
            /* PRINT_DBG("    add-back: U[j+n]+=c -> " CFX_PRI0xLIMB "\n", */
            /*        U.limb[j + n]); */
        }

        Q.limb[j] = qhat;
        /* PRINT_DBG("  => Q[%zd] = " CFX_PRI0xLIMB "\n", j, qhat); */
    }

    /* D8: Unnormalize remainder: r = (U[0..n-1] >> s) */
    if (r) {
        cfx_big_t Rn;
        cfx_big_init(&Rn);
        cfx_big_reserve(&Rn, n);
        Rn.n = n;
        for (size_t i = 0; i < n; ++i) Rn.limb[i] = U.limb[i];
#ifdef CFX_PRINT_DEBUG
        PRINT_BIG("R (norm pre>>s)", &Rn);
#endif
        cfx_big_shr_bits(r, &Rn, s);
        cfx_big_free(&Rn);
        cfx_big_trim(r);
#ifdef CFX_PRINT_DEBUG
        PRINT_BIG("r (unnormalized)", r);
#endif
    }

    /* Output quotient */
    cfx_big_trim(&Q);
    if (q) {
        cfx_big_swap(&Q, q);
    }
    cfx_big_free(&Q);

    /* Cleanup */
    cfx_big_free(&U);
    cfx_big_free(&V);

#ifdef CFX_PRINT_DEBUG
    if (q) PRINT_BIG("Q (final)", q);
#endif
    return 0;
}

int cfx_big_div(cfx_big_t *q, const cfx_big_t *u, const cfx_big_t *v) {
    return cfx_big_divrem(q, NULL, u, v);
}

int cfx_big_mod(cfx_big_t *r, const cfx_big_t *u, const cfx_big_t *v) {
    return cfx_big_divrem(NULL, r, u, v);
}


/* In-place: u := floor(u/v); optional remainder r. Alias-safe for any combination. */
int cfx_big_divrem_eq(cfx_big_t *u, const cfx_big_t *v, cfx_big_t *r /*nullable*/) {
    cfx_big_t qtmp, rtmp;
    cfx_big_init(&qtmp);
    cfx_big_init(&rtmp);

    int rc = cfx_big_divrem(&qtmp, r ? &rtmp : NULL, u, v);
    if (rc == 0) {
        cfx_big_swap(&qtmp, u);
        if (r) cfx_big_swap(&rtmp, r);
    }

    cfx_big_free(&qtmp);
    cfx_big_free(&rtmp);
    return rc;
}
