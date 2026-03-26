/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include "cfx/algo.h"

/* Montgomery */

/* if x >= n then x -= n */
static inline void cfx_big_cond_sub_n_(cfx_big_t *x, const cfx_big_t *n) {
    if (cfx_big_cmp(x, n) >= 0) cfx_big_sub_eq(x, n);
}

/* n0^{-1} mod 2^64 via Newton; n0 must be odd */
static inline cfx_limb_t inv64_mod2_64_(cfx_limb_t n0) {
    cfx_limb_t x = 1;
    for (int i = 0; i < 6; ++i) x *= (2 - x * n0);
    return x;
}
static inline cfx_limb_t mont_n0inv_(cfx_limb_t n0) {
    return (cfx_limb_t)(0 - inv64_mod2_64_(n0)); /* -n^{-1} mod 2^64 */
}

/* rr = R^2 mod n, with R=2^(64*k). Build rr = 2^(128*k) mod n by repeated doubling. */
static int compute_rr_(cfx_big_t *rr, const cfx_big_t *n, size_t k) {
    cfx_big_from_limb(rr, 1);
    for (size_t bit = 0; bit < 128ull * k; ++bit) {
        /* rr = (rr + rr) mod n (rr < n ⇒ rr+rr < 2n ⇒ one cond. subtract is enough) */
        cfx_big_reserve(rr, rr->n + 1);
        cfx_acc_t carry = 0;
        size_t i = 0;
        for (; i < rr->n; ++i) {
            cfx_acc_t s = (cfx_acc_t)rr->limb[i] * 2 + carry;
            rr->limb[i] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
        }
        if (carry) rr->limb[rr->n++] = (cfx_limb_t)carry;
        cfx_big_cond_sub_n_(rr, n);
    }
    return 1;
}

int cfx_big_mont_ctx_init(cfx_big_mont_ctx_t *ctx, const cfx_big_t *n_in) {
    if (!ctx || !n_in || n_in->n == 0) return 0;
    if ((n_in->limb[0] & 1ull) == 0) return 0; /* modulus must be odd */

    cfx_big_init(&ctx->n);
    cfx_big_init(&ctx->rr);
    cfx_big_init(&ctx->R1);

    cfx_big_assign(&ctx->n, n_in);
    cfx_big_trim(&ctx->n);

    ctx->k = ctx->n.n;
    if (ctx->k == 0) {
        cfx_big_mont_ctx_free(ctx); return 0;
    }

    ctx->n0inv = mont_n0inv_(ctx->n.limb[0]);


    if (!compute_rr_(&ctx->rr, &ctx->n, ctx->k)) {
        cfx_big_mont_ctx_free(ctx);
        return 0;
    }
    /* R1 calculation */
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_from_limb(&one, 1);
    int ok = cfx_big_mont_to(&ctx->R1, &one, ctx); /* R1 = 1*RR/R mod n */
    cfx_big_free(&one);

    if (!ok) return 0;
    return 1;
}

void cfx_big_mont_ctx_free(cfx_big_mont_ctx_t *ctx) {
    if (!ctx) return;
    cfx_big_free(&ctx->n);
    cfx_big_free(&ctx->rr);
    memset(ctx, 0, sizeof(*ctx));
}

int cfx_big_mont_mul(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b, const cfx_big_mont_ctx_t *ctx) {
    if (!out || !a || !b || !ctx) return 0;

    const cfx_big_t *n = &ctx->n;
    const size_t k = ctx->k;

    /* Handle aliasing: if out aliases a or b, we need a temporary.
     * Otherwise, we can work directly in out's buffer. */
    int need_temp = (out == a || out == b || out->limb == a->limb || out->limb == b->limb);

    cfx_big_t T_storage;
    cfx_big_t *T;

    if (need_temp) {
        cfx_big_init(&T_storage);
        cfx_big_reserve(&T_storage, k + 2);
        T = &T_storage;
    } else {
        /* Reuse out's buffer if possible */
        cfx_big_reserve(out, k + 2);
        T = out;
    }

    memset(T->limb, 0, (k + 2) * sizeof(cfx_limb_t));
    T->n = k + 2;

    /* Delegate core CIOS loop to backend implementation */
    cfx_big_mont_mul_impl(T->limb,
        a->limb, a->n,
        b->limb, b->n,
        n->limb,
        ctx->n0inv,
        k);

    /* Final normalization: result is in T[0..k], may be >= n */
    T->n = k + 1;
    cfx_big_trim(T);
    if (cfx_big_cmp(T, n) >= 0) {
        cfx_big_sub_eq(T, n);
    }

    if (need_temp) {
        cfx_big_move(out, &T_storage);
    }
    return 1;
}

/* out = a^e mod n (normal domain). n must be odd. */
int cfx_big_modexp_binary(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *e, const cfx_big_mont_ctx_t *ctx) {
    /* e == 0 -> 1 mod n */
    if (e->n == 0) {
        return cfx_big_mont_from(out, &ctx->R1, ctx); /* MontFrom(R1) == 1 */
    }

    cfx_big_t baseR, accR;
    cfx_big_init(&baseR);
    cfx_big_init(&accR);

    if (!cfx_big_mont_to(&baseR, a, ctx)) goto FAIL;   /* a in Montgomery */
    cfx_big_assign(&accR, &ctx->R1);                   /* acc = 1 (Montgomery) */

    /* msb index of e */
    size_t msb;
    {
        cfx_limb_t top = e->limb[e->n - 1];
        unsigned lz  = cfx_clz(top);
        msb = (e->n * CFX_LIMB_BITS) - 1u - lz;
    }

    /* LTR binary exponentiation */
    for (size_t i = msb + 1; i-- > 0; ) {
        /* square every round */
        if (!cfx_big_mont_mul(&accR, &accR, &accR, ctx)) goto FAIL;

        /* conditional multiply when bit is 1 */
        size_t limb = i / CFX_LIMB_BITS, bit = i % CFX_LIMB_BITS;
        if (limb < e->n && ((e->limb[limb] >> bit) & 1u)) {
            if (!cfx_big_mont_mul(&accR, &accR, &baseR, ctx)) goto FAIL;
        }
    }

    /* back to normal domain */
    if (!cfx_big_mont_from(out, &accR, ctx)) goto FAIL;
    cfx_big_free(&baseR);
    cfx_big_free(&accR);
    return 1;
FAIL:
    cfx_big_free(&baseR);
    cfx_big_free(&accR);
    return 0;
}


/* Convert to Montgomery: aR mod n = MontMul(a, R^2) */
int cfx_big_mont_to(cfx_big_t *out, const cfx_big_t *a, const cfx_big_mont_ctx_t *ctx) {
    /* fast path: if a fits and < n, skip the mod */
    if (a->n <= ctx->k) {
        if (cfx_big_cmp(a, &ctx->n) >= 0) {
            cfx_big_t t; cfx_big_init(&t);
            cfx_big_assign(&t, a);
            /* At most a couple subs when a.n == k */
            do {
                cfx_big_sub_eq(&t, &ctx->n);
            } while (cfx_big_cmp(&t, &ctx->n) >= 0);
            int ok = cfx_big_mont_mul(out, &t, &ctx->rr, ctx);
            cfx_big_free(&t);
            return ok;
        }
        return cfx_big_mont_mul(out, a, &ctx->rr, ctx);
    }

    /*  a is larger than k limbs, take remainder once */
    cfx_big_t rem;
    cfx_big_init(&rem);
    cfx_big_mod(&rem, a, &ctx->n);          /* Knuth D remainder (u mod n) */
    int ok = cfx_big_mont_mul(out, &rem, &ctx->rr, ctx);
    cfx_big_free(&rem);
    return ok;
}

/* Convert from Montgomery: aR * R^{-1} = a mod n = MontMul(aR, 1) */
int cfx_big_mont_from(cfx_big_t *out, const cfx_big_t *aR, const cfx_big_mont_ctx_t *ctx) {
    if (!out || !aR || !ctx) return 0;
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_from_limb(&one, 1);
    int ok = cfx_big_mont_mul(out, aR, &one, ctx);
    cfx_big_free(&one);
    return ok;
}

/* One-liners that hide the ctx */
/* these are pretty useless except for testing */
int cfx_big_mul_mod(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b, const cfx_big_t *n) {
    cfx_big_mont_ctx_t C;
    if (!cfx_big_mont_ctx_init(&C, n)) return 0;
    cfx_big_t aR, bR, r;
    cfx_big_init(&aR);
    cfx_big_init(&bR);
    cfx_big_init(&r);

    int ok = cfx_big_mont_to(&aR, a, &C)
             && cfx_big_mont_to(&bR, b, &C)
             && cfx_big_mont_mul(&r, &aR, &bR, &C)
             && cfx_big_mont_from(out, &r, &C);
    cfx_big_free(&aR);
    cfx_big_free(&bR);
    cfx_big_free(&r);
    cfx_big_mont_ctx_free(&C);
    return ok;
}

int cfx_big_sqr_mod(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *n) {
    cfx_big_mont_ctx_t C;
    if (!cfx_big_mont_ctx_init(&C, n)) return 0;
    cfx_big_t aR, r;
    cfx_big_init(&aR);
    cfx_big_init(&r);
    int ok = cfx_big_mont_to(&aR, a, &C)
             && cfx_big_mont_sqr(&r, &aR, &C)
             && cfx_big_mont_from(out, &r, &C);
    cfx_big_free(&aR);
    cfx_big_free(&r);
    cfx_big_mont_ctx_free(&C);
    return ok;
}

int cfx_big_modexp(cfx_big_t *out, const cfx_big_t *base, const cfx_big_t *exp, const cfx_big_t *n) {
    cfx_big_mont_ctx_t C;
    if (!cfx_big_mont_ctx_init(&C, n)) return 0;

    cfx_big_t A, X;
    cfx_big_init(&A);
    cfx_big_init(&X);

    int ok = cfx_big_mont_to(&A, base, &C);
    if (ok) {
        cfx_big_t one;
        cfx_big_init(&one);
        cfx_big_from_limb(&one, 1);
        ok = cfx_big_mont_to(&X, &one, &C);
        cfx_big_free(&one);
    }
    if (ok) {
        for (size_t i = 0; i < exp->n; ++i) {
            cfx_limb_t w = exp->limb[i];
            for (int b = 0; b < 64; ++b) {
                if (w & 1u) {
                    ok = cfx_big_mont_mul(&X, &X, &A, &C);
                    if (!ok) break;
                }
                w >>= 1;
                ok = cfx_big_mont_mul(&A, &A, &A, &C);
                if (!ok) break; /* square */
            }
            if (!ok) break;
        }
        if (ok) ok = cfx_big_mont_from(out, &X, &C);
    }

    cfx_big_free(&A);
    cfx_big_free(&X);
    cfx_big_mont_ctx_free(&C);
    return ok;
}
