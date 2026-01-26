/*
 * ecm.c - Elliptic Curve Method for factorization
 *
 * Uses Montgomery curves By^2 = x^3 + Ax^2 + x (mod n)
 * which lets us do point ops with just the x-coordinate.
 * Projective coords (X:Z) where x = X/Z avoids divisions.
 */

#include "cfx/ecm.h"
#include "cfx/algo.h"
#include "cfx/rand.h"
#include <string.h>
#include <stdlib.h>

/* store a24 = (A+2)/4 instead of A, all arithmetic in montgomery form */

void cfx_ecm_point_init(cfx_ecm_point_t *P) {
    cfx_big_init(&P->X);
    cfx_big_init(&P->Z);
}

void cfx_ecm_point_free(cfx_ecm_point_t *P) {
    cfx_big_free(&P->X);
    cfx_big_free(&P->Z);
}

void cfx_ecm_point_copy(cfx_ecm_point_t *dst, const cfx_ecm_point_t *src) {
    cfx_big_copy(&dst->X, &src->X);
    cfx_big_copy(&dst->Z, &src->Z);
}

/*
 * point doubling R = 2P
 *
 * u = (X + Z)^2,  v = (X - Z)^2
 * X2 = u * v
 * Z2 = (u - v) * (v + a24*(u - v))
 *
 * the (u-v) term equals 4XZ which is why it works out
 */
static void ecm_point_double(cfx_ecm_point_t *R, const cfx_ecm_point_t *P,
    const cfx_big_t *a24, const cfx_big_mont_ctx_t *ctx) {
    cfx_big_t u, v, diff, t1;
    cfx_big_init(&u);
    cfx_big_init(&v);
    cfx_big_init(&diff);
    cfx_big_init(&t1);

    /* u = (X + Z)^2 */
    cfx_big_copy(&u, &P->X);
    cfx_big_add_eq(&u, &P->Z);
    if (cfx_big_cmp(&u, &ctx->n) >= 0) cfx_big_sub_eq(&u, &ctx->n);
    cfx_big_mont_sqr(&u, &u, ctx);

    /* v = (X - Z)^2, handle underflow */
    cfx_big_copy(&v, &P->X);
    if (cfx_big_cmp(&v, &P->Z) >= 0) {
        cfx_big_sub_eq(&v, &P->Z);
    } else {
        cfx_big_copy(&t1, &P->Z);
        cfx_big_sub_eq(&t1, &P->X);
        cfx_big_copy(&v, &ctx->n);
        cfx_big_sub_eq(&v, &t1);
    }
    cfx_big_mont_sqr(&v, &v, ctx);

    cfx_big_mont_mul(&R->X, &u, &v, ctx);  /* X2 = u * v */

    /* diff = u - v = 4XZ */
    if (cfx_big_cmp(&u, &v) >= 0) {
        cfx_big_copy(&diff, &u);
        cfx_big_sub_eq(&diff, &v);
    } else {
        cfx_big_copy(&diff, &v);
        cfx_big_sub_eq(&diff, &u);
        cfx_big_copy(&t1, &ctx->n);
        cfx_big_sub_eq(&t1, &diff);
        cfx_big_copy(&diff, &t1);
    }

    /* Z2 = diff * (v + a24*diff) */
    cfx_big_mont_mul(&t1, a24, &diff, ctx);
    cfx_big_add_eq(&t1, &v);
    if (cfx_big_cmp(&t1, &ctx->n) >= 0) cfx_big_sub_eq(&t1, &ctx->n);
    cfx_big_mont_mul(&R->Z, &diff, &t1, ctx);

    cfx_big_free(&u);
    cfx_big_free(&v);
    cfx_big_free(&diff);
    cfx_big_free(&t1);
}

/*
 * differential addition: R = P + Q given we know P-Q
 *
 * montgomery's trick - you need the difference to compute the sum
 * when you only have x-coords. the ladder keeps this invariant.
 *
 *   u = (X1 - Z1)(X2 + Z2)
 *   v = (X1 + Z1)(X2 - Z2)
 *   X_out = Z_{P-Q} * (u + v)^2
 *   Z_out = X_{P-Q} * (u - v)^2
 */
static void ecm_point_add(cfx_ecm_point_t *R, const cfx_ecm_point_t *P, const cfx_ecm_point_t *Q,
    const cfx_ecm_point_t *PminusQ, const cfx_big_mont_ctx_t *ctx) {
    cfx_big_t u, v, t1, t2, t3, t4;
    cfx_big_init(&u);
    cfx_big_init(&v);
    cfx_big_init(&t1);
    cfx_big_init(&t2);
    cfx_big_init(&t3);
    cfx_big_init(&t4);

    /* t1 = X1 - Z1 */
    if (cfx_big_cmp(&P->X, &P->Z) >= 0) {
        cfx_big_copy(&t1, &P->X);
        cfx_big_sub_eq(&t1, &P->Z);
    } else {
        cfx_big_copy(&t1, &ctx->n);
        cfx_big_copy(&t3, &P->Z);
        cfx_big_sub_eq(&t3, &P->X);
        cfx_big_sub_eq(&t1, &t3);
    }

    /* t2 = X2 + Z2 */
    cfx_big_copy(&t2, &Q->X);
    cfx_big_add_eq(&t2, &Q->Z);
    if (cfx_big_cmp(&t2, &ctx->n) >= 0) cfx_big_sub_eq(&t2, &ctx->n);

    cfx_big_mont_mul(&u, &t1, &t2, ctx);  /* u = (X1-Z1)(X2+Z2) */

    /* t3 = X1 + Z1 */
    cfx_big_copy(&t3, &P->X);
    cfx_big_add_eq(&t3, &P->Z);
    if (cfx_big_cmp(&t3, &ctx->n) >= 0) cfx_big_sub_eq(&t3, &ctx->n);

    /* t4 = X2 - Z2 */
    if (cfx_big_cmp(&Q->X, &Q->Z) >= 0) {
        cfx_big_copy(&t4, &Q->X);
        cfx_big_sub_eq(&t4, &Q->Z);
    } else {
        cfx_big_copy(&t4, &ctx->n);
        cfx_big_copy(&t1, &Q->Z);
        cfx_big_sub_eq(&t1, &Q->X);
        cfx_big_sub_eq(&t4, &t1);
    }

    cfx_big_mont_mul(&v, &t3, &t4, ctx);  /* v = (X1+Z1)(X2-Z2) */

    /* now compute (u+v)^2 and (u-v)^2 */
    cfx_big_copy(&t1, &u);
    cfx_big_add_eq(&t1, &v);
    if (cfx_big_cmp(&t1, &ctx->n) >= 0) cfx_big_sub_eq(&t1, &ctx->n);

    if (cfx_big_cmp(&u, &v) >= 0) {
        cfx_big_copy(&t2, &u);
        cfx_big_sub_eq(&t2, &v);
    } else {
        cfx_big_copy(&t2, &ctx->n);
        cfx_big_copy(&t3, &v);
        cfx_big_sub_eq(&t3, &u);
        cfx_big_sub_eq(&t2, &t3);
    }

    cfx_big_mont_sqr(&t1, &t1, ctx);
    cfx_big_mont_sqr(&t2, &t2, ctx);

    cfx_big_mont_mul(&R->X, &PminusQ->Z, &t1, ctx);
    cfx_big_mont_mul(&R->Z, &PminusQ->X, &t2, ctx);

    cfx_big_free(&u);
    cfx_big_free(&v);
    cfx_big_free(&t1);
    cfx_big_free(&t2);
    cfx_big_free(&t3);
    cfx_big_free(&t4);
}

/*
 * montgomery ladder for scalar mult R = k * P
 *
 * keep R0 and R1 where R1 - R0 = P always. then we can use
 * differential add since we know the difference.
 *
 * for each bit of k (high to low):
 *   bit=0: R1 = R0+R1, R0 = 2*R0
 *   bit=1: R0 = R0+R1, R1 = 2*R1
 *
 * invariant preserved either way, result ends up in R0
 */
static void ecm_scalar_mul(cfx_ecm_point_t *R, const cfx_ecm_point_t *P, const cfx_big_t *k,
    const cfx_big_t *a24, const cfx_big_mont_ctx_t *ctx) {
    if (cfx_big_is_zero(k)) {
        /* point at infinity */
        cfx_big_from_limb(&R->X, 1);
        cfx_big_from_limb(&R->Z, 0);
        return;
    }

    cfx_ecm_point_t R0, R1, tmp;
    cfx_ecm_point_init(&R0);
    cfx_ecm_point_init(&R1);
    cfx_ecm_point_init(&tmp);

    cfx_ecm_point_copy(&R0, P);          /* R0 = P */
    ecm_point_double(&R1, P, a24, ctx);  /* R1 = 2P */

    /* find top bit */
    size_t top_limb = k->n - 1;
    cfx_limb_t top_val = k->limb[top_limb];
    int top_bit = CFX_LIMB_BITS - 1 - cfx_clz(top_val);

    /* process top limb, skip the MSB (already handled by R0=P, R1=2P) */
    for (int i = top_bit - 1; i >= 0; --i) {
        int bit = (top_val >> i) & 1;
        if (bit == 0) {
            ecm_point_add(&tmp, &R0, &R1, P, ctx);
            ecm_point_double(&R0, &R0, a24, ctx);
            cfx_ecm_point_copy(&R1, &tmp);
        } else {
            ecm_point_add(&tmp, &R0, &R1, P, ctx);
            ecm_point_double(&R1, &R1, a24, ctx);
            cfx_ecm_point_copy(&R0, &tmp);
        }
    }

    /* remaining limbs, all bits */
    for (size_t limb_idx = top_limb; limb_idx-- > 0; ) {
        cfx_limb_t val = k->limb[limb_idx];
        for (int i = CFX_LIMB_BITS - 1; i >= 0; --i) {
            int bit = (val >> i) & 1;
            if (bit == 0) {
                ecm_point_add(&tmp, &R0, &R1, P, ctx);
                ecm_point_double(&R0, &R0, a24, ctx);
                cfx_ecm_point_copy(&R1, &tmp);
            } else {
                ecm_point_add(&tmp, &R0, &R1, P, ctx);
                ecm_point_double(&R1, &R1, a24, ctx);
                cfx_ecm_point_copy(&R0, &tmp);
            }
        }
    }

    cfx_ecm_point_copy(R, &R0);

    cfx_ecm_point_free(&R0);
    cfx_ecm_point_free(&R1);
    cfx_ecm_point_free(&tmp);
}

/*
 * stage 1 - multiply point by product of prime powers up to B1
 *
 * basically computing lcm(1..B1) * P. if the curve's group order
 * (mod some prime factor of n) divides this, Z goes to 0 mod that
 * factor but not mod the other, so gcd(Z,n) finds it.
 *
 * we check gcd periodically, not every prime - faster that way
 */
static int ecm_stage1(cfx_big_t *factor,
    cfx_ecm_point_t *Q,
    uint64_t B1,
    const cfx_big_t *a24,
    const cfx_big_mont_ctx_t *ctx,
    const cfx_big_t *n){
    cfx_big_t k, g;
    cfx_big_init(&k);
    cfx_big_init(&g);

    cfx_ecm_point_t tmp;
    cfx_ecm_point_init(&tmp);

    uint64_t p = 2;
    while (p <= B1) {
        /* find largest p^e <= B1 */
        uint64_t pe = p;
        while (pe <= B1 / p) pe *= p;

        /* Q = pe * Q */
        cfx_big_from_u64(&k, pe);
        ecm_scalar_mul(&tmp, Q, &k, a24, ctx);
        cfx_ecm_point_copy(Q, &tmp);

        /* periodic gcd - every ~100 primes or near the end */
        if (p % 100 == 1 || p > B1 - 100) {
            cfx_big_gcd(&g, &Q->Z, n);
            if (!cfx_big_is_one(&g) && cfx_big_cmp(&g, n) != 0) {
                cfx_big_copy(factor, &g);
                cfx_big_free(&k);
                cfx_big_free(&g);
                cfx_ecm_point_free(&tmp);
                return 1;  /* found one! */
            }
        }

        /* next prime - simple but not fast for large B1 */
        if (p == 2) {
            p = 3;
        } else {
            do { p += 2; } while (p <= B1 && !cfx_is_prime_u64(p));
        }
    }

    /* one last check */
    cfx_big_gcd(&g, &Q->Z, n);
    if (!cfx_big_is_one(&g) && cfx_big_cmp(&g, n) != 0) {
        cfx_big_copy(factor, &g);
        cfx_big_free(&k);
        cfx_big_free(&g);
        cfx_ecm_point_free(&tmp);
        return 1;
    }

    cfx_big_free(&k);
    cfx_big_free(&g);
    cfx_ecm_point_free(&tmp);
    return 0;
}

/*
 * suyama parametrization - construct curve with group order divisible by 12
 *
 * given random sigma, we get:
 *   u = sigma^2 - 5
 *   v = 4*sigma
 *   x0 = u^3, z0 = v^3  (starting point)
 *   a24 = (v-u)^3 * (3u+v) / (16 * u^3 * v)
 *
 * returns 1 on success, 0 if we found a factor (stored in *factor), -1 for bad sigma
 */
static int ecm_suyama_curve(cfx_big_t *a24,
    cfx_ecm_point_t *P,
    cfx_big_t *factor,
    const cfx_big_mont_ctx_t *ctx,
    const cfx_big_t *n,
    uint64_t sigma){
    cfx_big_t u, v, t1, t2, t3, num, denom, g;
    cfx_big_init(&u);
    cfx_big_init(&v);
    cfx_big_init(&t1);
    cfx_big_init(&t2);
    cfx_big_init(&t3);
    cfx_big_init(&num);
    cfx_big_init(&denom);
    cfx_big_init(&g);

    /* u = sigma^2 - 5 */
    cfx_big_from_u64(&u, sigma);
    cfx_big_mul(&u, &u, &u);
    if (cfx_big_cmp_u64(&u, 5) >= 0) {
        cfx_big_sub_u64_eq(&u, 5);
    } else {
        /* sigma^2 < 5, add n and subtract */
        cfx_big_add(&u, &u, n);
        cfx_big_sub_u64_eq(&u, 5);
    }
    cfx_big_mod(&u, &u, n);

    /* v = 4*sigma */
    cfx_big_from_u64(&v, sigma);
    cfx_big_shl_bits_eq(&v, 2);
    cfx_big_mod(&v, &v, n);

    /* x0 = u^3 */
    cfx_big_mul(&t1, &u, &u);
    cfx_big_mod(&t1, &t1, n);
    cfx_big_mul(&P->X, &t1, &u);
    cfx_big_mod(&P->X, &P->X, n);

    /* z0 = v^3 */
    cfx_big_mul(&t1, &v, &v);
    cfx_big_mod(&t1, &t1, n);
    cfx_big_mul(&P->Z, &t1, &v);
    cfx_big_mod(&P->Z, &P->Z, n);

    /* (v - u) -- careful with underflow */
    if (cfx_big_cmp(&v, &u) >= 0) {
        cfx_big_copy(&t1, &v);
        cfx_big_sub_eq(&t1, &u);
    } else {
        cfx_big_copy(&t1, n);
        cfx_big_copy(&t2, &u);
        cfx_big_sub_eq(&t2, &v);
        cfx_big_sub_eq(&t1, &t2);
    }

    /* (v - u)^3 */
    cfx_big_mul(&t2, &t1, &t1);
    cfx_big_mod(&t2, &t2, n);
    cfx_big_mul(&t1, &t2, &t1);
    cfx_big_mod(&t1, &t1, n);  /* t1 = (v-u)^3 */

    /* 3u + v */
    cfx_big_copy(&t2, &u);
    cfx_big_add_eq(&t2, &u);
    cfx_big_add_eq(&t2, &u);
    cfx_big_add_eq(&t2, &v);
    cfx_big_mod(&t2, &t2, n);  /* t2 = 3u + v */

    /* numerator = (v-u)^3 * (3u+v) */
    cfx_big_mul(&num, &t1, &t2);
    cfx_big_mod(&num, &num, n);

    /* denominator = 16 * u^3 * v */
    cfx_big_mul(&t1, &u, &u);
    cfx_big_mod(&t1, &t1, n);
    cfx_big_mul(&t1, &t1, &u);
    cfx_big_mod(&t1, &t1, n);  /* t1 = u^3 */
    cfx_big_mul(&denom, &t1, &v);
    cfx_big_mod(&denom, &denom, n);
    cfx_big_shl_bits_eq(&denom, 4);  /* *16 */
    cfx_big_mod(&denom, &denom, n);

    /* check if we can invert denom */
    cfx_big_gcd(&g, &denom, n);
    if (!cfx_big_is_one(&g)) {
        if (cfx_big_cmp(&g, n) != 0) {
            /* found a factor! */
            cfx_big_copy(factor, &g);
            cfx_big_free(&u); cfx_big_free(&v);
            cfx_big_free(&t1); cfx_big_free(&t2); cfx_big_free(&t3);
            cfx_big_free(&num); cfx_big_free(&denom); cfx_big_free(&g);
            return 0;
        }
        /* g == n means bad sigma, denom is 0 mod n */
        cfx_big_free(&u); cfx_big_free(&v);
        cfx_big_free(&t1); cfx_big_free(&t2); cfx_big_free(&t3);
        cfx_big_free(&num); cfx_big_free(&denom); cfx_big_free(&g);
        return -1;
    }

    /* a24 = num / denom = num * denom^(-1) */
    cfx_big_t inv;
    cfx_big_init(&inv);
    cfx_big_modinv(&inv, &denom, n);
    cfx_big_mul(a24, &num, &inv);
    cfx_big_mod(a24, a24, n);
    cfx_big_free(&inv);

    /* convert to montgomery form */
    cfx_big_mont_to(a24, a24, ctx);
    cfx_big_mont_to(&P->X, &P->X, ctx);
    cfx_big_mont_to(&P->Z, &P->Z, ctx);

    cfx_big_free(&u); cfx_big_free(&v);
    cfx_big_free(&t1); cfx_big_free(&t2); cfx_big_free(&t3);
    cfx_big_free(&num); cfx_big_free(&denom); cfx_big_free(&g);
    return 1;
}

/*
 * main ecm routine - try several curves hoping one has smooth order
 */
int cfx_ecm_factor(cfx_big_t *factor, const cfx_big_t *n,
    uint64_t B1, unsigned curves){

    if (cfx_big_is_zero(n) || cfx_big_is_one(n)) {
        return 0;
    }
    if (cfx_big_is_even(n)) {
        cfx_big_from_limb(factor, 2);
        return 1;
    }

    cfx_big_mont_ctx_t ctx;
    if (!cfx_big_mont_ctx_init(&ctx, n)) {
        return 0;  /* n not odd? shouldn't happen */
    }

    cfx_big_t a24;
    cfx_big_init(&a24);

    cfx_ecm_point_t P, Q;
    cfx_ecm_point_init(&P);
    cfx_ecm_point_init(&Q);

    int found = 0;
    for (unsigned curve = 0; curve < curves && !found; ++curve) {
        /* deterministic sigma for reproducibility - avoid small values */
        uint64_t sigma = 6 + curve * 271828182ULL;

        int rc = ecm_suyama_curve(&a24, &P, factor, &ctx, n, sigma);
        if (rc == 0) {
            /* suyama found a factor during curve setup! */
            found = 1;
            break;
        }
        if (rc < 0) {
            /* bad sigma, skip to next */
            continue;
        }

        cfx_ecm_point_copy(&Q, &P);
        found = ecm_stage1(factor, &Q, B1, &a24, &ctx, n);
    }

    cfx_big_free(&a24);
    cfx_ecm_point_free(&P);
    cfx_ecm_point_free(&Q);
    cfx_big_mont_ctx_free(&ctx);

    return found;
}

/*
 * auto-tuned version - picks B1 and curve count based on size of n
 *
 * these are rough heuristics, nothing too scientific. bigger numbers
 * need more curves and higher B1 to find factors with reasonable prob.
 */
int cfx_ecm_factor_auto(cfx_big_t *factor, const cfx_big_t *n){
    if (n->n == 0) return 0;
    size_t bits = (n->n - 1) * CFX_LIMB_BITS + (CFX_LIMB_BITS - cfx_clz(n->limb[n->n - 1]));

    uint64_t B1;
    unsigned curves;

    if (bits <= 64) {
        B1 = 2000;
        curves = 10;
    } else if (bits <= 96) {
        B1 = 50000;
        curves = 25;
    } else if (bits <= 128) {
        B1 = 1000000;
        curves = 50;
    } else if (bits <= 192) {
        B1 = 10000000;
        curves = 100;
    } else {
        B1 = 50000000;
        curves = 200;
    }

    return cfx_ecm_factor(factor, n, B1, curves);
}
