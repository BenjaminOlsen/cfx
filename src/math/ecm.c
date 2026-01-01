/*
 * ecm.c - Elliptic Curve Method (ECM) implementation
 *
 * We use Montgomery curves: By^2 = x³ + Ax^2 + x (mod n)
 *
 * The beauty of Montgomery curves is that point operations can be done
 * using only the x-coordinate! This is because of the "Montgomery ladder"
 * which computes both P+Q and P-Q simultaneously.
 *
 * We use projective coordinates (X:Z) where x = X/Z to avoid divisions.
 */

#include "cfx/ecm.h"
#include "cfx/algo.h"
#include <string.h>
#include <stdlib.h>

/* ============================================================================
 * MONTGOMERY CURVE ARITHMETIC
 *
 * For a Montgomery curve By^2 = x³ + Ax^2 + x, we need the constant:
 *   a24 = (A + 2) / 4
 *
 * All arithmetic is done in Montgomery form (using cfx_big_mont_* functions)
 * for efficient modular multiplication.
 * ============================================================================
 */

void cfx_ecm_point_init(cfx_ecm_point_t* P) {
    cfx_big_init(&P->X);
    cfx_big_init(&P->Z);
}

void cfx_ecm_point_free(cfx_ecm_point_t* P) {
    cfx_big_free(&P->X);
    cfx_big_free(&P->Z);
}

void cfx_ecm_point_copy(cfx_ecm_point_t* dst, const cfx_ecm_point_t* src) {
    cfx_big_copy(&dst->X, &src->X);
    cfx_big_copy(&dst->Z, &src->Z);
}

/*
 * Point Doubling on Montgomery Curve
 * ===================================
 *
 * Given point P = (X1:Z1), compute 2P = (X2:Z2)
 *
 * The formulas (in projective coordinates):
 *
 *   u = (X1 + Z1)^2
 *   v = (X1 - Z1)^2
 *   diff = u - v                    [this equals 4*X1*Z1]
 *   X2 = u * v
 *   Z2 = diff * (v + a24*diff)
 *
 * where a24 = (A+2)/4 for the curve parameter A.
 *
 * Cost: 3 multiplications, 2 squarings, 4 add/sub
 *
 * Why does this work? On a Montgomery curve, if (x,y) is a point, then:
 *   x(2P) = (x^2 - 1)^2 / (4x(x^2 + Ax + 1))
 *
 * The projective formulas avoid the division by working with (X:Z).
 */
static void ecm_point_double(cfx_ecm_point_t* R, const cfx_ecm_point_t* P,
                             const cfx_big_t* a24, const cfx_big_mont_ctx_t* ctx) {
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

    /* v = (X - Z)^2 */
    cfx_big_copy(&v, &P->X);
    if (cfx_big_cmp(&v, &P->Z) >= 0) {
        cfx_big_sub_eq(&v, &P->Z);
    } else {
        cfx_big_copy(&t1, &P->Z);
        cfx_big_sub_eq(&t1, &P->X);
        cfx_big_copy(&v, &ctx->n);
        cfx_big_sub_eq(&v, &t1);  /* v = n - (Z - X) = X - Z (mod n) */
    }
    cfx_big_mont_sqr(&v, &v, ctx);

    /* X2 = u * v */
    cfx_big_mont_mul(&R->X, &u, &v, ctx);

    /* diff = u - v = (X+Z)^2 - (X-Z)^2 = 4XZ */
    if (cfx_big_cmp(&u, &v) >= 0) {
        cfx_big_copy(&diff, &u);
        cfx_big_sub_eq(&diff, &v);
    } else {
        cfx_big_copy(&diff, &v);
        cfx_big_sub_eq(&diff, &u);
        cfx_big_copy(&t1, &ctx->n);
        cfx_big_sub_eq(&t1, &diff);
        cfx_big_copy(&diff, &t1);  /* diff = n - (v-u) */
    }

    /* t1 = a24 * diff */
    cfx_big_mont_mul(&t1, a24, &diff, ctx);

    /* t1 = v + a24*diff */
    cfx_big_add_eq(&t1, &v);
    if (cfx_big_cmp(&t1, &ctx->n) >= 0) cfx_big_sub_eq(&t1, &ctx->n);

    /* Z2 = diff * (v + a24*diff) */
    cfx_big_mont_mul(&R->Z, &diff, &t1, ctx);

    cfx_big_free(&u);
    cfx_big_free(&v);
    cfx_big_free(&diff);
    cfx_big_free(&t1);
}

/*
 * Differential Point Addition on Montgomery Curve
 * ================================================
 *
 * Given P = (X1:Z1), Q = (X2:Z2), and P-Q = (X4:Z4),
 * compute P+Q = (X4:Z4)
 *
 * Note: We need to know P-Q! This is the "differential" part.
 * The Montgomery ladder maintains both P and P-Q at all times.
 *
 * Formulas:
 *   u = (X1 - Z1)(X2 + Z2)
 *   v = (X1 + Z1)(X2 - Z2)
 *   X4 = Z4 * (u + v)^2
 *   Z4 = X4 * (u - v)^2
 *
 * Cost: 2 multiplications, 2 squarings, 4 add/sub
 *
 * Why do we need P-Q? Montgomery's insight was that x(P+Q) and x(P-Q)
 * are related by a simple formula involving x(P) and x(Q). If we know
 * any three of {P, Q, P+Q, P-Q}, we can compute the fourth.
 */
static void ecm_point_add(cfx_ecm_point_t* R, const cfx_ecm_point_t* P, const cfx_ecm_point_t* Q,
                          const cfx_ecm_point_t* PminusQ, const cfx_big_mont_ctx_t* ctx) {
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

    /* u = (X1 - Z1)(X2 + Z2) */
    cfx_big_mont_mul(&u, &t1, &t2, ctx);

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

    /* v = (X1 + Z1)(X2 - Z2) */
    cfx_big_mont_mul(&v, &t3, &t4, ctx);

    /* t1 = u + v */
    cfx_big_copy(&t1, &u);
    cfx_big_add_eq(&t1, &v);
    if (cfx_big_cmp(&t1, &ctx->n) >= 0) cfx_big_sub_eq(&t1, &ctx->n);

    /* t2 = u - v */
    if (cfx_big_cmp(&u, &v) >= 0) {
        cfx_big_copy(&t2, &u);
        cfx_big_sub_eq(&t2, &v);
    } else {
        cfx_big_copy(&t2, &ctx->n);
        cfx_big_copy(&t3, &v);
        cfx_big_sub_eq(&t3, &u);
        cfx_big_sub_eq(&t2, &t3);
    }

    /* t1 = (u + v)^2 */
    cfx_big_mont_sqr(&t1, &t1, ctx);

    /* t2 = (u - v)^2 */
    cfx_big_mont_sqr(&t2, &t2, ctx);

    /* X4 = Z4 * (u + v)^2 */
    cfx_big_mont_mul(&R->X, &PminusQ->Z, &t1, ctx);

    /* Z4 = X4 * (u - v)^2 */
    cfx_big_mont_mul(&R->Z, &PminusQ->X, &t2, ctx);

    cfx_big_free(&u);
    cfx_big_free(&v);
    cfx_big_free(&t1);
    cfx_big_free(&t2);
    cfx_big_free(&t3);
    cfx_big_free(&t4);
}

/*
 * Montgomery Ladder for Scalar Multiplication
 * ===========================================
 *
 * Compute k*P for scalar k and point P.
 *
 * The Montgomery ladder is elegant: it maintains two points (R0, R1)
 * such that R1 - R0 = P at all times. This invariant lets us use
 * differential addition!
 *
 * Algorithm (for k with bits b_{n-1}, b_{n-2}, ..., b_1, b_0):
 *   R0 = P
 *   R1 = 2P
 *   for i from n-2 down to 0:
 *       if bit i of k is 0:
 *           R1 = R0 + R1  (diff = P)
 *           R0 = 2*R0
 *       else:
 *           R0 = R0 + R1  (diff = P)
 *           R1 = 2*R1
 *   return R0
 *
 * The invariant is maintained because:
 *   - If we double R0 and add to get R1: new R1 - new R0 = (R0+R1) - 2R0 = R1-R0 = P ✓
 *   - If we double R1 and add to get R0: new R1 - new R0 = 2R1 - (R0+R1) = R1-R0 = P ✓
 */
static void ecm_scalar_mul(cfx_ecm_point_t* R, const cfx_ecm_point_t* P, const cfx_big_t* k,
                           const cfx_big_t* a24, const cfx_big_mont_ctx_t* ctx) {
    if (cfx_big_is_zero(k)) {
        /* k = 0: return point at infinity (Z = 0) */
        cfx_big_from_limb(&R->X, 1);
        cfx_big_from_limb(&R->Z, 0);
        return;
    }

    cfx_ecm_point_t R0, R1, tmp;
    cfx_ecm_point_init(&R0);
    cfx_ecm_point_init(&R1);
    cfx_ecm_point_init(&tmp);

    /* R0 = P */
    cfx_ecm_point_copy(&R0, P);

    /* R1 = 2P */
    ecm_point_double(&R1, P, a24, ctx);

    /* Find the highest bit of k */
    size_t top_limb = k->n - 1;
    cfx_limb_t top_val = k->limb[top_limb];
    int top_bit = CFX_LIMB_BITS - 1 - cfx_clz(top_val);

    /* Process bits from second-highest down to 0 */
    for (int i = top_bit - 1; i >= 0; --i) {
        int bit = (top_val >> i) & 1;

        if (bit == 0) {
            /* R1 = R0 + R1, R0 = 2*R0 */
            ecm_point_add(&tmp, &R0, &R1, P, ctx);
            ecm_point_double(&R0, &R0, a24, ctx);
            cfx_ecm_point_copy(&R1, &tmp);
        } else {
            /* R0 = R0 + R1, R1 = 2*R1 */
            ecm_point_add(&tmp, &R0, &R1, P, ctx);
            ecm_point_double(&R1, &R1, a24, ctx);
            cfx_ecm_point_copy(&R0, &tmp);
        }
    }

    /* Continue with remaining limbs */
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

    /* Result is in R0 */
    cfx_ecm_point_copy(R, &R0);

    cfx_ecm_point_free(&R0);
    cfx_ecm_point_free(&R1);
    cfx_ecm_point_free(&tmp);
}

/*
 * ECM Stage 1
 * ===========
 *
 * Compute Q = k*P where k = ∏(p^e) for all prime powers p^e ≤ B1.
 *
 * This is equivalent to computing lcm(1,2,...,B1)*P.
 *
 * For each prime p ≤ B1:
 *   - Find the largest e such that p^e ≤ B1
 *   - Multiply the point by p^e
 *
 * After stage 1, if gcd(Z, n) is non-trivial, we found a factor!
 *
 * Why this works: If the group order mod p is B1-smooth (all prime
 * factors ≤ B1), then k is a multiple of the group order, so k*P = O,
 * meaning Z becomes 0 mod p (but probably not mod q).
 */
static int ecm_stage1(cfx_big_t* factor,
                      cfx_ecm_point_t* Q,
                      uint64_t B1,
                      const cfx_big_t* a24,
                      const cfx_big_mont_ctx_t* ctx,
                      const cfx_big_t* n)
{
    cfx_big_t k, g;
    cfx_big_init(&k);
    cfx_big_init(&g);

    cfx_ecm_point_t tmp;
    cfx_ecm_point_init(&tmp);

    /* For each prime p up to B1 */
    /* We use the precomputed primes if available, otherwise generate */
    uint64_t p = 2;

    while (p <= B1) {
        /* Find largest e such that p^e <= B1 */
        uint64_t pe = p;
        while (pe <= B1 / p) {
            pe *= p;
        }

        /* Q = pe * Q */
        cfx_big_from_u64(&k, pe);
        ecm_scalar_mul(&tmp, Q, &k, a24, ctx);
        cfx_ecm_point_copy(Q, &tmp);

        /* Check if Z became zero (mod some factor) */
        /* We do this periodically, not every prime, for efficiency */
        if (p % 100 == 1 || p > B1 - 100) {
            cfx_big_gcd(&g, &Q->Z, n);
            if (!cfx_big_is_one(&g) && cfx_big_cmp(&g, n) != 0) {
                /* Found a factor! */
                cfx_big_copy(factor, &g);
                cfx_big_free(&k);
                cfx_big_free(&g);
                cfx_ecm_point_free(&tmp);
                return 1;
            }
        }

        /* Next prime (simple sieve for small primes, then increment by 2) */
        if (p == 2) {
            p = 3;
        } else {
            /* Simple prime finding - not optimal but works */
            do {
                p += 2;
            } while (p <= B1 && !cfx_is_prime_u64(p));
        }
    }

    /* Final GCD check */
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
    return 0;  /* No factor found */
}

/*
 * Generate a random curve and starting point
 * ==========================================
 *
 * We use Suyama's parametrization to generate curves with group order
 * divisible by 12, which improves the probability of finding smooth orders.
 *
 * Given a random σ > 5:
 *   u = σ^2 - 5
 *   v = 4σ
 *   x₀ = u³ / v³  (our starting point's x-coordinate)
 *   A = (v-u)³(3u+v) / (4u³v) - 2
 *
 * The resulting curve has group order divisible by 12.
 *
 * For simplicity, we use a simpler method: pick random A and x₀.
 */
static void ecm_random_curve(cfx_big_t* a24,
                             cfx_ecm_point_t* P,
                             const cfx_big_mont_ctx_t* ctx,
                             uint64_t seed)
{
    /* Simple PRNG for deterministic curve generation */
    uint64_t state = seed * 6364136223846793005ULL + 1442695040888963407ULL;

    /* Generate curve parameter A (we store a24 = (A+2)/4) */
    /* For simplicity, pick a24 directly as a random value */
    cfx_big_from_u64(a24, state);
    state = state * 6364136223846793005ULL + 1442695040888963407ULL;

    /* Reduce mod n and convert to Montgomery form */
    if (cfx_big_cmp(a24, &ctx->n) >= 0) {
        cfx_big_mod(a24, a24, &ctx->n);
    }
    cfx_big_mont_to(a24, a24, ctx);

    /* Generate starting point x-coordinate */
    cfx_big_from_u64(&P->X, state);
    if (cfx_big_cmp(&P->X, &ctx->n) >= 0) {
        cfx_big_mod(&P->X, &P->X, &ctx->n);
    }
    cfx_big_mont_to(&P->X, &P->X, ctx);

    /* Z = 1 (in Montgomery form) */
    cfx_big_copy(&P->Z, &ctx->R1);
}

/*
 * Main ECM factorization
 */
int cfx_ecm_factor(cfx_big_t* factor, const cfx_big_t* n,
                   uint64_t B1, unsigned curves)
{
    /* Basic checks */
    if (cfx_big_is_zero(n) || cfx_big_is_one(n)) {
        return 0;
    }

    /* Check for even n */
    if (cfx_big_is_even(n)) {
        cfx_big_from_limb(factor, 2);
        return 1;
    }

    /* Initialize Montgomery context */
    cfx_big_mont_ctx_t ctx;
    if (!cfx_big_mont_ctx_init(&ctx, n)) {
        return 0;
    }

    cfx_big_t a24;
    cfx_big_init(&a24);

    cfx_ecm_point_t P, Q;
    cfx_ecm_point_init(&P);
    cfx_ecm_point_init(&Q);

    int found = 0;

    for (unsigned curve = 0; curve < curves && !found; ++curve) {
        /* Generate random curve and point */
        uint64_t seed = 314159265ULL + curve * 271828182ULL;
        ecm_random_curve(&a24, &P, &ctx, seed);

        /* Copy P to Q (we'll modify Q in stage 1) */
        cfx_ecm_point_copy(&Q, &P);

        /* Run stage 1 */
        found = ecm_stage1(factor, &Q, B1, &a24, &ctx, n);
    }

    cfx_big_free(&a24);
    cfx_ecm_point_free(&P);
    cfx_ecm_point_free(&Q);
    cfx_big_mont_ctx_free(&ctx);

    return found;
}

/*
 * Auto-tuned ECM
 *
 * Picks B1 and curve count based on n's size.
 * These are rough heuristics based on the expected factor size.
 */
int cfx_ecm_factor_auto(cfx_big_t* factor, const cfx_big_t* n)
{
    size_t bits = n->n * CFX_LIMB_BITS;

    /* Heuristic: assume smallest factor is about bits/2 */
    /* Adjust B1 and curves based on expected factor size */

    uint64_t B1;
    unsigned curves;

    if (bits <= 64) {
        /* Small number - use modest parameters */
        B1 = 2000;
        curves = 10;
    } else if (bits <= 96) {
        /* Up to ~48 bit factors */
        B1 = 50000;
        curves = 25;
    } else if (bits <= 128) {
        /* Up to ~64 bit factors */
        B1 = 1000000;
        curves = 50;
    } else if (bits <= 192) {
        /* Up to ~96 bit factors */
        B1 = 10000000;
        curves = 100;
    } else {
        /* Large numbers - use aggressive parameters */
        B1 = 50000000;
        curves = 200;
    }

    return cfx_ecm_factor(factor, n, B1, curves);
}
