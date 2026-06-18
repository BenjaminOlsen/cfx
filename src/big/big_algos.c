/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include "cfx/sbig.h"
#include "cfx/algo.h"


void cfx_big_exp(cfx_big_t *out, const cfx_big_t *n, const cfx_big_t *p) {
    if (cfx_big_is_zero(p)) {
        cfx_big_from_limb(out, 1); return;
    }
    if (cfx_big_is_zero(n)) {
        PRINT_BIG(">>>>>>>>>> n is zero!", n); cfx_big_from_limb(out, 0); return;
    }
    if (cfx_big_eq_u64(n, 1)) {
        cfx_big_from_limb(out, 1); return;
    }

    cfx_big_t acc, pp, np; /* accumulator, p copy, n^p*/
    cfx_big_init(&acc);
    cfx_big_init(&pp);
    cfx_big_init(&np);

    cfx_big_from_limb(&np, 1);
    cfx_big_copy(&pp, p);
    cfx_big_copy(&acc, n);

    while (!cfx_big_is_zero(&pp)) {
        if (pp.n && (pp.limb[0] & 1)) {
            cfx_big_mul_auto(&np, &acc);
        }
        cfx_big_shr_bits(&pp, &pp, 1);
        if (!cfx_big_is_zero(&pp)) {
            cfx_big_mul_auto(&acc, &acc);
        }
    }
    cfx_big_move(out, &np);
    cfx_big_free(&np);
    cfx_big_free(&pp);
    cfx_big_free(&acc);
}

void cfx_big_exp_u64(cfx_big_t *out, const cfx_big_t *n, cfx_limb_t p) {
    if (p == 0) {
        cfx_big_from_limb(out, 1); return;
    }
    if (cfx_big_is_zero(n)) {
        cfx_big_from_limb(out, 0); return;
    }
    if (cfx_big_eq_u64(n, 1)) {
        cfx_big_from_limb(out, 1); return;
    }

    cfx_big_t acc, np; /* accumulator, p copy, n^p*/
    cfx_big_init(&acc);
    cfx_big_init(&np);
    cfx_big_from_limb(&np, 1);
    cfx_big_copy(&acc, n);
    while (p) {
        if (p & 1) {
            cfx_big_mul_auto(&np, &acc);
        }
        p >>= 1;
        if (p) cfx_big_mul_auto(&acc, &acc);
    }
    cfx_big_move(out, &np);
    cfx_big_free(&np);
    cfx_big_free(&acc);
}

/* if p is secret, this leaks info */
void cfx_big_mod_exp(cfx_big_t *out, const cfx_big_t *n, const cfx_big_t *p, const cfx_big_t *m) {
    if (cfx_big_is_zero(p)) {
        cfx_big_from_limb(out, 1); return;
    }
    if (cfx_big_is_zero(n)) {
        cfx_big_from_limb(out, 0); return;
    }
    if (cfx_big_eq_u64(n, 1)) {
        cfx_big_from_limb(out, 1); return;
    }

    cfx_big_t acc, pp, np; /* accumulator, p copy, n^p*/
    cfx_big_init(&acc);
    cfx_big_init(&pp);
    cfx_big_init(&np);
    cfx_big_from_limb(&np, 1);
    cfx_big_copy(&pp, p);
    cfx_big_copy(&acc, n);
    while (!cfx_big_is_zero(&pp)) {
        if (pp.n && (pp.limb[0] & 1)) {
            cfx_big_mul_auto(&np, &acc);
            cfx_big_mod(&np, &np, m);
        }
        cfx_big_shr_bits(&pp, &pp, 1);
        if (!cfx_big_is_zero(&pp)) cfx_big_mul_auto(&acc, &acc);
    }
    cfx_big_move(out, &np);
    cfx_big_free(&np);
    cfx_big_free(&pp);
    cfx_big_free(&acc);
}


/* out = (a^e) mod m */
static void powmod_u64_base(cfx_big_t *out, cfx_limb_t a, const cfx_big_t *e, const cfx_big_t *mod) {
    /* out = a^e mod mod, with small base a and big exponent e */
    cfx_big_t base, res, ee;
    cfx_big_init(&base);
    cfx_big_init(&res);
    cfx_big_init(&ee);
    cfx_big_assign_sm(&base, a);
    cfx_big_mod(&base, &base, mod);

    cfx_big_assign_sm(&res, 1);
    cfx_big_copy(&ee, e);

    /* square & multiply */
    while (!cfx_big_is_zero(&ee)) {
        if (!cfx_big_is_even(&ee)) {  /* if (ee.limb[0] & 1) */
            cfx_big_mulmod(&res, &res, &base, mod);
        }
        cfx_big_mulmod(&base, &base, &base, mod);
        cfx_big_shr_bits_eq(&ee, 1);
    }
    cfx_big_copy(out, &res);
    cfx_big_free(&ee);
    cfx_big_free(&res);
    cfx_big_free(&base);

}

/* Single Miller–Rabin round.
 * n must be odd > 2.  n-1 = d * 2^s, with d odd.
 * a is a small base (2,3,5,7,...) < 2^64.
 */
static int miller_rabin_once(const cfx_big_t *n, cfx_limb_t a, const cfx_big_t *d, cfx_limb_t s) {
    cfx_big_t x, nm1;
    int rc = 0;

    cfx_big_init(&x);
    cfx_big_init(&nm1);

    /* x = a^d mod n */
    powmod_u64_base(&x, a, d, n);

    /* nm1 = n-1 */
    cfx_big_copy(&nm1, n);
    cfx_big_sub_sm_eq(&nm1, 1);

    if (cfx_big_is_one(&x)) {
        rc = 1;
        goto cleanup;
    }
    if (cfx_big_cmp(&x, &nm1) == 0) {
        rc = 1;
        goto cleanup;
    }

    /* TODO: use montgomery modular exponentiation here */
    for (cfx_limb_t r = 1; r < s; ++r) {
        cfx_big_mulmod(&x, &x, &x, n);
        if (cfx_big_cmp(&x, &nm1) == 0) {
            rc = 1;
            goto cleanup;
        }
    }

    rc = 0;     /* witness: composite */

cleanup:
    cfx_big_free(&x);
    cfx_big_free(&nm1);
    return rc;
}



/* Does miller - rabin's primality test on a big */
int cfx_big_is_prime(const cfx_big_t *n) {

    if (cfx_big_cmp_sm(n, 2) < 0)
        return 0;

    static const uint32_t small_primes[] = {
        2,3,5,7,11,13,17,19,23,29,31,37,
        41,43,47,53,59,61,67,71,73,79,83,89,97,0
    };

    /* Small trial division first */
    for (int i = 0; small_primes[i]; ++i) {
        uint32_t p = small_primes[i];
        int cmp = cfx_big_cmp_sm(n, p);
        if (cmp == 0)
            return 1;
        if (cmp > 0 && cfx_big_mod_sm(n, p) == 0)
            return 0;
    }

    if (cfx_big_is_even(n))
        return 0;

    /* Write n-1 = d * 2^s with d odd */
    cfx_big_t d;
    cfx_big_init(&d);
    cfx_limb_t s = 0;

    cfx_big_copy(&d, n);
    cfx_big_sub_sm_eq(&d, 1);          /* d = n-1 */

    while (!cfx_big_is_zero(&d) && cfx_big_is_even(&d)) {
        cfx_big_shr_bits_eq(&d, 1);
        ++s;
    }

    if (cfx_big_is_zero(&d)) {
        cfx_big_free(&d);
        return 0;
    }

    /* for numbers up to 64 bits, use the "Sinclair 7" bases which are
       deterministic for all 64-bit inputs. For larger numbers, use
       additional small prime bases for higher confidence. */
    size_t bitlen = cfx_big_bitlen(n);

    if (bitlen <= 64) {
        static const uint64_t sinclair7[] = {
            2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL, 0
        };
        for (size_t i = 0; sinclair7[i]; ++i) {
            if (!miller_rabin_once(n, (cfx_limb_t)sinclair7[i], &d, s)) {
                cfx_big_free(&d);
                return 0;
            }
        }
    } else {
        /* For larger numbers, use first 12 primes as witnesses.
           This is probabilistic but with very low false positive rate. */
        static const cfx_limb_t bases[] = {
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 0
        };
        for (size_t i = 0; bases[i]; ++i) {
            if (!miller_rabin_once(n, bases[i], &d, s)) {
                cfx_big_free(&d);
                return 0;
            }
        }
    }

    cfx_big_free(&d);
    return 1;
}

/* Binary GCD algorithm for big integers */
void cfx_big_gcd(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
    if (cfx_big_is_zero(a)) {
        cfx_big_copy(out, b);
        return;
    }
    if (cfx_big_is_zero(b)) {
        cfx_big_copy(out, a);
        return;
    }

    cfx_big_t u, v;
    cfx_big_init(&u);
    cfx_big_init(&v);
    cfx_big_copy(&u, a);
    cfx_big_copy(&v, b);

    /* cnt common factors of 2 */
    size_t shift = 0;
    while (cfx_big_is_even(&u) && cfx_big_is_even(&v)) {
        cfx_big_shr_bits_eq(&u, 1);
        cfx_big_shr_bits_eq(&v, 1);
        shift++;
    }

    /* remove remaining factors of 2 from u */
    while (cfx_big_is_even(&u)) {
        cfx_big_shr_bits_eq(&u, 1);
    }

    do {
        /* remove factors of 2 from v */
        while (cfx_big_is_even(&v)) {
            cfx_big_shr_bits_eq(&v, 1);
        }

        /* make sure u <= v */
        if (cfx_big_cmp(&u, &v) > 0) {
            cfx_big_swap(&u, &v);
        }

        cfx_big_sub_eq(&v, &u);
        cfx_big_trim(&v);

    } while (!cfx_big_is_zero(&v));

    /* restore common factors of 2 */
    cfx_big_shl_bits_eq(&u, shift);

    cfx_big_swap(out, &u);
    cfx_big_free(&u);
    cfx_big_free(&v);
}

/* Extended GCD: computes g = gcd(a, b) and Bézout coefficients x, y
 * such that a*x + b*y = g.
 * Uses the iterative extended Euclidean algorithm. */
void cfx_big_xgcd(cfx_big_t *g, cfx_sbig_t *x, cfx_sbig_t *y,
    const cfx_big_t *a, const cfx_big_t *b) {
    /* Handle zero cases */
    if (cfx_big_is_zero(a)) {
        if (g) cfx_big_copy(g, b);
        if (x) cfx_sbig_from_i64(x, 0);
        if (y) cfx_sbig_from_i64(y, cfx_big_is_zero(b) ? 0 : 1);
        return;
    }
    if (cfx_big_is_zero(b)) {
        if (g) cfx_big_copy(g, a);
        if (x) cfx_sbig_from_i64(x, 1);
        if (y) cfx_sbig_from_i64(y, 0);
        return;
    }

    /* Working copies for the Euclidean algorithm */
    cfx_big_t r0, r1, q, r;
    cfx_big_init(&r0);
    cfx_big_init(&r1);
    cfx_big_init(&q);
    cfx_big_init(&r);
    cfx_big_copy(&r0, a);
    cfx_big_copy(&r1, b);

    /* Bézout coefficients: x0, x1, y0, y1 (signed) */
    cfx_sbig_t x0, x1, y0, y1;
    cfx_sbig_init(&x0);
    cfx_sbig_init(&x1);
    cfx_sbig_init(&y0);
    cfx_sbig_init(&y1);
    cfx_sbig_from_i64(&x0, 1);  /* x0 = 1 */
    cfx_sbig_from_i64(&x1, 0);  /* x1 = 0 */
    cfx_sbig_from_i64(&y0, 0);  /* y0 = 0 */
    cfx_sbig_from_i64(&y1, 1);  /* y1 = 1 */

    cfx_sbig_t q_sbig, tmp;
    cfx_sbig_init(&q_sbig);
    cfx_sbig_init(&tmp);

    while (!cfx_big_is_zero(&r1)) {
        /* q = r0 / r1, r = r0 % r1 */
        cfx_big_divrem(&q, &r, &r0, &r1);

        /* r0 = r1, r1 = r */
        cfx_big_swap(&r0, &r1);
        cfx_big_swap(&r1, &r);

        /* Convert q to signed for coefficient updates */
        cfx_sbig_assign_big(&q_sbig, &q, 1);

        /* x0, x1 = x1, x0 - q*x1 */
        cfx_sbig_mul(&tmp, &q_sbig, &x1);
        cfx_sbig_sub(&tmp, &x0, &tmp);
        cfx_sbig_swap(&x0, &x1);
        cfx_sbig_swap(&x1, &tmp);

        /* y0, y1 = y1, y0 - q*y1 */
        cfx_sbig_mul(&tmp, &q_sbig, &y1);
        cfx_sbig_sub(&tmp, &y0, &tmp);
        cfx_sbig_swap(&y0, &y1);
        cfx_sbig_swap(&y1, &tmp);
    }

    /* r0 = gcd, x0 and y0 are the coefficients */
    if (g) cfx_big_swap(g, &r0);
    if (x) cfx_sbig_swap(x, &x0);
    if (y) cfx_sbig_swap(y, &y0);

    cfx_big_free(&r0);
    cfx_big_free(&r1);
    cfx_big_free(&q);
    cfx_big_free(&r);
    cfx_sbig_free(&x0);
    cfx_sbig_free(&x1);
    cfx_sbig_free(&y0);
    cfx_sbig_free(&y1);
    cfx_sbig_free(&q_sbig);
    cfx_sbig_free(&tmp);
}

/* modular inverse via xgcd: out = a^(-1) mod n
 * returns 1 on success, 0 if gcd(a,n) != 1 */
int cfx_big_modinv(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *n) {
    cfx_big_t g;
    cfx_sbig_t x;
    cfx_big_init(&g);
    cfx_sbig_init(&x);

    cfx_big_xgcd(&g, &x, NULL, a, n);

    if (!cfx_big_is_one(&g)) {
        /* no inverse exists */
        cfx_big_free(&g);
        cfx_sbig_free(&x);
        return 0;
    }

    /* x might be negative, need to reduce mod n */
    if (x.sign < 0) {
        /* out = n - |x| */
        cfx_big_copy(out, n);
        cfx_big_sub_eq(out, &x.mag);
    } else {
        cfx_big_copy(out, &x.mag);
    }

    /* ensure out < n (shouldn't be needed but just in case) */
    if (cfx_big_cmp(out, n) >= 0) {
        cfx_big_mod(out, out, n);
    }

    cfx_big_free(&g);
    cfx_sbig_free(&x);
    return 1;
}

/* Pollard-Rho factorization using Montgomery multiplication (Brent's improvement).
 * Returns a non-trivial factor in 'factor', or copies n if n is prime/unfactorable. */
void cfx_big_pollard_rho(cfx_big_t *factor, const cfx_big_t *n) {
    if (cfx_big_cmp_sm(n, 2) < 0) {
        cfx_big_copy(factor, n);
        return;
    }
    if (cfx_big_is_even(n)) {
        cfx_big_from_limb(factor, 2);
        return;
    }
    if (cfx_big_cmp_sm(n, 3) == 0) {
        cfx_big_from_limb(factor, 3);
        return;
    }

    /* monty context for modular arithmetic mod n */
    cfx_big_mont_ctx_t ctx;
    if (cfx_big_mont_ctx_init(&ctx, n) == 0) {
        cfx_big_copy(factor, n);
        return;
    }

    cfx_big_t y, c, x, ys, q, g, one, diff, temp, scratch;
    cfx_big_init(&y);
    cfx_big_init(&c);
    cfx_big_init(&x);
    cfx_big_init(&ys);
    cfx_big_init(&q);
    cfx_big_init(&g);
    cfx_big_init(&one);
    cfx_big_init(&diff);
    cfx_big_init(&temp);
    cfx_big_init(&scratch);

    /* Pre-allocate scratch buffer for squaring operations to avoid repeated mallocs */
    cfx_big_reserve(&scratch, ctx.k + 2);

    /* init one in Montgomery form */
    cfx_big_from_limb(&one, 1);
    cfx_big_mont_to(&one, &one, &ctx);

    /* Try different c values if we fail with the first one */
    static const cfx_limb_t c_vals[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23};
    const size_t num_c_vals = sizeof(c_vals) / sizeof(c_vals[0]);

    for (size_t c_idx = 0; c_idx < num_c_vals; ++c_idx) {
        /* Initialize with new c value */
        cfx_big_from_limb(&y, 2);
        cfx_big_mont_to(&y, &y, &ctx);
        cfx_big_from_limb(&c, c_vals[c_idx]);
        cfx_big_mont_to(&c, &c, &ctx);

        cfx_big_copy(&x, &y);
        cfx_big_copy(&ys, &y);
        cfx_big_copy(&q, &one);
        cfx_big_from_limb(&g, 1);

        size_t nbits = n->n * CFX_LIMB_BITS;
        size_t max_iters = (size_t)1 << (nbits / 4 + 4);  /* e.g., 128-bit -> 2^36 */
#if SIZE_MAX > 0xFFFFFFFFUL
        /* 64-bit: cap at 4 billion iterations */
        if (max_iters > 0x100000000ULL) max_iters = 0x100000000ULL;
#endif

        size_t r = 1;
        const size_t batch = 128;
        size_t iters = 0;

        while (cfx_big_is_one(&g)) {
            cfx_big_copy(&x, &y);

            /* Phase 2 of Brent: advance y by r steps WITHOUT accumulation */
            for (size_t i = 0; i < r; ++i) {
                cfx_big_mont_mul(&scratch, &y, &y, &ctx);
                cfx_big_swap(&y, &scratch);
                cfx_big_add_eq(&y, &c);
                if (cfx_big_cmp(&y, n) >= 0) {
                    cfx_big_sub_eq(&y, n);
                }
            }
            iters += r;

            /* Phase 3 of Brent: advance y by r more steps WITH accumulation */
            for (size_t k = 0; k < r && cfx_big_is_one(&g); k += batch) {
                cfx_big_copy(&ys, &y);
                size_t lim = (batch < r - k) ? batch : (r - k);

                for (size_t i = 0; i < lim; ++i) {
                    /* y = y^2 + c mod n (use scratch to avoid allocation) */
                    cfx_big_mont_mul(&scratch, &y, &y, &ctx);
                    cfx_big_swap(&y, &scratch);
                    cfx_big_add_eq(&y, &c);
                    if (cfx_big_cmp(&y, n) >= 0) {
                        cfx_big_sub_eq(&y, n);
                    }

                    /* diff = |x - y| */
                    if (cfx_big_cmp(&x, &y) >= 0) {
                        cfx_big_copy(&diff, &x);
                        cfx_big_sub_eq(&diff, &y);
                    } else {
                        cfx_big_copy(&diff, &y);
                        cfx_big_sub_eq(&diff, &x);
                    }

                    /* q = q * diff mod n (in Montgomery form) */
                    if (!cfx_big_is_zero(&diff)) {
                        cfx_big_mont_mul(&q, &q, &diff, &ctx);
                    }
                }
                iters += lim;

                /* g = gcd(q, n) */
                cfx_big_mont_from(&temp, &q, &ctx);
                cfx_big_gcd(&g, &temp, n);
            }

            r <<= 1;
            if (iters > max_iters) {
                /* Try next c value */
                break;
            }
        }

        /* If g == n, backtrack step-by-step */
        if (cfx_big_cmp(&g, n) == 0) {
            cfx_big_from_limb(&g, 1);
            cfx_big_mont_from(&x, &x, &ctx);  /* Convert x back from Montgomery */

            for (size_t i = 0; i < batch && cfx_big_is_one(&g); ++i) {
                /* ys = ys^2 + c mod n (use scratch to avoid allocation) */
                cfx_big_mont_mul(&scratch, &ys, &ys, &ctx);
                cfx_big_swap(&ys, &scratch);
                cfx_big_add_eq(&ys, &c);
                if (cfx_big_cmp(&ys, n) >= 0) {
                    cfx_big_sub_eq(&ys, n);
                }

                /* Convert ys from Montgomery form for GCD */
                cfx_big_mont_from(&temp, &ys, &ctx);

                /* diff = |x - ys| */
                if (cfx_big_cmp(&x, &temp) >= 0) {
                    cfx_big_copy(&diff, &x);
                    cfx_big_sub_eq(&diff, &temp);
                } else {
                    cfx_big_copy(&diff, &temp);
                    cfx_big_sub_eq(&diff, &x);
                }

                cfx_big_gcd(&g, &diff, n);
            }
        }

        /* Check if we found a non-trivial factor */
        if (cfx_big_cmp(&g, n) != 0 && !cfx_big_is_one(&g)) {
            cfx_big_copy(factor, &g);
            goto cleanup;
        }
    }

    /* All c values exhausted - return n (failure) */
    cfx_big_copy(factor, n);

cleanup:
    cfx_big_free(&y);
    cfx_big_free(&c);
    cfx_big_free(&x);
    cfx_big_free(&ys);
    cfx_big_free(&q);
    cfx_big_free(&g);
    cfx_big_free(&one);
    cfx_big_free(&diff);
    cfx_big_free(&temp);
    cfx_big_free(&scratch);
    cfx_big_mont_ctx_free(&ctx);
}
