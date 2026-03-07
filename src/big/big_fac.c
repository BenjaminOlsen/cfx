/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include "cfx/fac.h"
#include "cfx/primes.h"
#include "cfx/ecm.h"


/* Materialize factorization into cfx_big_t */
void cfx_big_from_fac(cfx_big_t *b, const cfx_fac_t *f) {
    cfx_big_from_limb(b, 1);
    for (size_t i = 0; i < f->len; i++){
        cfx_big_expmul_prime(b, f->data[i].p, f->data[i].e);
    }
}

/* ******************************************************************************* */
/* bucket-carry strategy for making bigs out of fac : */
#define CFX_FAC_BUCKETS 128 /* plenty: log2(bitlen(n!)) */

static void bucket_init(cfx_big_t buckets[CFX_FAC_BUCKETS], uint8_t used[CFX_FAC_BUCKETS]) {
    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        cfx_big_init(&buckets[i]);
        used[i] = 0;
    }
}

static void bucket_free(cfx_big_t buckets[CFX_FAC_BUCKETS], uint8_t used[CFX_FAC_BUCKETS]) {
    (void)used;
    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        cfx_big_free(&buckets[i]);
    }
}

/* Inserts x into buckets by multiplying upward while occupied. x is consumed (moved from). */
static void bucket_insert(cfx_big_t buckets[CFX_FAC_BUCKETS], uint8_t used[CFX_FAC_BUCKETS], cfx_big_t *x) {
    size_t level = 0;

    if (cfx_big_is_one(x)) {
        cfx_big_free(x);
        cfx_big_init(x);
        return;
    }

    while (level < CFX_FAC_BUCKETS) {
        if (!used[level]) {
            cfx_big_move(&buckets[level], x);  /* buckets[level] takes ownership */
            used[level] = 1;
            return;
        }

        /* x *= buckets[level]; clear bucket; carry upward */
        cfx_big_mul_auto(x, &buckets[level]);    /* x = x * buckets[level] */
        cfx_big_free(&buckets[level]);
        cfx_big_init(&buckets[level]);
        used[level] = 0;
        level++;
    }

    /* ---------> todo: increase CFX_FAC_BUCKETS. */
}

void cfx_big_from_fac_fast(cfx_big_t *out, const cfx_fac_t *f) {
    cfx_big_t buckets[CFX_FAC_BUCKETS];
    uint8_t used[CFX_FAC_BUCKETS];

    bucket_init(buckets, used);


    for (size_t i = 0; i < f->len; i++) {
        const cfx_limb_t p = f->data[i].p;
        const cfx_limb_t e = f->data[i].e;
        if (!e) continue;

        cfx_big_t x;
        cfx_big_init(&x);
        cfx_big_from_limb(&x, 1);
        cfx_big_expmul_prime(&x, p, e);  /* x = p^e */

        bucket_insert(buckets, used, &x);

        cfx_big_free(&x);   /* x was moved into a bucket or freed; ensure no leak if insert didn't consume */
    }


    cfx_big_from_limb(out, 1);
    for (size_t level = 0; level < CFX_FAC_BUCKETS; level++) {
        if (used[level]) {
            cfx_big_mul_auto(out, &buckets[level]); /* out *= buckets[level] */
        }
    }

    bucket_free(buckets, used);
}

/* x > 0 */
static size_t floor_log2_size_t(size_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return (size_t)(8u * sizeof(size_t) - 1u - (unsigned)__builtin_clzl(x));
#else
    size_t r = 0;
    while (x >>= 1) r++;
    return r;
#endif
}

static void bucket_clear(cfx_big_t *b) {
    cfx_big_free(b); cfx_big_init(b);
}


void cfx_big_from_fac_faster(cfx_big_t *out, const cfx_fac_t *f) {
    cfx_big_t buckets[CFX_FAC_BUCKETS];
    uint8_t used[CFX_FAC_BUCKETS];

    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        cfx_big_init(&buckets[i]);
        used[i] = 0;
    }

    for (size_t i = 0; i < f->len; i++) {
        cfx_limb_t p = f->data[i].p;
        cfx_limb_t e = f->data[i].e;
        if (!e) continue;

        cfx_big_t x;
        cfx_big_init(&x);
        cfx_big_pow_sm(&x, p, e);

        if (cfx_big_is_one(&x)) {
            cfx_big_free(&x);
            continue;
        }

        size_t bl = cfx_big_bitlen(&x);
        size_t level = (bl > 0) ? floor_log2_size_t(bl) : 0;
        if (level >= CFX_FAC_BUCKETS) level = CFX_FAC_BUCKETS - 1;

        for (;;) {
            if (!used[level]) {
                cfx_big_move(&buckets[level], &x);
                used[level] = 1;
                break;
            }
            /* x *= buckets[level], clear bucket, carry up */
            cfx_big_mul_auto(&x, &buckets[level]);
            bucket_clear(&buckets[level]);
            used[level] = 0;

            if (level + 1 < CFX_FAC_BUCKETS) {
                ++level;
            } else {
                /* extremely unlikely: fold into last bucket */
                if (!used[level]) {
                    cfx_big_move(&buckets[level], &x);
                    used[level] = 1;
                } else {
                    cfx_big_mul_auto(&buckets[level], &x);
                }
                break;
            }
        }

        cfx_big_free(&x);
    }

    cfx_big_from_limb(out, 1);
    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        if (used[i]) cfx_big_mul_auto(out, &buckets[i]);
        cfx_big_free(&buckets[i]);
    }
}


/* Helper: check if a big integer fits in 64 bits */
static int big_fits_in_64(const cfx_big_t *b) {
#if CFX_LIMB_BITS == 64
    return b->n <= 1;
#elif CFX_LIMB_BITS == 32
    return b->n <= 2;
#endif
}

/* Helper: convert big integer to uint64 (assumes it fits) */
static uint64_t big_to_u64(const cfx_big_t *b) {
    if (b->n == 0) return 0;
#if CFX_LIMB_BITS == 64
    return b->limb[0];
#elif CFX_LIMB_BITS == 32
    uint64_t val = b->limb[0];
    if (b->n >= 2) val |= ((uint64_t)b->limb[1] << 32);
    return val;
#endif
}

/*
 * Factorize a big integer into prime factors.
 * Strategy:
 *   1. Trial division by small primes from cfx_primes[]
 *   2. If remainder fits in 64 bits, delegate to cfx_fac_from_u64
 *   3. Otherwise, use Pollard-Rho to find factors recursively
 *   4. If we encounter a prime > 64 bits, return incomplete
 *
 * Returns:
 *   0  - complete factorization
 *   1  - incomplete (remainder holds unfactored prime > 64 bits)
 *  -1  - error
 */
int cfx_big_to_fac(cfx_fac_t *f, const cfx_big_t *b, cfx_big_t *remainder) {
    cfx_fac_init(f);

    if (remainder) {
        cfx_big_from_limb(remainder, 1);
    }

    if (cfx_big_is_zero(b)) return 0;
    if (b->n == 1 && b->limb[0] == 1) return 0;

    /* work on a copy */
    cfx_big_t work;
    cfx_big_init(&work);
    cfx_big_copy(&work, b);

    /* Trial division by small primes */
    for (size_t i = 0; i < cfx_primes_len && !cfx_big_is_zero(&work); ++i) {
        uint32_t p = cfx_primes[i];

        /* If work < p, we're done with trial division */
        if (work.n == 1 && work.limb[0] < p) break;

        /* Count exponent */
        uint32_t e = 0;
        while (cfx_big_mod_sm(&work, p) == 0) {
            cfx_big_div_sm_eq(&work, p);
            cfx_big_trim(&work);
            e++;
        }
        if (e > 0) {
            cfx_fac_push(f, (uint64_t)p, e);
        }

        if (cfx_big_is_one(&work)) break;
    }

    /* If work == 1, factorization is complete */
    if (cfx_big_is_zero(&work) || cfx_big_is_one(&work)) {
        cfx_big_free(&work);
        return 0;
    }

    /* use a stack for iterative factorization */
    cfx_big_t stack[64];  /* Enough for any reasonable factorization */
    size_t stack_top = 0;
    stack_top = 1;
    cfx_big_init(&stack[0]);
    cfx_big_swap(&stack[0], &work);

    int result = 0;  /* 0 = complete */

    while (stack_top > 0) {
        cfx_big_t *cur = &stack[--stack_top];

        if (cfx_big_is_one(cur) || cfx_big_is_zero(cur)) {
            cfx_big_free(cur);
            continue;
        }

        /* If it fits in 64 bits, use fast 64-bit factorization */
        if (big_fits_in_64(cur)) {
            uint64_t val = big_to_u64(cur);
            cfx_fac_t fac64;
            if (cfx_fac_from_u64(&fac64, val) == 0) {
                for (size_t i = 0; i < fac64.len; ++i) {
                    cfx_fac_push(f, fac64.data[i].p, fac64.data[i].e);
                }
            }
            cfx_fac_free(&fac64);
            cfx_big_free(cur);
            continue;
        }

        /* Check if it's prime */
        if (cfx_big_is_prime(cur)) {
            /* Prime > 64 bits - store in fac_t's big_primes array */
            cfx_fac_push_big(f, cur, 1);
            cfx_big_free(cur);
            continue;
        }

        /* Use Pollard-Rho to find a factor */
        cfx_big_t factor;
        cfx_big_init(&factor);
        cfx_big_pollard_rho(&factor, cur);

        /* If Pollard-Rho failed (returned n), try ECM */
        if (cfx_big_cmp(&factor, cur) == 0) {
            cfx_big_from_limb(&factor, 0);
            if (cfx_ecm_factor_auto(&factor, cur)) {
                /* ECM found a factor - fall through to push it */
            } else {
                /* ECM also failed - give up on this composite */
                if (remainder) {
                    cfx_big_mul_eq(remainder, cur);
                }
                result = 1;  /* incomplete */
                cfx_big_free(&factor);
                cfx_big_free(cur);
                continue;
            }
        }

        /* Push both factors onto the stack */
        cfx_big_t quotient;
        cfx_big_init(&quotient);
        cfx_big_copy(&quotient, cur);
        cfx_big_divrem_eq(&quotient, &factor, NULL);
        cfx_big_trim(&quotient);

        /* Free cur BEFORE pushing, since cur points into the stack array
        * and we're about to reuse that slot. This avoids double-free. */
        cfx_big_free(cur);

        /* Push factor and quotient */
        if (stack_top < 62) {
            cfx_big_init(&stack[stack_top]);
            cfx_big_swap(&stack[stack_top], &factor);
            stack_top++;

            cfx_big_init(&stack[stack_top]);
            cfx_big_swap(&stack[stack_top], &quotient);
            stack_top++;
        }

        cfx_big_free(&factor);
        cfx_big_free(&quotient);
    }

    cfx_big_free(&work);

    /* Sort and coalesce the factorization */
    /* (cfx_fac_push maintains sorted order, but we might have duplicates) */
    if (f->len > 1) {
        /* Coalesce duplicate primes */
        cfx_fac_t coalesced;
        cfx_fac_init(&coalesced);

        /* Sort by prime (bubble sort - factors are usually few) */
        for (size_t i = 0; i < f->len; ++i) {
            for (size_t j = i + 1; j < f->len; ++j) {
                if (f->data[i].p > f->data[j].p) {
                    cfx_pf_t tmp = f->data[i];
                    f->data[i] = f->data[j];
                    f->data[j] = tmp;
                }
            }
        }

        /* Coalesce adjacent equal primes */
        uint64_t cur_p = f->data[0].p;
        uint32_t cur_e = f->data[0].e;
        for (size_t i = 1; i < f->len; ++i) {
            if (f->data[i].p == cur_p) {
                cur_e += f->data[i].e;
            } else {
                cfx_fac_push(&coalesced, cur_p, cur_e);
                cur_p = f->data[i].p;
                cur_e = f->data[i].e;
            }
        }
        cfx_fac_push(&coalesced, cur_p, cur_e);

        /* Transfer big_primes from f to coalesced before freeing */
        coalesced.big_primes = f->big_primes;
        coalesced.big_len = f->big_len;
        coalesced.big_cap = f->big_cap;
        /* Clear f's big_primes pointers so cfx_fac_free doesn't free them */
        f->big_primes = NULL;
        f->big_len = 0;
        f->big_cap = 0;

        /* Swap */
        cfx_fac_free(f);
        *f = coalesced;
    }

    return result;
}

int cfx_fac_from_big(cfx_fac_t *fac, const cfx_big_t *in) {
    (void)fac;
    (void)in;
    /* TODO: implement big integer factorization into cfx_fac_t
     * See cfx_big_to_fac() for a working implementation that uses
     * trial division + Pollard-Rho. This function would need helper
     * functions like cfx_big_strip_small, cfx_big_fits_u64, etc. */
    return -1;
}
