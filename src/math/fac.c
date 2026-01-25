/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/fac.h"
#include "cfx/algo.h"
#include "cfx/vec.h"
#include "cfx/macros.h"
#include "cfx/primes.h"
#include "cfx/rand.h"
#include "cfx/compat.h"
#include "cfx/big.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

void cfx_fac_print(cfx_fac_t *f) {
    if (f->len == 0 && f->big_len == 0) {
        printf("(empty)\n");
        return;
    }
    /* Print small primes */
    for (size_t k = 0; k < f->len; ++k) {
        cfx_pf_t *pf = &f->data[k];
        printf("%" PRIu64 "^%" PRIu32, pf->p, pf->e);
        if (k + 1 < f->len || f->big_len > 0) {
            printf(" * ");
        }
    }
    /* Print big primes */
    for (size_t k = 0; k < f->big_len; ++k) {
        cfx_big_pf_t *pf = &f->big_primes[k];
        char *p_str = cfx_big_dec_alloc(pf->p, NULL);
        printf("%s^%" PRIu32, p_str, pf->e);
        free(p_str);
        if (k + 1 < f->big_len) {
            printf(" * ");
        }
    }
    printf("\n");
}

void cfx_fac_init(cfx_fac_t *f) {
    f->data = NULL;
    f->len = 0;
    f->cap = 0;
    f->big_primes = NULL;
    f->big_len = 0;
    f->big_cap = 0;
}

void cfx_fac_clear(cfx_fac_t *f) {
    f->len = 0;
}

void cfx_fac_free(cfx_fac_t *f) {
    free(f->data);
    f->data = NULL;
    f->len = 0;
    f->cap = 0;
    /* Free big primes */
    for (size_t i = 0; i < f->big_len; ++i) {
        cfx_big_free(f->big_primes[i].p);
        free(f->big_primes[i].p);
    }
    free(f->big_primes);
    f->big_primes = NULL;
    f->big_len = 0;
    f->big_cap = 0;
}

static inline int _cfx_mul_zu_ok(size_t a, size_t b, size_t *out) {
    if (a == 0 || b == 0) {
        *out = 0; return -1;
    }
    if (a > SIZE_MAX / b) return -1;
    *out = a * b;
    return 0;
}

void cfx_fac_reserve(cfx_fac_t *f, size_t req_cap) {
    if (req_cap <= f->cap) {
        return;
    }
    size_t new_cap = f->cap ? 2 * f->cap : 16;
    if (new_cap < req_cap) {
        new_cap = req_cap;
    }
    size_t bytes;
    if (_cfx_mul_zu_ok(new_cap, sizeof(cfx_pf_t), &bytes) != 0) {
        fprintf(stderr, "cfx_fac_reserve: allocation size overflow\n");
        abort();
    }
    void *tmp = realloc(f->data, bytes);
    if (!tmp) {
        fprintf(stderr, "cfx_fac_reserve: out of memory (requested %zu bytes)\n", bytes);
        abort();
    }
    f->data = (cfx_pf_t *)tmp;
    f->cap = new_cap;
}

int cfx_fac_push(cfx_fac_t *f, uint64_t p, uint32_t e) {
    if (e == 0) return -1;
    cfx_fac_reserve(f, f->len + 1);
    ++f->len;
    cfx_pf_t *pf = &f->data[f->len - 1];
    pf->p = p;
    pf->e = e;
    return 0;
}

int cfx_fac_push_big(cfx_fac_t *f, const struct cfx_big *p, uint32_t e) {
    if (e == 0) return -1;
    /* Reserve space for big primes */
    if (f->big_len >= f->big_cap) {
        size_t new_cap = f->big_cap ? 2 * f->big_cap : 4;
        void *tmp = realloc(f->big_primes, new_cap * sizeof(cfx_big_pf_t));
        if (!tmp) return -1;
        f->big_primes = (cfx_big_pf_t *)tmp;
        f->big_cap = new_cap;
    }
    /* Allocate and copy the big prime */
    cfx_big_t *pcopy = (cfx_big_t *)malloc(sizeof(cfx_big_t));
    if (!pcopy) return -1;
    cfx_big_init(pcopy);
    cfx_big_copy(pcopy, (const cfx_big_t *)p);
    f->big_primes[f->big_len].p = pcopy;
    f->big_primes[f->big_len].e = e;
    ++f->big_len;
    return 0;
}

/* Deep copy: dst becomes a copy of src, using cfx_fac_push for each element. */
void cfx_fac_copy(cfx_fac_t *dst, const cfx_fac_t *src) {
    if (dst == src) return;
    cfx_fac_free(dst);
    cfx_fac_init(dst);
    cfx_fac_reserve(dst, src->len);

    for (size_t i = 0; i < src->len; ++i) {
        cfx_fac_push(dst, src->data[i].p, src->data[i].e);
    }
}

/* dst += src */
void cfx_fac_add(cfx_fac_t *dst, const cfx_fac_t *src) {
    cfx_fac_t out;
    cfx_fac_init(&out);
    cfx_fac_reserve(&out, src->len + dst->len);

    /* Assumes data[] is sorted by p */
    size_t i = 0;
    size_t j = 0;

    while (i < dst->len || j < src->len) {
        cfx_pf_t *pf1 = &dst->data[i];
        cfx_pf_t *pf2 = &src->data[j];

        if (j == src->len || (i < dst->len && pf1->p < pf2->p)) {
            /* p is in dst, but not in src */
            cfx_fac_push(&out, pf1->p, pf1->e);
            ++i;
        } else if (i == dst->len || (j < src->len && pf1->p > pf2->p)) {
            /* p is in src, but not in dst */
            cfx_fac_push(&out, pf2->p, pf2->e);
            ++j;
        } else {
            /* p is in src and dst */
            uint64_t p = pf1->p;
            uint32_t e = pf1->e + pf2->e;
            if (e < pf1->e) {
                fprintf(stderr, "cfx_fac_add: exponent overflow for prime %" PRIu64 "\n", p);
                abort();
            }
            if (e) cfx_fac_push(&out, p, e);
            ++i;
            ++j;
        }
    }
    free(dst->data);
    *dst = out;
}

/* dst -= src */
void cfx_fac_sub(cfx_fac_t *dst, const cfx_fac_t *src) {
    cfx_fac_t out;

    cfx_fac_init(&out);
    cfx_fac_reserve(&out, dst->len);

    size_t i = 0;
    size_t j = 0;

    while (i < dst->len || j < src->len) {
        cfx_pf_t *pf1 = &dst->data[i];
        cfx_pf_t *pf2 = &src->data[j];

        if (j == src->len || (i < dst->len && pf1->p < pf2->p)) {
            /* p is in dst but not src - no change */
            cfx_fac_push(&out, pf1->p, pf1->e);
            ++i;
        } else if (i == dst->len || (j < src->len && pf1->p > pf2->p)) {
            /* p is in src but not dst - dst does not divide src */
            fprintf(stderr, "cfx_fac_sub: underflow (prime %" PRIu64 " in src but not dst)\n", pf2->p);
            abort();
        } else {
            /* p is in both src and dst: */
            uint64_t p = pf1->p;
            uint32_t e;
            if (pf1->e > pf2->e) { /* dst divides src */
                e = pf1->e - pf2->e;
                cfx_fac_push(&out, p, e);
            } else if (pf2->e > pf1->e) { /* dst does not divide src */
                fprintf(stderr, "cfx_fac_sub: underflow (exponent for prime %" PRIu64 ")\n", p);
                abort();
            }
            /* else, e == 0, p is cancelled out! */
            ++i;
            ++j;
        }
    }
    free(dst->data);
    *dst = out;
}

/** calculate the factorial of n.
 * we pass in a list of primes to not have to calculate it on every call of this function -
 * precondition: primes is sorted strictly increasing!
 */
int cfx_fac_factorial(cfx_fac_t *f, cfx_limb_t n, const cfx_vec_t *primes) {
    int ret = 0;
    if (n == 0) {
        cfx_fac_init(f);
        goto end;
    }
    for (size_t i = 0; i < primes->size; ++i) {
        cfx_limb_t p = primes->data[i];
        if (p > n) break;
        cfx_limb_t e = cfx_legendre(n, p);
        if (e) {
            ret = cfx_fac_push(f, p, e);
            if (ret != 0) goto end;
        }
    }
end:
    return ret;
}

/* Factorization of C(n,k) = n! / (k! (n-k)!) */
cfx_fac_t cfx_fac_binom(cfx_limb_t n, cfx_limb_t k){
    if (k > n) {
        cfx_fac_t z;
        cfx_fac_init(&z);
        return z;
    }
    if (k > n-k) k = n-k;
    cfx_vec_t primes = cfx_sieve_primes(n);
    cfx_fac_t fn;
    cfx_fac_init(&fn);
    cfx_fac_factorial(&fn, n, &primes);
    cfx_fac_t fk;
    cfx_fac_init(&fk);
    cfx_fac_factorial(&fk, k, &primes);
    cfx_fac_t fnk;
    cfx_fac_init(&fnk);
    cfx_fac_factorial(&fnk, n-k, &primes);
    cfx_fac_sub(&fn, &fk);
    cfx_fac_sub(&fn, &fnk);
    cfx_vec_free(&primes);
    return fn; /* sorted, coalesced */
}

/* Compare function for qsort */
static int cmp_u64_fac(const void *a, const void *b) {
    uint64_t x = *(const uint64_t *)a;
    uint64_t y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

/* Public entry: factor n (uint64_t) into fac (coalesced).
 * Returns 0 on success, -1 on error (n==0). */
int cfx_fac_from_u64(cfx_fac_t *fac, uint64_t n) {
    cfx_fac_init(fac);
    if (n == 0) return -1;
    if (n == 1) return 0;

    /* Use a local array for 64-bit prime storage during computation */
    uint64_t plist[256];  /* enough for any 64-bit factorization */
    size_t pcount = 0;

    /* 1) Strip small primes using precomputed list */
    for (size_t i = 0; i < cfx_primes_len && i < 54; ++i) {
        uint64_t p = cfx_primes[i];
        while ((n % p) == 0) {
            if (pcount < 256) plist[pcount++] = p;
            n /= p;
        }
        if (n == 1) break;
    }

    /* 2) Iterative stack of composites to split */
    uint64_t st[64]; /* plenty for 64-bit factoring */
    size_t top = 0;
    if (n > 1) st[top++] = n;

    while (top) {
        uint64_t m = st[--top];
        if (m == 1) continue;

        /* If small enough, finish by trial division quickly */
        if (m < 1ull<<20) {
            for (uint64_t p = 3; p*p <= m; p += 2) {
                while ((m % p) == 0) {
                    if (pcount < 256) plist[pcount++] = p;
                    m /= p;
                }
            }
            if (m > 1 && pcount < 256) plist[pcount++] = m;
            continue;
        }

        if (cfx_is_prime_u64(m)) {
            if (pcount < 256) plist[pcount++] = m;
            continue;
        }

        /* Split with Brent–Rho; allow a few retries with different seeds */
        uint64_t d = 0;
        for (int tries = 0; tries < 8; ++tries) {
            cfx_srand(0xC0FFEEu + (uint32_t)tries);
            d = cfx_pollard_rho_brent(m);
            if (d > 1 && d < m && (m % d) == 0) break;
            d = 0;
        }

        if (d == 0) {
            /* Last resort: bounded trial division */
            uint64_t found = 0;
            for (uint64_t p = 3; p*p <= m; p += 2) {
                if (m % p == 0) {
                    found = p; break;
                }
            }
            if (!found) {
                if (pcount < 256) plist[pcount++] = m;
                continue;
            }
            d = found;
        }

        /* Push the two factors back for further processing */
        uint64_t a = d;
        uint64_t b = m / d;
        if (a < b) {
            uint64_t t=a; a=b; b=t;
        }
        st[top++] = a;
        st[top++] = b;
    }

    /* 3) Sort & coalesce -> push to fac */
    qsort(plist, pcount, sizeof(uint64_t), cmp_u64_fac);

    for (size_t i = 0; i < pcount; ) {
        uint64_t p = plist[i];
        uint32_t e = 1;
        size_t j = i + 1;
        while (j < pcount && plist[j] == p) {
            ++e; ++j;
        }
        cfx_fac_push(fac, p, e);
        i = j;
    }
    return 0;
}

