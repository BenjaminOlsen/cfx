#include "cfx/fac.h"
#include "cfx/vector.h"
#include "cfx/algo.h"
#include "cfx/error.h"
#include "cfx/macros.h"
#include "cfx/primes.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

void cfx_fac_print(cfx_fac_t* f) {
    CFX_PRINT_DBG("f: %p, cap: %zu, len: %zu\n", (void*)f, f->cap, f->len);
    for (size_t k = 0; k < f->len-1; ++k) {
        cfx_pf_t* pf = &f->data[k];
        CFX_PRINT_DBG("%llu^%llu * ", pf->p, pf->e);
    }
    cfx_pf_t* pf = &f->data[f->len-1];
    CFX_PRINT_DBG("%llu^%llu\n", pf->p, pf->e);
}

void cfx_fac_init(cfx_fac_t* f) {
    f->data = NULL;
    f->len = 0;
    f->cap = 0;
}

void cfx_fac_free(cfx_fac_t* f) {
    free(f->data);
    cfx_fac_init(f);
}

static inline int _cfx_mul_zu_ok(size_t a, size_t b, size_t* out) {
    if (a == 0 || b == 0) { *out = 0; return CFX_OK; }
    if (a > SIZE_MAX / b) return CFX_EOVERFLOW;
    *out = a * b;
    return CFX_OK;
}

int cfx_fac_reserve(cfx_fac_t* f, size_t req_cap) {
    if (req_cap <= f->cap) {
        return CFX_OK;
    }
    size_t new_cap = f->cap ? 2 * f->cap : 16;
    if (new_cap < req_cap) {
        new_cap = req_cap;
    }
    // CFX_PRINT_DBG("cfx_fac_reserve %p requested cap: %zu, new cap: %zu\n", f, req_cap, new_cap);
    size_t bytes;
    if (_cfx_mul_zu_ok(new_cap, sizeof(cfx_pf_t), &bytes) != CFX_OK) {
        return CFX_ENOMEM;
    }
    f->data = (cfx_pf_t*)realloc(f->data, bytes);
    f->cap = new_cap;
    return CFX_OK;
}

int cfx_fac_push(cfx_fac_t* f, uint64_t p, uint64_t e) {
    int ret = CFX_OK;
    if (e == 0) goto end;
    ret = cfx_fac_reserve(f, f->len + 1);
    if (ret != CFX_OK) goto end;
    ++f->len;
    cfx_pf_t* pf = &f->data[f->len - 1];
    pf->p = p;
    pf->e = e;
end:
    return ret;
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

// dst += src
void cfx_fac_add(cfx_fac_t* dst, cfx_fac_t* src) {
    cfx_fac_t out;
    cfx_fac_init(&out);
    cfx_fac_reserve(&out, src->len + dst->len);

    // Assumes data[] is sorted by p
    size_t i = 0;
    size_t j = 0;

    while (i < dst->len || j < src->len) {
        cfx_pf_t* pf1 = &dst->data[i];
        cfx_pf_t* pf2 = &src->data[j];

        if (j == src->len || (i < dst->len && pf1->p < pf2->p)) {
            // p is in dst, but not in src
            cfx_fac_push(&out, pf1->p, pf1->e);
            ++i;
        } else if (i == dst->len || (j < src->len && pf1->p > pf2->p)) {
            // p is in src, but not in dst
            cfx_fac_push(&out, pf2->p, pf2->e);
            ++j;
        } else {
            // p is in src and dst
            uint64_t p = pf1->p;
            uint64_t e = pf1->e + pf2->e;
            if (e) cfx_fac_push(&out, p, e);
            ++i; 
            ++j;
        }
    }
    free(dst->data);
    *dst = out;
}

/* dst -= src */
void cfx_fac_sub(cfx_fac_t* dst, cfx_fac_t* src) {
    cfx_fac_t out;
    
    cfx_fac_init(&out);
    cfx_fac_reserve(&out, dst->len);

    size_t i = 0;
    size_t j = 0;

    while (i < dst->len || j < src->len) {
        cfx_pf_t* pf1 = &dst->data[i];
        cfx_pf_t* pf2 = &src->data[j];

        if (j == src->len || (i < dst->len && pf1->p < pf2->p)) {
            // p is in dst but not src - no change
            cfx_fac_push(&out, pf1->p, pf1->e);
            ++i;
        } else if (i == dst->len || (j < src->len && pf1->p > pf2->p)) {
            // p is in src but not dst - dst does not divide src!!!!
            assert(0 && "cfx_fac_sub underflow");
            j++;
        } else {
            // p is in both src and dst:
            uint64_t p = pf1->p;
            uint64_t e;
            if (pf1->e > pf2->e) { // dst divides src
                e = pf1->e - pf2->e;
                cfx_fac_push(&out, p, e);
            } else if (pf2->e > pf1->e) { // dst does not divide src!
                assert(0 && "cfx_fac_sub underflow");
            } 
            // else, e == 0, p is cancelled out!ç
            ++i;
            ++j;
        }
    }
    free(dst->data);
    *dst = out;
}

/* calculate the factorization of n.
we pass in a list of primes to not have to calculate it on every call of this function */
int cfx_fac_factorial(cfx_fac_t *f, uint64_t n, const cfx_vec_t *primes) {
    int ret = CFX_OK;
    for (size_t i = 0; i < primes->size; ++i) {
        uint64_t p = primes->data[i];
        if (p > n) break;
        uint64_t e = cfx_legendre(n, p); 
        if (e) {
            ret = cfx_fac_push(f, p, e);
            if (ret != CFX_OK) goto end;
        }
    }
end:
    return ret;
}

/* Factorization of C(n,k) = n! / (k! (n-k)!) */
cfx_fac_t cfx_fac_binom(uint64_t n, uint64_t k){
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


static int cmp_u64(const void* a, const void* b) {
    uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x < y) ? -1 : (x > y);
}

// Public entry: factor n (uint64_t) into fac (coalesced).
// Returns 1 on success; 0 for n==0 (degenerate) or allocation failures.
int cfx_fac_from_u64(cfx_fac_t* fac, uint64_t n) {
    if (n == 0) return 0;
    if (n == 1) return 1;

    cfx_vec_t primes; 
    cfx_vec_init(&primes);

    /* 1) Strip tiny (<255) primes fast */
    for (size_t i = 0; i < 54; ++i) {
        uint64_t p = cfx_primes[i];
        while ((n % p) == 0) {
            cfx_vec_push(&primes, p);
            n /= p;
        }
        if (n == 1) break;
    }

    /* 2) Iterative stack of composites to split */
    uint64_t st[64]; // plenty for 64-bit factoring
    size_t top = 0;
    if (n > 1) st[top++] = n;

    while (top) {
        uint64_t m = st[--top];
        if (m == 1) continue;

        // If small enough, finish by trial division quickly
        if (m < 1ull<<20) {
            // trial divide from 3 upward (skip even—shouldn’t be even here often)
            for (uint64_t p = 3; p*p <= m; p += 2) {
                while ((m % p) == 0) {
                    cfx_vec_push(&primes, p); m /= p;
                }
            }
            if (m > 1) cfx_vec_push(&primes, m);
            continue;
        }

        if (cfx_is_prime_u64(m)) {
            cfx_vec_push(&primes, m);
            continue;
        }

        // 3) Split with Brent–Rho; allow a few retries with different seeds
        uint64_t d = 0;
        for (int tries = 0; tries < 8; ++tries) {
            // If your rho_brent uses rand(), reseed to diversify:
            srand(0xC0FFEEu + tries);
            d = cfx_rho_brent(m);
            if (d > 1 && d < m && (m % d) == 0) break;
            d = 0;
        }

        if (d == 0) {
            // Very rare fallback: as a last resort, do a bounded trial division up to 2^32
            // (For 64-bit inputs this almost never triggers with a decent rho.)
            // You could also call SQUFOF here if you have it.
            uint64_t found = 0;
            for (uint64_t p = 3; p*p <= m; p += 2) {
                if (m % p == 0) { found = p; break; }
            }
            if (!found) { // if we still failed, treat as prime to avoid infinite loop
                cfx_vec_push(&primes, m);
                continue;
            }
            d = found;
        }

        // Push the two factors back for further processing (order doesn’t matter)
        uint64_t a = d;
        uint64_t b = m / d;
        // Prefer to push the larger first (optional), purely heuristic
        if (a < b) { uint64_t t=a; a=b; b=t; }
        st[top++] = a;
        st[top++] = b;
    }

    // 4) Sort & coalesce → emit (p, e)
    qsort(primes.data, primes.size, sizeof(uint64_t), cmp_u64);

    for (size_t i = 0; i < primes.size; ) {
        uint64_t p = primes.data[i];
        uint32_t e = 1;
        size_t j = i + 1;
        while (j < primes.size && primes.data[j] == p) { ++e; ++j; }
        cfx_fac_push(fac, p, e);
        i = j;
    }

    cfx_vec_free(&primes);
    return 1;
}

