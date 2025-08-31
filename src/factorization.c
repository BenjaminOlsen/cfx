#include "cfx/factorization.h"
#include "cfx/vector.h"
#include "cfx/primes.h"

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

void cfx_fac_print(cfx_factorization_t* f) {
    printf("f: %p, cap: %zu, len: %zu\n", f, f->cap, f->len);
    for (size_t k = 0; k < f->len-1; ++k) {
        cfx_pf_t* pf = &f->data[k];
        printf("%u^%u * ", pf->p, pf->e);
    }
    cfx_pf_t* pf = &f->data[f->len-1];
    printf("%u^%u\n", pf->p, pf->e);
}

void cfx_fac_init(cfx_factorization_t* f) {
    f->data = NULL;
    f->len = 0;
    f->cap = 0;
}

void cfx_fac_free(cfx_factorization_t* f) {
    free(f->data);
    cfx_fac_init(f);
}

void cfx_fac_reserve(cfx_factorization_t* f, size_t req_cap) {
    if (req_cap <= f->cap) {
        return;
    }
    size_t new_cap = f->cap ? 2 * f->cap : 16;
    if (new_cap < req_cap) {
        new_cap = req_cap;
    }
    // printf("cfx_fac_reserve %p requested cap: %zu, new cap: %zu\n", f, req_cap, new_cap);
    f->data = (cfx_pf_t*)realloc(f->data, new_cap * sizeof(cfx_pf_t));
    f->cap = new_cap;
}

void cfx_fac_push(cfx_factorization_t* f, uint32_t p, uint32_t e) {
    if (e == 0) return;
    cfx_fac_reserve(f, f->len + 1);
    ++f->len;
    cfx_pf_t* pf = &f->data[f->len - 1];
    pf->p = p;
    pf->e = e;
    
}

/* Deep copy: dst becomes a copy of src, using cfx_fac_push for each element. */
void cfx_fac_copy(cfx_factorization_t *dst, const cfx_factorization_t *src) {
    if (dst == src) return;  
    cfx_fac_free(dst);
    cfx_fac_init(dst);
    cfx_fac_reserve(dst, src->len);

    for (size_t i = 0; i < src->len; ++i) {
        cfx_fac_push(dst, src->data[i].p, src->data[i].e);
    }
}

// dst += src
void cfx_fac_add(cfx_factorization_t* dst, cfx_factorization_t* src) {
    cfx_factorization_t out;
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
            uint32_t p = pf1->p;
            uint32_t e = pf1->e + pf2->e;
            if (e) cfx_fac_push(&out, p, e);
            ++i; 
            ++j;
        }
    }
    free(dst->data);
    *dst = out;
}

/* dst -= src */
void cfx_fac_sub(cfx_factorization_t* dst, cfx_factorization_t* src) {
    cfx_factorization_t out;
    
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
            uint32_t p = pf1->p;
            uint32_t e;
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
cfx_factorization_t cfx_fac_factorial(uint32_t n, const cfx_vec_t *primes) {
    cfx_factorization_t f;
    cfx_fac_init(&f);
    for (size_t i = 0; i < primes->size; ++i) {
        uint32_t p = primes->data[i];
        if (p > n) break;
        uint32_t e = cfx_legendre(n, p);
        if (e) cfx_fac_push(&f, p, e);
    }
    return f;
}

/* Factorization of C(n,k) = n! / (k! (n-k)!) */
cfx_factorization_t cfx_fac_binom(uint32_t n, uint32_t k){
    if (k > n) { 
        cfx_factorization_t z;
        cfx_fac_init(&z);
        return z; 
    }
    if (k > n-k) k = n-k;
    cfx_vec_t primes = cfx_sieve_primes(n);
    cfx_factorization_t fn = cfx_fac_factorial(n, &primes);
    cfx_factorization_t fk = cfx_fac_factorial(k, &primes);
    cfx_factorization_t fnk = cfx_fac_factorial(n-k, &primes);
    cfx_fac_sub(&fn, &fk);
    cfx_fac_sub(&fn, &fnk);
    cfx_vec_free(&primes);
    return fn; /* sorted, coalesced */
}
