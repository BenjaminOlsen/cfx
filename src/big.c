#include "cfx/big.h"
#include "cfx/fac.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void cfx_big_init(cfx_big_t* b) {
    b->limb = NULL;
    b->n = 0;
    b->cap = 0;
}

void cfx_big_free(cfx_big_t* b) {
    free(b->limb);
    cfx_big_init(b);
}

void cfx_big_reserve(cfx_big_t* b, size_t need) {
    if (need < b->cap) return;
    size_t new_cap = b->cap ? 2*b->cap : 16;
    if (new_cap < need) {
        new_cap = need;
    }
    b->limb = (uint64_t*)realloc(b->limb, new_cap*sizeof(uint64_t));
    b->cap = new_cap;
}

void cfx_big_set_val(cfx_big_t* b, uint64_t v) {
    cfx_big_reserve(b, 1);
    b->n = v ? 1:0;
    if (v) b->limb[0] = v;
}

void cfx_big_mul(cfx_big_t* b, uint64_t m) {
    if (m == 1) return;
    if (b->n == 0) {
        cfx_big_set_val(b, 0);
        return;
    }
    uint64_t carry = 0;
    for (size_t i = 0; i < b->n; ++i) {
        uint64_t t = b->limb[i] * m + carry;
        b->limb[i] = (uint64_t)(t % CFX_BIG_BASE);
        carry = t / CFX_BIG_BASE;
    }
    while (carry) {
        // have to add a digit 
        cfx_big_reserve(b, b->n+1);
        b->limb[b->n] = (uint64_t)(carry % CFX_BIG_BASE);
        b->n++;
        carry /= CFX_BIG_BASE;
    }

}

/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
void cfx_big_powmul_prime(cfx_big_t *b, uint64_t p, uint64_t e) {
    /* Simple robust version: multiply by p, e times, but group in powers that fit in 32-bit. */
    /* Find largest t so that p^t fits in uint64_t */
    uint64_t t = 1;
    uint64_t acc = p;
    const uint64_t lim = 0xFFFFFFFFllu;
    printf("multiplying py %llu ^ %llu\n", p, e);
    while (acc <= lim / acc) {
        acc *= acc;
        t *= 2u;
    }
    printf("t: %llu\n", t); /* t is the largest s.t. p^t <= lim */

    /* Now multiply by (p^t)^(e/t) and then the remainder */
    /* Compute p^t */
    /* compute power safely */
    uint64_t pow_t = p;
    uint64_t rem_t = t;
    /* fast power to get pow_t = p^t by using binary expansion of t:
    p^t = p^(b_i*2^i) * p^(b_(i-1)*2^(i-1)*/
    pow_t = 1;
    acc = p;
    uint64_t tt = t;
    while (tt) {
        if (tt & 1u) pow_t = (uint64_t)(pow_t * acc); 
        tt >>= 1u;
        acc = acc*acc;
    }
    uint64_t q = e / t;
    uint64_t r = e % t;

    for (uint64_t i = 0; i < q; i++) cfx_big_mul(b, pow_t);

    /* multiply remainder by binary exponentiation (still fits u32) */
    uint64_t rempow = 1;
    acc = p;
    uint64_t rr = r;

    while (rr) {
        if (rr & 1u) rempow = (uint64_t)(rempow * acc);
        rr >>= 1u;
        acc = acc*acc;
    }
    if (rempow != 1) cfx_big_mul(b, rempow);
}

/* Materialize factorization into cfx_big_t */
void cfx_big_from_factorization(cfx_big_t *b, const cfx_fac_t *f){
    cfx_big_set_val(b, 1);
    for (size_t i = 0; i < f->len; i++){
        cfx_big_powmul_prime(b, f->data[i].p, f->data[i].e);
    }
}

/* Convert cfx_big_t to decimal string */
char* cfx_big_to_string(const cfx_big_t *b, size_t* sz_out){
    if (b->n == 0) {
        char *s = (char*)malloc(2);
        s[0]='0';
        s[1]='\0';
        return s;
    }

    /* optional check */
    size_t dec_per_limb = 0;
    size_t base = CFX_BIG_BASE;
    while (base /= 10) {
        ++dec_per_limb;
    }
    assert(dec_per_limb == CFX_DIGITS_PER_LIMB);

    /* worst-case digits ~ CFX_DIGITS_PER_LIMB*n + 1 */
    size_t cap = (CFX_DIGITS_PER_LIMB * b->n) + 1;
    char *buf = (char*)malloc(cap + 1);
    char *p = buf + cap;
    *p = '\0';

    /* least significant limbs with CFX_DIGITS_PER_LIMB-digit zero padding */
    for (ssize_t i = 0; i < (ssize_t)(b->n - 1); ++i) {
        uint64_t x = b->limb[i];
        char tmp[CFX_DIGITS_PER_LIMB];
        for (int k = 8; k >= 0; --k) {
            tmp[k] = (char)('0' + (x % 10));
            x /= 10;
        }
        p -= CFX_DIGITS_PER_LIMB;
        memcpy(p, tmp, CFX_DIGITS_PER_LIMB);

        /* debug */        
        // for (char* pt = p+CFX_DIGITS_PER_LIMB-1; pt >= p; --pt) {
        //     printf("printed ls: %c (ord: 0x%x)\n", *pt, *pt);
        // }
    }
    
    /* most significant limb without zero padding */
    uint64_t ms = b->limb[b->n-1];
    do {
        /* debug */
        // printf("printed ms: %c\n", (char)('0' + (ms % 10)));
        *(--p) = (char)('0' + (ms % 10));
        ms /= 10;
    } while (ms);

    size_t len = buf + cap - p;
    memmove(buf, p, len + 1);
    if (sz_out) *sz_out = len;
    return buf;
}