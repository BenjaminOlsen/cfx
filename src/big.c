#include "cfx/big.h"
#include "cfx/factorization.h"
#include <string.h>
#include <stdlib.h>

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
    b->limb = (uint32_t*)realloc(b->limb, new_cap*sizeof(uint32_t));
    b->cap = new_cap;
}

void cfx_big_set_u32(cfx_big_t* b, uint32_t v) {
    cfx_big_reserve(b, 1);
    b->n = v ? 1:0;
    if (v) b->limb[0] = v;
}

void cfx_big_mul_u32(cfx_big_t* b, uint32_t m) {
    if (m == 1) return;
    if (b->n == 0) {
        cfx_big_set_u32(b, 0);
        return;
    }
    uint64_t carry = 0;
    for (size_t i = 0; i < b->n; ++i) {
        uint64_t t = b->limb[i] * m + carry;
        b->limb[i] = (uint32_t)(t % CFX_BIG_BASE);
        carry = t / CFX_BIG_BASE;
    }
    while (carry) {
        // have to add a digit 
        cfx_big_reserve(b, b->n+1);
        b->limb[b->n] = (uint32_t)(carry % CFX_BIG_BASE);
        b->n++;
        carry /= CFX_BIG_BASE;
    }

}

/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
void cfx_big_powmul_prime(cfx_big_t *b, uint32_t p, uint32_t e){
    /* Simple robust version: multiply by p, e times, but group in powers that fit in 32-bit. */
    /* Find largest t so that p^t fits in uint32_t */
    uint32_t t = 1;
    uint64_t acc = p;
    while ((acc*acc) <= 0xFFFFFFFFull) { 
        acc *= acc;
        t *= 2;
    }

    /* Now multiply by (p^t)^(e/t) and then the remainder */
    /* Compute p^t */
    /* compute power safely */
    uint32_t pow_t = p;
    uint32_t rem_t = t;
    /* fast power to get pow_t = p^t */
    pow_t = 1;
    acc = p;
    uint32_t tt = t;
    while (tt) {
        if (tt & 1u) {
            pow_t = (uint32_t)(pow_t * acc); 
        }
        tt >>= 1u;
        if (tt) {
            acc = acc*acc;
        }
    }
    uint32_t q = e / t;
    uint32_t r = e % t;

    for (uint32_t i = 0; i < q; i++) cfx_big_mul_u32(b, pow_t);

    /* multiply remainder by binary exponentiation (still fits u32) */
    uint32_t rempow = 1;
    acc = p;
    uint32_t rr = r;

    while (rr) {
        if (rr & 1u) rempow = (uint32_t)(rempow * acc);
        rr >>= 1u;
        if (rr) acc = acc*acc;
    }
    if (rempow != 1) cfx_big_mul_u32(b, rempow);
}

/* Materialize factorization into cfx_big_t */
void cfx_big_from_factorization(cfx_big_t *b, const cfx_factorization_t *f){
    cfx_big_set_u32(b, 1);
    for (size_t i=0;i<f->len;i++){
        cfx_big_powmul_prime(b, f->data[i].p, f->data[i].e);
    }
}

/* Convert cfx_big_t to decimal string */
char* cfx_big_to_string(const cfx_big_t *b){
    if (b->n == 0) {
        char *s = (char*)malloc(2);
        s[0]='0';
        s[1]='\0';
        return s;
    }

    /* worst-case digits ~ 9*n + 1 */
    size_t cap = b->n*9 + 1;
    char *buf = (char*)malloc(cap+1);
    char *p = buf + cap; *p = '\0';
    // print most significant limb without zero padding
    uint32_t ms = b->limb[b->n-1];
    do { *--p = (char)('0' + (ms % 10)); ms/=10; } while (ms);
    // remaining limbs with 9-digit zero padding
    for (ssize_t i=(ssize_t)b->n-2; i>=0; --i){
        uint32_t x = b->limb[i];
        for (int k=0;k<9;k++){ *--p = (char)('0' + (x % 10)); x/=10; }
    }
    size_t len = buf+cap - p;
    char *out = (char*)malloc(len+1);
    memcpy(out, p, len+1);
    free(buf);
    return out;
}