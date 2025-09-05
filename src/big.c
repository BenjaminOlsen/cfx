#include "cfx/big.h"
#include "cfx/fac.h"
#include "cfx/types.h"
#include "cfx/macros.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>

#if defined(__SIZEOF_INT128__)
#define acc_t __uint128_t
#define MAX_ACC 0xFFFFFFFFFFFFFFFFllu
#else
#define acc_t __uint128_t
#define MAX_ACC 0xFFFFFFFFllu
#endif

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
    size_t new_cap = b->cap ? 2*b->cap : 32;
    if (new_cap < need) {
        new_cap = need;
    }
    size_t old_cap = b->cap;
    b->limb = (uint64_t*)realloc(b->limb, new_cap*sizeof(uint64_t));
    memset(b->limb + old_cap, 0, (new_cap - old_cap) * sizeof(uint64_t));
    b->cap = new_cap;
    CFX_PRINT_DBG("realloc new cap: %zu\n", new_cap);
}

static void _cfx_big_trim_leading_zeros(cfx_big_t* b) {
    if (b->n == 0) return;
    size_t i = b->n - 1;
    while (i >= 0 && b->limb[i] == 0 && b->n > 0) {
        --b->n;
        --i;
    }
}

void cfx_big_set_val(cfx_big_t* b, uint64_t v) {
    cfx_big_reserve(b, 1);
    b->n = v ? 1:0;
    if (v) b->limb[0] = v;
}


    
static void cfx_big_mul_sm_fast(cfx_big_t* b, uint64_t m) {
    uint64_t carry = 0;
    acc_t t;
    for (size_t i = 0; i < b->n; ++i) {
        t = b->limb[i];
        t *= m;
        t += carry;
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
    // _cfx_big_trim_leading_zeros(b);
}

/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
void cfx_big_powmul_prime(cfx_big_t* b, uint64_t p, uint64_t e) {
    
    /* Find largest t so that p^t fits in 32 bits -> p^2t fits in 64 */
    uint64_t t = 1;
    acc_t acc = p;
    const acc_t lim = MAX_ACC;
    
    while (acc <= lim / acc) {
        acc *= acc;
        t *= 2u;
    }
    /* now t is the largest s.t. p^t <= lim */
    CFX_PRINT_DBG("multiplying by %llu^%llu by breaking it into (%llu^%llu)^(%llu/%llu) \n", p, e, p, t, e, t);

    /* Now multiply by (p^t)^(e/t) and then the remainder */
    /* Compute p^t */
    /* compute power safely */
    uint64_t pow_t = p;
    
    /* fast power to get pow_t = p^t by using binary expansion of t:
    p^t = p^(b_i*2^i) * p^(b_(i-1)*2^(i-1) * ... */
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

    for (uint64_t i = 0; i < q; i++) cfx_big_mul_sm_fast(b, pow_t);

    /* multiply remainder by binary exponentiation (still fits u32) */
    uint64_t rempow = 1;
    acc = p;
    uint64_t rr = r;

    while (rr) {
        if (rr & 1u) rempow = (uint64_t)(rempow * acc);
        rr >>= 1u;
        acc = acc*acc;
    }
    if (rempow != 1) cfx_big_mul_sm_fast(b, rempow);
}

void cfx_big_add_sm(cfx_big_t* b, uint64_t n) {
    if (n == 0) return;
    if (b->n == 0) {
        cfx_big_init(b);
        cfx_big_set_val(b, n);
        return;
    }

    size_t i = 0;
    while (n) {
        cfx_big_reserve(b, b->n + 1); // doesn't do anything if b->n + 1 < b->cap
        uint64_t s = b->limb[i] + n; // sum
        b->limb[i] = s % CFX_BIG_BASE;
        n = s / CFX_BIG_BASE; // use n as the carry
        ++i;
    }
    if (i > b->n) b->n = i;
}

void cfx_big_sub_sm(cfx_big_t* b, uint64_t n) {
    if (n == 0 || b->n == 0) return;
    if (b->limb[0] < n) {
        if (b->n == 1) {
            /* underflow */
            return;
        } else {
            b->limb[1] -= 1;
            b->limb[0] += (CFX_BIG_BASE - n);
            b->limb[0] %= CFX_BIG_BASE;
        }
    } else {
        b->limb[0] -= n;
    }
    _cfx_big_trim_leading_zeros(b);
}

void cfx_big_mul_sm(cfx_big_t* b, uint64_t m) {
    if (m == 1) return;
    if (m == 0 || b->n == 0) {
        cfx_big_init(b);
        cfx_big_set_val(b, 0);
        return;
    }
    cfx_big_mul_sm_fast(b, m);
}

/* Materialize factorization into cfx_big_t */
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t *f) {
    cfx_big_set_val(b, 1);
    for (size_t i = 0; i < f->len; i++){
        cfx_big_powmul_prime(b, f->data[i].p, f->data[i].e);
    }
}


/* Convert cfx_big_t to decimal string */
char* cfx_big_to_str(const cfx_big_t* b, size_t* sz_out) {
    if (b->n == 0) {
        char *s = (char*)malloc(2);
        s[0]='0';
        s[1]='\0';
        if (sz_out) *sz_out = 1;
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
        //     CFX_PRINT_DBG("printed ls: %c (ord: 0x%x)\n", *pt, *pt);
        // }
    }
    
    /* most significant limb without zero padding */
    uint64_t ms = b->limb[b->n-1];
    do {
        /* debug */
        // CFX_PRINT_DBG("printed ms: %c\n", (char)('0' + (ms % 10)));
        *(--p) = (char)('0' + (ms % 10));
        ms /= 10;
    } while (ms);

    size_t len = buf + cap - p;
    memmove(buf, p, len + 1);
    if (sz_out) *sz_out = len;
    return buf;
}

int cfx_big_from_str(cfx_big_t* out, const char* s) {
    cfx_big_init(out);
    cfx_big_set_val(out, 0);

    while (isspace((unsigned char)*s)) s++;

    for (; *s; s++) {
        if (!isdigit((unsigned char)*s)) return -1;
        uint32_t digit = (uint32_t)(*s - '0');
        cfx_big_mul_sm(out, 10); // out = out * 10
        cfx_big_add_sm(out, digit); // out = out + digit
    }
    return 0;
}
