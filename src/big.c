#include "cfx/big.h"
#include "cfx/fac.h"
#include "cfx/types.h"
#include "cfx/macros.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>

/* accumulator type of X bits and max value of multiplicands to fit therein */
#if defined(__SIZEOF_INT128__)
#define acc_t __uint128_t
#define MAX_ACC_MUL 0xFFFFFFFFFFFFFFFFllu   /* floor(sqrt(UINT128_MAX+1)) - 1 */
#else
#define acc_t uint64_t
#define MAX_ACC_MUL 0xFFFFFFFFllu           /* floor(sqrt(UINT64_MAX+1)) - 1 */
#endif

void cfx_big_init(cfx_big_t* b) {
    b->limb = NULL;
    b->cache = NULL;
    b->n = 0;
    b->cap = 0;
}

void cfx_big_free(cfx_big_t* b) {
    free(b->limb);
    free(b->cache);
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

void cfx_big_enable_fac(cfx_big_t* b) {
    if (!b->cache) {
        b->cache = malloc(sizeof(struct cfx_fac_cache));
        cfx_fac_init(&b->cache->primes);
        cfx_big_init(&b->cache->cofactor);
        b->cache->state = CFX_FAC_NONE;
        b->cache->bnd = 0;
    }
}

void cfx_big_disable_fac(cfx_big_t* b) {
    if (b->cache) {
        cfx_fac_free(&b->cache->primes);
        cfx_big_free(&b->cache->cofactor);
        free(b->cache);
        b->cache = NULL;
    }
}

static inline void _cfx_big_trim_leading_zeros(cfx_big_t* b) {
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
    
    if (b->cache) {
        /* factor v */
    }
}


    
static void cfx_big_mul_sm_fast(cfx_big_t* b, uint64_t m) {
    __uint128_t carry = 0;
    acc_t t;
    for (size_t i = 0; i < b->n; ++i) {
        t = b->limb[i];
        t *= m;
        t += carry;
        b->limb[i] = (uint64_t)t;
        carry = t >> 64;
    }
    while (carry) {
        // have to add a digit 
        cfx_big_reserve(b, b->n+1);
        b->limb[b->n] = (uint64_t)carry;
        b->n++;
        carry >>= 64;
    }
}

/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
void cfx_big_powmul_prime(cfx_big_t* b, uint64_t p, uint64_t e) {
    
    /* Find largest t so that p^t fits in 32 bits -> p^2t fits in 64 */
    uint64_t t = 1;
    acc_t acc = p;
    const acc_t lim = MAX_ACC_MUL;
    
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
        __uint128_t s = (__uint128_t)b->limb[i] + n; // sum
        b->limb[i] = (uint64_t)s;
        n = (uint64_t)(s >> 64);
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
        }
    }
    b->limb[0] -= n;
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

uint64_t cfx_big_mod_small(const cfx_big_t* n, uint64_t p) {
    __uint128_t acc = 0;
    for (size_t i = n->n; i-- > 0;) {
        acc = ( (acc << 64) + n->limb[i] ) % p;
    }
    return (uint64_t)acc;
}


/* Materialize factorization into cfx_big_t */
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t *f) {
    cfx_big_set_val(b, 1);
    for (size_t i = 0; i < f->len; i++){
        cfx_big_powmul_prime(b, f->data[i].p, f->data[i].e);
    }
}

void cfx_big_to_fac(cfx_fac_t *f, const cfx_big_t *b) {
    
}

uint64_t cfx_big_div_sm(cfx_big_t* b, uint64_t d) {
    __uint128_t rem = 0;
    for (ssize_t i = (ssize_t)b->n - 1; i >=0; --i) {
        __uint128_t cur = ((__uint128_t)rem << 64) | b->limb[i];
        uint64_t q = (uint64_t)(cur / d);
        rem = (uint64_t)(cur % d);
        b->limb[i] = q;
    }
    return rem;
}

// Divides x (base 2^64) by uint32_t d, returns remainder.
uint32_t cfx_big_div_sm_u32(cfx_big_t* b, uint32_t d) {
    uint64_t rem = 0;
    for (ssize_t i = (ssize_t)b->n - 1; i >= 0; --i) {
        __uint128_t cur = ((__uint128_t)rem << 64) | b->limb[i];
        uint64_t q = (uint64_t)(cur / d);   // 128/32 -> 128/64 promoted; OK
        rem = (uint64_t)(cur % d);
        b->limb[i] = q;
    }
    _cfx_big_trim_leading_zeros(b);
    return (uint32_t)rem;
}

/* uses Horner's rule to evaluate the polynomial with x = B=2^64 */
/* P(x)=(...((an​x+an−1​)x+an−2​)x+...+a1​)x+a0 */
uint64_t cfx_big_mod_sm(cfx_big_t* b, uint64_t m) {
    if (b->n == 0) return 0;
    __uint128_t acc = 0;
    for (ssize_t i = (ssize_t)b->n - 1; i >= 0; --i) {
        /* acc << 64 is acc * 2^64, which is (a_n * B); limb[i] is a_(n-1)*/
        acc = ((acc << 64) + b->limb[i]) % m;
    }
    return (uint64_t)acc;
}

/* Convert cfx_big_t to decimal string */
char* cfx_big_to_str(const cfx_big_t* src, size_t *sz_out) {
    if (src->n == 0) {
        char* s = malloc(2);
        s[0]='0';
        s[1]='\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    cfx_big_t tmp = *src;
    tmp.cap = tmp.n;
    tmp.limb = malloc(tmp.n * sizeof(uint64_t));
    memcpy(tmp.limb, src->limb, tmp.n * sizeof(uint64_t));

    enum { CHUNK_BASE = 1000000000u, CHUNK_DIGS = 9 };
    size_t maxdig = src->n * 20; /* log10(2^64) == 19.2659... */
    uint32_t *chunks = (uint32_t*)malloc(maxdig);
    size_t k = 0;

    while (tmp.n) {
        // printf("tmp.n:%zu; k: %zu/%zu\n", tmp.n, k, maxdig);
        chunks[k++] = cfx_big_div_sm_u32(&tmp, CHUNK_BASE);
    }
    cfx_big_free(&tmp);

    // build string
    // first chunk has no zero-padding, others are 9-digit padded
    size_t len = 0;
    {
        uint32_t first = chunks[k-1];
        char buf[16];
        int l = snprintf(buf, sizeof buf, "%u", first);
        len += l + (k-1) * CHUNK_DIGS;
    }
    char* s = malloc(len+1);
    char* p = s;

    // write first
    p += sprintf(p, "%u", chunks[k-1]);
    // write rest padded
    for (ssize_t i = (ssize_t)k-2; i >= 0; --i) {
        p += sprintf(p, "%09u", chunks[i]);
    }
    s[len] = '\0';
    if (sz_out) *sz_out = len;
    free(chunks);
    return s;
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
