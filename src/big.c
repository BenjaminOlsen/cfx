#include "cfx/big.h"
#include "cfx/fac.h"
#include "cfx/algo.h"
#include "cfx/types.h"
#include "cfx/macros.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>

/* accumulator type of X bits and max value of multiplicands to fit therein */
#if defined(__SIZEOF_INT128__)
#define acc_t uint128_t
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

void cfx_big_copy(cfx_big_t* dst, const cfx_big_t* src) {
    cfx_big_init(dst);
    if (src->n) {
        cfx_big_reserve(dst, src->n);
        memcpy(dst->limb, src->limb, src->n*sizeof(uint64_t));
    }
    dst->n = src->n;
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
    // CFX_PRINT_DBG("realloc new cap: %zu\n", new_cap);
}

void cfx_big_enable_fac(cfx_big_t* b) {
    if (!b->cache) {
        b->cache = (cfx_fac_cache_t*)malloc(sizeof(struct cfx_fac_cache));
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

static inline void cfx_big_trim_leading_zeros(cfx_big_t* b) {
    while (b->n && b->limb[b->n-1] == 0) --b->n;
    if (b->n == 0 && b->cap){ b->limb[0] = 0; } /* dont trim last zero */
}

void cfx_big_set_val(cfx_big_t* b, uint64_t v) {
    cfx_big_init(b);
    cfx_big_reserve(b, 1);
    b->n = v ? 1 : 0;
    b->limb[0] = v;
    
    if (b->cache) {
        /* factor v */
    }
}

static void cfx_big_mul_sm_fast(cfx_big_t* b, uint64_t m) {
    uint128_t carry = 0;
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

int cfx_big_is_zero(const cfx_big_t* b) {
    return b->n <= 1 && b->limb[0] == 0;
}

int cfx_big_eq_sm(const cfx_big_t* b, uint64_t n) {
    return b->n == 1 && b->limb[0] == n;
}

int cfx_big_eq(const cfx_big_t* b1, const cfx_big_t* b2) {
    if (b1->n != b2->n) return 0;
    uint64_t diff = 0;
    for (size_t i = 0; i < b1->n; ++i) {
        diff |= (b1->limb[i] ^ b2->limb[i]);
    }
    return diff == 0;
}

void cfx_big_swap(cfx_big_t* a, cfx_big_t* b) {
    if (a == b) return;
    cfx_big_t tmp = *a;
    *a = *b;
    *b = tmp;
}

void cfx_big_sq(cfx_big_t* b) {
#if 1
    if (b->n == 0) return;
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    size_t n = b->n;
    size_t szout = 2*n;
    cfx_big_reserve(&tmp, szout);
    tmp.n = szout;

    for (size_t i = 0; i < n; ++i) {
        uint128_t carry = 0;
        uint128_t bi = b->limb[i];
        /* terms ij == terms ji -> iterate over half the terms but double them*/
        for (size_t j = i+1; j < n; ++j) {
            uint128_t prod = bi * b->limb[j];
            prod <<= 1; /* double */
            prod |= carry;
            prod += tmp.limb[i+j];
            tmp.limb[i+j] = (uint64_t)prod;
            carry = prod >> 64;
            // printf("doubling term i: %zu, j: %zu; prod: %llu, carry: %llu\n", i, j, prod, (uint64_t)carry);
        }
        uint128_t sq = bi*bi;
        sq += tmp.limb[2*i];
        uint64_t lo = (uint64_t)sq;
        uint128_t c2 = sq >> 64;
        // printf("squaring term i: %zu, lo: %llu, carry: %llu\n", i, lo, (uint64_t)carry);

        // propagate carry from cross terms into next limb
        uint128_t u = (uint128_t)tmp.limb[i + n] + carry + c2;
        tmp.limb[i + n] = (uint64_t)u;
        // store fixed lower
        tmp.limb[2 * i] = lo;
        // carry continues if u overflowed 64 bits
        uint128_t k = u >> 64;
        if (k) {
            size_t idx = i + n + 1;
            while (k) {
                printf("carry continues %llu\n", (uint64_t)k);
                if (idx >= szout) break;
                uint128_t w = (uint128_t)tmp.limb[idx] + (uint64_t)k;
                tmp.limb[idx] = (uint64_t)w;
                k = w >> 64;
                ++idx;
            }
        }
    }
    cfx_big_trim_leading_zeros(&tmp);
    cfx_big_swap(&tmp, b);
    cfx_big_free(&tmp);
#else
    if (b->n == 0) return;

    size_t n = b->n;
    size_t out_n = 2 * n;
    // temp buffer initialized to zero
    uint64_t* tmp = (uint64_t*)calloc(out_n, sizeof(uint64_t));
    if (!tmp) { /* OOM handling as per your policy */ return; }

    for (size_t i = 0; i < n; ++i) {
        uint128_t carry = 0;

        // cross terms j > i counted twice
        size_t j = i + 1;
        for (; j < n; ++j) {
            uint128_t t = (uint128_t)b->limb[i] * (uint128_t)b->limb[j];
            t <<= 1; // times 2
            t += tmp[i + j];
            t += carry;
            tmp[i + j] = (uint64_t)t;
            carry = t >> 64;
        }

        // add the square term (i,i)
        uint128_t t2 = (uint128_t)b->limb[i] * (uint128_t)b->limb[i];
        t2 += tmp[2 * i];
        uint64_t lo = (uint64_t)t2;
        uint128_t c2 = t2 >> 64;

        // propagate carry from cross terms into next limb
        uint128_t u = (uint128_t)tmp[i + n] + carry + c2;
        tmp[i + n] = (uint64_t)u;
        // store fixed lower
        tmp[2 * i] = lo;
        // carry continues if u overflowed 64 bits
        uint128_t k = u >> 64;
        if (k) {
            size_t idx = i + n + 1;
            while (k) {
                if (idx >= out_n) break;
                uint128_t w = (uint128_t)tmp[idx] + (uint64_t)k;
                tmp[idx] = (uint64_t)w;
                k = w >> 64;
                ++idx;
            }
        }
    }

    // move tmp into b
    // reserve and copy
    cfx_big_reserve(b, out_n);
    memcpy(b->limb, tmp, out_n * sizeof(uint64_t));
    b->n = out_n;
    cfx_big_trim_leading_zeros(b);
    free(tmp);
#endif
}

void cfx_big_mul(cfx_big_t* b, const cfx_big_t* m) {
    if(cfx_big_is_zero(m)) {
        printf("multiplying b by zero!\n");
        cfx_big_free(b);
        cfx_big_set_val(b, 0);
        return;
    }
    if (cfx_big_eq(b, m)) {
        cfx_big_sq(b);
        return;
    }

    cfx_big_t tmp;
    cfx_big_init(&tmp);
    size_t szb = b->n;
    size_t szm = m->n;

    cfx_big_reserve(&tmp, szb + szm);
    tmp.n = szm + szb;

    for (size_t i = 0; i < szm; ++i) {
        uint128_t carry = 0;
        uint128_t mi = (uint128_t)m->limb[i];
        for (size_t j = 0; j < szb; ++j) {
            uint128_t prod = mi * b->limb[j];
            prod += carry;
            tmp.limb[i+j] += (uint64_t)prod;
            carry = prod >> 64;
            // printf("i: %zu, j: %zu; prod: %llu, carry: %llu\n", i, j, (uint64_t)prod, (uint64_t)carry);
        }
        uint128_t l = carry + tmp.limb[i+szb];
        if (carry) {
            tmp.limb[i+szb] = (uint64_t)l;
        }
    }
    cfx_big_trim_leading_zeros(&tmp);
    cfx_big_swap(&tmp, b);
    cfx_big_free(&tmp);
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
        uint128_t s = (uint128_t)b->limb[i] + n; // sum
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
    cfx_big_trim_leading_zeros(b);
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
    uint128_t acc = 0;
    for (size_t i = n->n; i-- > 0;) {
        acc = ((acc << 64) + n->limb[i] ) % p;
    }
    return (uint64_t)acc;
}

void cfx_big_from_limbs(cfx_big_t* b, const uint64_t* limbs, size_t n) {
    cfx_big_init(b);
    cfx_big_reserve(b, n);
    memcpy(b->limb, limbs, n * sizeof(uint64_t));
    b->n = n;
}

/* Materialize factorization into cfx_big_t */
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t *f) {
    cfx_big_set_val(b, 1);
    for (size_t i = 0; i < f->len; i++){
        cfx_big_powmul_prime(b, f->data[i].p, f->data[i].e);
    }
}


int cfx_big_div(cfx_big_t* b, const cfx_big_t* d, cfx_big_t* r) {

}

void cfx_big_to_fac(cfx_fac_t *f, const cfx_big_t *b) {
    
}

uint64_t cfx_big_div_sm(cfx_big_t* b, uint64_t d) {
    uint128_t rem = 0;
    for (ssize_t i = (ssize_t)b->n - 1; i >=0; --i) {
        uint128_t cur = ((uint128_t)rem << 64) | b->limb[i];
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
        uint128_t cur = ((uint128_t)rem << 64) | b->limb[i];
        uint64_t q = (uint64_t)(cur / d);   // 128/32 -> 128/64 promoted; OK
        rem = (uint64_t)(cur % d);
        b->limb[i] = q;
    }
    cfx_big_trim_leading_zeros(b);
    return (uint32_t)rem;
}

/* uses Horner's rule to evaluate the polynomial with x = B=2^64 */
/* P(x)=(...((an​x+an−1​)x+an−2​)x+...+a1​)x+a0 */
uint64_t cfx_big_mod_sm(cfx_big_t* b, uint64_t m) {
    if (b->n == 0) return 0;
    uint128_t acc = 0;
    for (ssize_t i = (ssize_t)b->n - 1; i >= 0; --i) {
        /* acc << 64 is acc * 2^64, which is (a_n * B); limb[i] is a_(n-1)*/
        acc = ((acc << 64) + b->limb[i]) % m;
    }
    return (uint64_t)acc;
}

int cfx_fac_from_big(cfx_fac_t* fac, const cfx_big_t* in) {
#if 0
    cfx_big_t N = *in; // make a working copy
    cfx_fac_init(fac);

    // 1) strip small primes up to B (say, 10k–100k)
    cfx_big_strip_small(&N, fac, /*B=*/100000);

    // stack of pending cofactors to split (skip 1)
    cfx_big_t stk[MAX_STACK]; size_t top=0;
    if (!cfx_big_is_one(&N)) stk[top++] = N;

    while (top) {
        cfx_big_t m = stk[--top];
        if (cfx_big_is_one(&m)) continue;

        uint64_t u64;
        if (cfx_big_fits_u64(&m, &u64)) {
            cfx_fac_from_u64(fac, u64);
            continue;
        }

        if (cfx_big_is_probable_prime(&m)) { // big-PRP (MR or BPSW)
            // emit m as a prime factor (store big primes in fac as u64? or separate big-prime path)
            // If you require u64 primes only, keep splitting until fit; otherwise extend fac to hold big p.
            cfx_fac_push_big(fac, &m, 1);
            continue;
        }

        cfx_big_t d; cfx_big_init(&d);
        if (!cfx_big_rho_brent_split(&m, &d)) {
            // optional: try p-1; else increase small-prime bound and repeat
            // as a fallback, treat as prime (to avoid infinite loops)
            cfx_fac_push_big(fac, &m, 1);
            continue;
        }

        cfx_big_t q;
        cfx_big_div_big(&q, &m, &d); // q = m / d (exact)
        stk[top++] = d;
        stk[top++] = q;
    }

    cfx_fac_coalesce_sort(fac);
    return 0;
#endif
}


/* Convert cfx_big_t to decimal string */
char* cfx_big_to_str(const cfx_big_t* src, size_t *sz_out) {
    if (src->n == 0) {
        char* s = (char*)malloc(2);
        s[0]='0';
        s[1]='\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    cfx_big_t tmp = *src;
    tmp.cap = tmp.n;
    tmp.limb = (uint64_t*)malloc(tmp.n * sizeof(uint64_t));
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
    char* s = (char*)malloc(len+1);
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
