#include "cfx/big.h"
#include "cfx/fac.h"
#include "cfx/algo.h"
#include "cfx/types.h"
#include "cfx/macros.h"

#include <pthread.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <ctype.h>
#include <inttypes.h>

// x86-64 BMI2/ADX fast path
// compile with -march=native (or -03 -mbmi2 -madx)
#include <immintrin.h>

/* accumulator type of X bits and max value of multiplicands to fit therein */
#if defined(__SIZEOF_INT128__)
#define acc_t uint128_t
#define MAX_ACC_MUL 0xFFFFFFFFFFFFFFFFllu   /* floor(sqrt(UINT128_MAX+1)) - 1 */
#else
#define acc_t uint64_t
#define MAX_ACC_MUL 0xFFFFFFFFllu           /* floor(sqrt(UINT64_MAX+1)) - 1 */
#endif

void cfx_big_init(cfx_big_t* b) {
    memset(b, 0, sizeof(*b));
}

void cfx_big_clear(cfx_big_t* b) {
    b->n = 0;
}

void cfx_big_free(cfx_big_t* b) {
    b->n = 0;
    b->cap = 0;
    if (b->limb) free(b->limb);
    b->limb = NULL;
    if (b->cache) free(b->cache);
    b->cache = NULL;
}

int cfx_big_copy(cfx_big_t* dst, const cfx_big_t* src) {
    if (dst == src) return 0;

    cfx_big_t tmp;
    cfx_big_init(&tmp);

    if (src->n) {
        if (cfx_big_reserve(&tmp, src->n) != 0) {
            cfx_big_free(&tmp);
            return -1; // OOM
        }
        memcpy(tmp.limb, src->limb, src->n * sizeof(*src->limb));
        tmp.n = src->n;
    } else {
        tmp.n = 0;
    }
    cfx_big_swap(&tmp, dst);
    cfx_big_free(&tmp);
    return 0;
}

void cfx_big_assign(cfx_big_t* dst, const cfx_big_t* src) {
    if (dst == src) return;
    cfx_big_reserve(dst, src->n);
    memcpy(dst->limb, src->limb, src->n * sizeof(uint64_t));
    dst->n = src->n;
}

void cfx_big_move(cfx_big_t* dst, cfx_big_t* src) {
    if (dst == src) return;
    cfx_big_free(dst);
    *dst = *src;
    cfx_big_init(src);
}

int cfx_big_is_zero(const cfx_big_t* b) {
    return (b->n == 0) || (b->n == 1 && b->limb[0] == 0);
}

int cfx_big_eq_sm(const cfx_big_t* b, uint64_t n) {
    return b->n == 1 && b->limb[0] == n;
}

int cfx_big_eq(const cfx_big_t* a, const cfx_big_t* b) {
    if (a->n != b->n) return 0;
    uint64_t diff = 0;
    for (size_t i = 0; i < a->n; ++i) {
        diff |= (a->limb[i] ^ b->limb[i]);
    }
    return diff == 0;
}

/** compare two bigs - 
 * returns: 
 * -1 if a < b
 * 0 if a == b
 * 1 if a > b
 **/
int cfx_big_cmp(const cfx_big_t* a, const cfx_big_t* b) {
    if (a->n != b->n) return (a->n < b->n) ? -1 : 1;
    for (size_t i = a->n; i-- > 0; ) {
        if (a->limb[i] != b->limb[i]) return (a->limb[i] < b->limb[i]) ? -1 : 1;
    }
    return 0;
}

void cfx_big_swap(cfx_big_t* a, cfx_big_t* b) {
    if (a == b) return;
    cfx_big_t tmp = *a;
    *a = *b;
    *b = tmp;
}

/* returns 0 on success, -1 on failure (errno set) */
int cfx_big_reserve(cfx_big_t* b, size_t need) {
    if (need <= b->cap) return 0;

    size_t old_cap = b->cap;
    size_t new_cap = old_cap ? old_cap : 32;

    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) {  // doubling would overflow
            new_cap = need;            // fall back to exactly need
            break;
        }
        new_cap *= 2;
    }

    if (new_cap > SIZE_MAX / sizeof(uint64_t)) {
        // errno = ENOMEM;
        return -1;
    }

    size_t new_bytes = new_cap * sizeof(uint64_t);

    // Use a temp so we don't clobber b->limb on failure
    void* tmp = realloc(b->limb, new_bytes);
    if (!tmp) {
        // b->limb is still valid; caller can continue using old buffer
        // or handle the error.
        return -1;
    }

    b->limb = (uint64_t*)tmp;

    // Zero the newly added region only (realloc doesn't do this for you)
    if (new_cap > old_cap) {
        size_t add = new_cap - old_cap;
        memset(b->limb + old_cap, 0, add * sizeof(uint64_t));
    }

    b->cap = new_cap;
    // CFX_PRINT_DBG("realloc new cap: %zu\n", new_cap);
    return 0;
}


void cfx_big_enable_cache(cfx_big_t* b) {
    if (!b->cache) {
        b->cache = (cfx_fac_cache_t*)malloc(sizeof(struct cfx_fac_cache));
        cfx_fac_init(&b->cache->primes);
        cfx_big_init(&b->cache->cofactor);
        b->cache->state = CFX_FAC_NONE;
        b->cache->bnd = 0;
    }
}

void cfx_big_disable_cache(cfx_big_t* b) {
    if (b->cache) {
        cfx_fac_free(&b->cache->primes);
        cfx_big_free(&b->cache->cofactor);
        free(b->cache);
        b->cache = NULL;
    }
}

static inline void cfx_big_trim(cfx_big_t* b) {
    while (b->n && b->limb[b->n-1] == 0) --b->n;
    if (b->n == 0 && b->cap){ b->limb[0] = 0; } /* dont trim last zero */
}
static inline void cfx_clear_cache(cfx_big_t* b) {
    if (b->cache) {
        free(b->cache);
        b->cache = NULL;
    }
}

/* assumes b is already initted */
void cfx_big_set_val(cfx_big_t* b, uint64_t v) {
    if (v == 0) {
        b->n = 0;
        if (b->cap) b->limb[0] = 0;
        cfx_clear_cache(b);
        return;
    }

    if (b->cap < 1) {
        int rc = cfx_big_reserve(b, 1);
        if (rc) /* error */ return;
    }

    b->limb[0] = v;
    if (b->n > 1) { /* hygiene: zero out any higher limbs, if they exist */
        memset(b->limb + 1, 0, (b->n-1) * sizeof(*b->limb));
    }
    b->n = 1;
    cfx_clear_cache(b);
}


static inline void _mul_sm_fast(cfx_big_t* b, uint64_t m) {
    uint64_t* p = b->limb;
    size_t    n = b->n;
#if defined(__x86_64__) && defined(__BMI2__) || defined(__INTELLISENSE__) /* FU!!!! */
#pragma message "Compiling with instrinsics!! "

    if (b->cap <= n) {
        cfx_big_reserve(b, n + 1);
        p = b->limb;
    }

    uint64_t carry = 0;
    for (size_t i = 0; i < n; ++i) {
        uint64_t lo, hi;
        // lo = (p[i] * m) low 64, hi = high 64 — flags NOT clobbered
        lo = _mulx_u64(p[i], m, &hi);

        // add previous carry into lo, get carry-out in c1
        unsigned char c1 = _addcarry_u64(0, lo, carry, &lo);
        // hi += c1
        hi += c1;

        p[i] = lo;
        carry = hi;
    }

    p[n] = carry;
    b->n = n + (carry != 0);
#else // Portable fallback (uint128_t)
    uint128_t carry = 0;
    for (size_t i = 0; i < n; ++i) {
        uint128_t t = (uint128_t)p[i] * m + carry;
        p[i] = (uint64_t)t;
        carry = t >> 64;
    }
    p[n] = (uint64_t)carry;
    b->n = n + ((uint64_t)carry != 0);
#endif
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

    for (uint64_t i = 0; i < q; i++) _mul_sm_fast(b, pow_t);

    /* multiply remainder by binary exponentiation (still fits u32) */
    uint64_t rempow = 1;
    acc = p;
    uint64_t rr = r;

    while (rr) {
        if (rr & 1u) rempow = (uint64_t)(rempow * acc);
        rr >>= 1u;
        acc = acc*acc;
    }
    if (rempow != 1) _mul_sm_fast(b, rempow);
}



void cfx_big_sq(cfx_big_t* b) {
#if 0
    if (b->n == 0) return;
    cfx_big_t ret;
    cfx_big_init(&ret);
    size_t n = b->n;
    size_t szout = 2*n;
    cfx_big_reserve(&ret, szout);
    ret.n = szout;

    for (size_t i = 0; i < n; ++i) {
        uint128_t carry = 0;
        uint128_t bi = b->limb[i];
        /* terms ij == terms ji -> iterate over half the terms but double them*/
        for (size_t j = i+1; j < n; ++j) {
            uint128_t prod = bi * b->limb[j];
            prod <<= 1; /* double */
            prod |= carry;
            prod += ret.limb[i+j];
            ret.limb[i+j] = (uint64_t)prod;
            carry = prod >> 64;
            // printf("doubling term i: %zu, j: %zu; prod: %llu, carry: %llu\n", i, j, prod, (uint64_t)carry);
        }
        uint128_t sq = bi*bi;
        sq += ret.limb[2*i];
        uint64_t lo = (uint64_t)sq;
        uint128_t c2 = sq >> 64;
        // printf("squaring term i: %zu, lo: %llu, carry: %llu\n", i, lo, (uint64_t)carry);

        // propagate carry from cross terms into next limb
        uint128_t u = (uint128_t)ret.limb[i + n] + carry + c2;
        ret.limb[i + n] = (uint64_t)u;
        // store fixed lower
        ret.limb[2 * i] = lo;
        // carry continues if u overflowed 64 bits
        uint128_t k = u >> 64;
        if (k) {
            size_t idx = i + n + 1;
            while (k) {
                printf("carry continues %llu\n", (uint64_t)k);
                if (idx >= szout) break;
                uint128_t w = (uint128_t)ret.limb[idx] + (uint64_t)k;
                ret.limb[idx] = (uint64_t)w;
                k = w >> 64;
                ++idx;
            }
        }
    }
    cfx_big_trim(&ret);
    cfx_big_swap(&ret, b);
    cfx_big_free(&ret);
    
#elif 0
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

    // move tmp into ret
    cfx_big_t ret;
    cfx_big_reserve(ret, out_n);
    memcpy(ret->limb, tmp, out_n * sizeof(uint64_t));
    ret->n = out_n;
    cfx_big_trim(ret);
    free(tmp);
    return ret;
#else 
    const size_t n = b->n;
    cfx_big_t ret;
    cfx_big_init(&ret);

    if (n == 0) {
        return;
    }

    size_t cnt = 0;
    
    cfx_big_reserve(&ret, 2*n);
    memset(ret.limb, 0, 2*n * sizeof(uint64_t));
    ret.n = 2*n;

    // 1) Cross terms: for i < j, add 2*b[i]*b[j] into ret[i+j]
    for (size_t i = 0; i < n; ++i) {
        uint128_t carry = 0;
        for (size_t j = i + 1; j < n; ++j) {
            uint128_t p = (uint128_t)b->limb[i] * b->limb[j];

            // add p once
            uint128_t t = (uint128_t)ret.limb[i + j]
                            + (uint64_t)p
                            + (uint64_t)carry;
            ret.limb[i + j] = (uint64_t)t;
            carry = (carry >> 64) + (t >> 64) + (p >> 64);

            // add p again (to double) -- same carry rule
            t = (uint128_t)ret.limb[i + j]
                + (uint64_t)p
                + (uint64_t)carry;
            ret.limb[i + j] = (uint64_t)t;
            carry = (carry >> 64) + (t >> 64) + (p >> 64);
            cnt++;
        }

        // propagate whatever is left in 'carry'
        size_t k = i + n; // next column after the last updated (i + (n-1))
        while (carry) {
            uint128_t t = (uint128_t)ret.limb[k] + (uint64_t)carry;
            ret.limb[k] = (uint64_t)t;
            carry = (carry >> 64) + (t >> 64);
            ++k;
            cnt++;
        }
    }

    // 2) Diagonals: add b[i]^2 once at ret[2*i]
    for (size_t i = 0; i < n; ++i) {
        uint128_t sq = (uint128_t)b->limb[i] * b->limb[i];

        uint128_t t = (uint128_t)ret.limb[2*i] + (uint64_t)sq;
        ret.limb[2*i] = (uint64_t)t;
        uint128_t c = (t >> 64) + (sq >> 64);

        size_t k = 2*i + 1;
        cnt++;
        while (c) {
            t = (uint128_t)ret.limb[k] + (uint64_t)c;
            ret.limb[k] = (uint64_t)t;
            c = (c >> 64) + (t >> 64);
            ++k;
        }
    }

    // trim
    while (ret.n && ret.limb[ret.n - 1] == 0) --ret.n;

    // printf("%s loop cnt: %zu\n", __func__, cnt);
    cfx_big_swap(&ret, b);
    cfx_big_free(&ret);
#endif
}

/* TODO */
void cfx_big_mul_fft(cfx_big_t* b, const cfx_big_t* m) {
    cfx_big_mul(b, m);
}

void cfx_big_mul_csa(cfx_big_t* b, const cfx_big_t* m) {
    if (cfx_big_is_zero(b) || cfx_big_is_zero(m)) {
        cfx_big_set_val(b, 0);
        return;
    }

    const size_t nb = b->n;
    const size_t nm = m->n;
    const size_t nout = nm + nb;

    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_reserve(&tmp, nout);

    cfx_mul_csa_portable(b->limb, nb, m->limb, nm, tmp.limb);
    tmp.n = nout;
    cfx_big_trim(&tmp);
    cfx_big_swap(&tmp, b);
    cfx_big_free(&tmp);
}

/* assumes scratch is allocated with the appropriate size b->n + m->n already. */
void cfx_big_mul_csa_scratch(cfx_big_t* b, const cfx_big_t* m, cfx_mul_scratch_t* scratch) {
    const size_t nb = b->n;
    const size_t nm = m->n;
    size_t nout = nm + nb;
    // cfx_mul_scratch_alloc(scratch, nout);
    cfx_mul_scratch_zero(scratch, nout);
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_reserve(&tmp, nout);
    cfx_mul_csa_portable_fast(b->limb, nb, m->limb, nm, tmp.limb, scratch);
    tmp.n = nout;
    cfx_big_trim(&tmp);
    cfx_big_swap(&tmp, b);
    cfx_big_free(&tmp);
}

void cfx_big_mul(cfx_big_t* b, const cfx_big_t* m) {
    
    if (cfx_big_is_zero(b) || cfx_big_is_zero(m)) {
        cfx_big_set_val(b, 0);
        return;
    }
    /* not sure this is worth it:*/
    /* if (cfx_big_eq(b, m)) {
        cfx_big_sq(b);
        return;
    }
    */

    size_t nb = b->n;
    size_t nm = m->n;

    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_reserve(&tmp, nb + nm);
    tmp.n = nb + nm;

    for (size_t i = 0; i < nm; ++i) {
        uint128_t carry = 0;
        uint128_t mi = (uint128_t)m->limb[i];

        size_t k = i; // index into tmp
        for (size_t j = 0; j < nb; ++j, ++k) {
            uint128_t s = (uint128_t)tmp.limb[k]
                        + mi * (uint128_t)b->limb[j]
                        + carry;
            tmp.limb[k] = (uint64_t)s;
            carry = s >> 64;
        }

        // propagate any remaining carry
        while (carry) {
            uint128_t s = (uint128_t)tmp.limb[k] + carry;
            tmp.limb[k] = (uint64_t)s;
            carry = s >> 64;
            ++k; // safe: nb+nm is enough for two nb/nm-digit numbers
        }
    }

    cfx_big_trim(&tmp);
    cfx_big_swap(&tmp, b);
    cfx_big_free(&tmp);
}

/* b += a */
void cfx_big_add(cfx_big_t* b, const cfx_big_t* a) {
    uint64_t carry = 0;
    size_t i = 0;

    while (i < a->n || carry) {
        cfx_big_reserve(b, i + 1);

        uint128_t bi = (i < b->n) ? b->limb[i] : 0;
        uint128_t ai = (i < a->n) ? a->limb[i] : 0;
        uint128_t s  = bi + ai + carry;

        if (i >= b->n) b->n = i + 1;      // we’re extending b
        b->limb[i] = (uint64_t)s;
        carry      = (uint64_t)(s >> 64);
        ++i;
    }

    // If a->n > b->n and carry ended at 0, we may still need to bump b->n
    if (i > b->n) b->n = i;
}


void cfx_big_add_sm(cfx_big_t* b, uint64_t n) {
    if (n == 0) return;
    if (b->n == 0) {
        cfx_big_set_val(b, n);
        return;
    }

    size_t i = 0;
    while (i < b->n) {
        uint128_t s = (uint128_t)b->limb[i] + n;
        b->limb[i] = (uint64_t)s;
        n = (uint64_t)(s >> 64);
        if (n == 0) return;
        ++i;
    }

    // carry left; append one limb
    cfx_big_reserve(b, b->n + 1);
    b->limb[b->n++] = n;      // n < 2^64 here
}

void cfx_big_sub(cfx_big_t* a, const cfx_big_t* b) {
    /* assumes a >= b; subtract b from a */
    __uint128_t borrow = 0;
    size_t i = 0, n = b->n;
    for (; i < n; ++i) {
        __uint128_t ai = a->limb[i];
        __uint128_t s  = (__uint128_t)b->limb[i] + borrow;
        __uint128_t r  = ai - s;
        a->limb[i] = (uint64_t)r;
        borrow = (ai < s);
    }
    while (borrow && i < a->n) {
        uint64_t ai = a->limb[i];
        uint64_t r  = ai - 1u;
        a->limb[i++] = r;
        borrow = (ai == 0u);
    }
    cfx_big_trim(a);
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
    cfx_big_trim(b);
}


void cfx_big_mul_sm(cfx_big_t* b, uint64_t m) {
    if (m == 1) return;
    if (m == 0 || b->n == 0) {
        cfx_big_init(b);
        cfx_big_set_val(b, 0);
        return;
    }
    _mul_sm_fast(b, m);
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
    cfx_big_trim(b);
}

/* Materialize factorization into cfx_big_t */
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t *f) {
    cfx_big_set_val(b, 1);
    for (size_t i = 0; i < f->len; i++){
        cfx_big_powmul_prime(b, f->data[i].p, f->data[i].e);
    }
}


void cfx_big_to_fac(cfx_fac_t *f, const cfx_big_t *b) {
    (void)f;
    (void)b;
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
    cfx_big_trim(b);
    return (uint32_t)rem;
}

/* uses Horner's rule to evaluate the polynomial with x = B=2^64 */
/* P(x)=(...((an​x+an−1​)x+an−2​)x+...+a1​)x+a0 as a running sum:
acc0 = (an*B + an-1)
acc1 = acc0*B + an-2 = (an*B + an-1)B + an-2 = anB^2 + an-1B + an-2
acc2 = acc1*B + an-3 .. etc

and each step, take % m because the remainder of the sum div m 
is the sum of the remainders modulo m
*/
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
    (void)fac;
    (void)in;
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
#endif
    return 0;
}

static inline size_t hex_digits_u64(uint64_t v) {
    if (!v) return 1;
#if defined(__GNUC__) || defined(__clang__)
    unsigned lead = (unsigned)__builtin_clzll(v);
    unsigned bits = 64u - lead;
    return (bits + 3u) / 4u; // ceil(bits/4)
#else
    size_t d = 0;
    while (v) { v >>= 4; ++d; }
    return d;
#endif
}
static size_t bits_in_u64(uint64_t x) {
    if (!x) return 0;
#if defined(__GNUC__) || defined(__clang__)
    return 64u - (size_t)__builtin_clzll(x);
#else
    size_t n = 0;
    while (x) { ++n; x >>= 1; }
    return n;
#endif
}

char* cfx_big_to_bin(const cfx_big_t* src, size_t* sz_out) {
    if (!src || src->n == 0) {
        char* s = (char*)malloc(2);
        if (!s) return NULL;
        s[0] = '0'; s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    const size_t ms_idx  = src->n - 1;              /* most-significant limb index */
    const uint64_t msval  = src->limb[ms_idx];
    const size_t ms_bits  = bits_in_u64(msval);      /* 1..64 */
    const size_t total_len = ms_bits + 64u * ms_idx; /* total bits as characters */

    char* s = (char*)malloc(total_len + 1);
    if (!s) return NULL;

    size_t pos = 0;
    const char bch[2] = {'0', '1'};

    for (size_t b = ms_bits; b-- > 0; ) {
        s[pos++] = bch[(size_t)((msval >> b) & 1u)];
    }


    for (size_t i = ms_idx; i-- > 0; ) {
        uint64_t limb = src->limb[i];
        for (size_t b = 64; b-- > 0; ) {
            s[pos++] = bch[(size_t)((limb >> b) & 1u)];
        }
    }

    s[total_len] = '\0';
    if (sz_out) *sz_out = total_len;
    return s;
}

char* cfx_big_to_hex(const cfx_big_t* src, size_t* sz_out) {
    // Treat empty/zero as "0"
    if (!src || src->n == 0) {
        char* s = (char*)malloc(2);
        if (!s) return NULL;
        s[0] = '0'; s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    /* trim leading zeros */
    size_t ms = src->n;

    const uint64_t ms_val = src->limb[ms - 1];
    const size_t ms_digits = hex_digits_u64(ms_val);
    const size_t total_len = ms_digits + (ms - 1) * 16;  // 16 hex chars per remaining limb

    char* s = (char*)malloc(total_len + 1); // +1 for NUL
    if (!s) return NULL;

    char* p = s;
    size_t rem = total_len + 1;

    // Most-significant limb without leading zeros
    int written = snprintf(p, rem, "%" PRIx64, (uint64_t)ms_val);
    assert(written > 0 && (size_t)written == ms_digits);
    p   += written;
    rem -= (size_t)written;

    // Remaining limbs, zero-padded to 16 hex chars each
    size_t k = 0;
    size_t cnt = 0;
    for (ssize_t i = (ssize_t)ms - 2; i >= 0; --i) {
        written = snprintf(p, rem, "%016" PRIx64, (uint64_t)src->limb[i]);
        assert(written == 16);
        p   += written;
        cnt += written;
        rem -= (size_t)written;

        if (!(cnt % 7)) {  /* some number.. */
            const char spinner[] = "|/-\\";
            ++k;
            printf("%zu hex digits done... %zu/%zu limbs remain... %c        \r",
                cnt, rem, total_len, spinner[k % 4]);
            fflush(stdout);
        }
    }

    // `snprintf` already wrote the final '\0' on the last call
    if (sz_out) *sz_out = total_len;
    return s;
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

    // printf("[cfx_big_to_str] building base %u representation... max digits: %zu\n",CHUNK_BASE, maxdig);
    int cnt = 0;
    const size_t n0 = tmp.n;
    while (tmp.n) {
        
        if (!(tmp.n % 100)) { 
            const char spinner[] = "|/-\\";
            printf("%zu decimal digits done... %zu/%zu limbs remain... %c        \r",
                k*9, tmp.n, n0, spinner[cnt++ % 4]);
            fflush(stdout);
        }
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

// static inline int hexval(unsigned char c) {
//     unsigned v = c - '0';
//     if (v <= 9) return (int)v;
//     v = (c | 32) - 'a'; // to lowercase trick, then to [0,5]: 'A' | 32 -> 'a'; 'a' | 32 -> 'a'
//     if (v <= 5) return (int)(v + 10);
//     return -1;                   // not a hex digit
// }

static const int8_t hex_table[256] = {
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, -1,-1,-1,-1,-1,-1,
    -1,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,10,11,12,13,14,15,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1
};

int cfx_big_from_hex(cfx_big_t* out, const char* s) {
    cfx_big_init(out);
    cfx_big_set_val(out, 0);

    while (*s && isspace((unsigned char)*s)) ++s;

    if (s[0]=='0' && (s[1]=='x' || s[1]=='X')) s += 2;

    int any = 0;
    for (; *s; ++s) {
        int d = hex_table[(unsigned char)*s];   // -1 if not hex
        if (d < 0) break;                       // stop on first non-hex
        // out = (out << 4) + d;
        cfx_big_shl_bits(out, out, 4);               // or cfx_big_mul_sm(out, 16)
        cfx_big_add_sm(out, (uint64_t)d);
        any = 1;
    }

    // allow trailing spaces only
    while (*s && isspace((unsigned char)*s)) ++s;

    if (!any || *s != '\0') return -1;          // no digits or junk trailing
    return 0;
}

int cfx_big_from_str(cfx_big_t* out, const char* s) {
    cfx_big_init(out);
    cfx_big_set_val(out, 0);

    while (isspace((unsigned char)*s)) s++;

    for (; *s; ++s) {
        if (!isdigit((unsigned char)*s)) return -1;
        uint32_t digit = (uint32_t)(*s - '0');
        cfx_big_mul_sm(out, 10); // out = out * 10
        cfx_big_add_sm(out, digit); // out = out + digit
    }
    return 0;
}

/* ------------------------------------------------------------- */
static inline unsigned _clz64(uint64_t x) { return x ? __builtin_clzll(x) : 64u; }

/* b <<= s*/
void cfx_big_shl_bits_eq(cfx_big_t* b, unsigned s) {
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_shl_bits(&tmp, b, s);
    cfx_big_swap(b, &tmp);
    cfx_big_free(&tmp);
}

/* b >>= s */
void cfx_big_shr_bits_eq(cfx_big_t* b, unsigned s) {
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_shr_bits(&tmp, b, s);
    cfx_big_swap(b, &tmp);
    cfx_big_free(&tmp);
}

/* out = x << s (0..63)  */
void cfx_big_shl_bits(cfx_big_t* out, const cfx_big_t* x, unsigned s) {
    cfx_big_init(out);
    if (x->n == 0 || s == 0) { cfx_big_copy(out, x); return; }

    cfx_big_reserve(out, x->n + 1);
    uint64_t carry = 0;
    for (size_t i = 0; i < x->n; ++i) {
        uint128_t t = ((uint128_t)x->limb[i] << s) | carry;
        out->limb[i] = (uint64_t)t;
        carry = (uint64_t)(t >> 64);
    }
    out->n = x->n;
    if (carry) { out->limb[out->n++] = carry; }
}

/* out = x >> s (0..63)  */
void cfx_big_shr_bits(cfx_big_t* out, const cfx_big_t* x, unsigned s) {
    cfx_big_init(out);
    if (x->n == 0 || s == 0) { cfx_big_copy(out, x); return; }

    cfx_big_reserve(out, x->n);
    uint64_t carry = 0;
    for (ssize_t i = (ssize_t)x->n - 1; i >= 0; --i) {
        uint128_t t = ((uint128_t)carry << 64) | x->limb[i];
        out->limb[i] = (uint64_t)(t >> s);
        carry = (uint64_t)(x->limb[i] << (64 - s));
    }
    out->n = x->n;
    cfx_big_trim(out);
}

/* Core: q = n / d; r = n % d; any of q or r may be NULL. Returns 0, or -1 if d==0. */
int cfx_big_divrem(cfx_big_t* q, cfx_big_t* r,
                   const cfx_big_t* n, const cfx_big_t* d) {

    if (cfx_big_is_zero(d)) return -1;

    /* n == 0 */
    if (cfx_big_is_zero(n)) {
        cfx_big_set_val(q, 0);
        cfx_big_set_val(r, 0);
        return 0;
    }

    /* n < d */
    if (cfx_big_cmp(n, d) < 0) {
        cfx_big_set_val(q, 0);
        cfx_big_copy(r, n);
        return 0;
    }

    /* d is single-limb -> fast path */
    if (d->n == 1) {
        uint64_t div = d->limb[0];
        uint128_t rem = 0;
        { cfx_big_reserve(q, n->n); q->n = n->n; }
        for (ssize_t i = (ssize_t)n->n - 1; i >= 0; --i) {
            uint128_t cur = (rem << 64) | n->limb[i];
            uint64_t qi = (uint64_t)(cur / div);
            rem = cur % div;
            q->limb[i] = qi;
        }
        cfx_big_trim(q);
        cfx_big_set_val(r, (uint64_t)rem);
        return 0;
    }

    /* Knuth Algorithm D (normalized) */

    /* Normalize: shift so msb of d's top limb is set */
    unsigned s = _clz64(d->limb[d->n - 1]);
    cfx_big_t V, U;
    cfx_big_shl_bits(&V, d, s);
    cfx_big_shl_bits(&U, n, s);

    /* Ensure U has one extra top limb for convenience */
    cfx_big_reserve(&U, U.n + 1);
    U.limb[U.n] = 0;
    U.n += 1;

    size_t m = V.n;           /* divisor length after norm */
    size_t nn = U.n;          /* dividend length after norm (incl extra top) */
    size_t qlen = (nn >= m) ? (nn - m) : 0; /* expected quotient length */

    cfx_big_t QT;

    cfx_big_init(&QT);
    cfx_big_reserve(&QT, qlen);
    QT.n = qlen;
    for (size_t i = 0; i < qlen; ++i) QT.limb[i] = 0;

    const uint64_t v1 = V.limb[m - 1];
    const uint64_t v2 = V.limb[m - 2]; /* safe because m>=2 in this branch */

    for (ssize_t j = (ssize_t)qlen - 1; j >= 0; --j) {
        /* Estimate qhat from the top 2 (or 3) limbs */
        uint128_t top = ((uint128_t)U.limb[j + m] << 64) | U.limb[j + m - 1];
        uint64_t qhat = (uint64_t)(top / v1);
        uint64_t rhat = (uint64_t)(top % v1);

        /* Adjust qhat if necessary (Knuth step D3) */
        if (qhat == UINT64_MAX || (uint128_t)qhat * v2 > (((uint128_t)rhat << 64) | U.limb[j + m - 2])) {
            qhat--;
            rhat += v1;
            if (rhat < v1 && (uint128_t)qhat * v2 > (((uint128_t)rhat << 64) | U.limb[j + m - 2])) {
                qhat--;
                rhat += v1;
            }
        }

        /* Multiply-subtract U[j..j+m] -= qhat * V[0..m-1] */
        uint64_t carry = 0;
        for (size_t i = 0; i < m; ++i) {
            uint128_t prod = (uint128_t)qhat * V.limb[i] + carry;
            uint64_t plo = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);

            uint64_t uj = U.limb[j + i];
            uint64_t uj_new = uj - plo;
            uint64_t borrow = (uj_new > uj); /* underflow -> borrow 1 */
            U.limb[j + i] = uj_new;
            carry += borrow;
        }
        /* subtract final carry from U[j+m] */
        uint64_t ujm = U.limb[j + m];
        uint64_t ujm_new = ujm - carry;
        uint64_t borrow_out = (ujm_new > ujm);
        U.limb[j + m] = ujm_new;

        if (borrow_out) {
            /* Too big: qhat--, add V back */
            qhat--;
            uint64_t c = 0;
            for (size_t i = 0; i < m; ++i) {
                uint128_t ssum = (uint128_t)U.limb[j + i] + V.limb[i] + c;
                U.limb[j + i] = (uint64_t)ssum;
                c = (uint64_t)(ssum >> 64);
            }
            U.limb[j + m] += c; /* cannot overflow because we just fixed an over-subtraction */
        }

        QT.limb[j] = qhat;
    }

    /* Unnormalize remainder: R = (U[0..m-1] >> s) */
    cfx_big_t Rn;
    cfx_big_init(&Rn);
    cfx_big_reserve(&Rn, m);
    Rn.n = m;
    for (size_t i = 0; i < m; ++i) Rn.limb[i] = U.limb[i];
    cfx_big_shr_bits(r, &Rn, s);
    cfx_big_free(&Rn);

    /* Output quotient (already normalized) */
    cfx_big_trim(&QT);
    cfx_big_swap(&QT, q);
    cfx_big_free(&QT);


    cfx_big_free(&U);
    cfx_big_free(&V);
    return 0;
}

/* Convenience wrappers */
int cfx_big_div_out(cfx_big_t* q, const cfx_big_t* n, const cfx_big_t* d) {
    return cfx_big_divrem(q, NULL, n, d);
}
int cfx_big_mod_out(cfx_big_t* r, const cfx_big_t* n, const cfx_big_t* d) {
    return cfx_big_divrem(NULL, r, n, d);
}

/* In-place: b := floor(b/d); optional remainder r. Alias-safe for any combination. */
int cfx_big_div_eq(cfx_big_t* b, const cfx_big_t* d, cfx_big_t* r /*nullable*/) {
    cfx_big_t qtmp, rtmp;
    cfx_big_init(&qtmp);
    cfx_big_init(&rtmp);

    int rc = cfx_big_divrem(&qtmp, r ? &rtmp : NULL, b, d);
    if (rc == 0) {
        cfx_big_swap(&qtmp, b);
        if (r) cfx_big_swap(&rtmp, r);
    }

    cfx_big_free(&qtmp);
    cfx_big_free(&rtmp);
    return rc;
}


// 128-bit + wrap counter accumulator per column.
// Value represented = hi * 2^128 + lo
typedef struct {
    uint128_t lo;
    uint64_t hi;
    uint64_t pad; // pad to 32 bytes to reduce 0 sharing during reductions
} acc128p_t;

// add 64-bit x into accumulator (as a 128-bit add)
static inline void acc_add_u64(acc128p_t* a, uint64_t x) {
    uint128_t old = a->lo;
    a->lo = old + (uint128_t)x;
    a->hi += (a->lo < old); // wrap over 2^128
}

// add 128-bit x into accumulator
// static inline void acc_add_u128(acc128p_t* a, uint128_t x) {
//     uint128_t old = a->lo;
//     a->lo = old + x;
//     a->hi += (a->lo < old);
// }

// dst[k] += src[k] for k in [0..n)
static inline void acc_vec_add(acc128p_t* dst, const acc128p_t* src, size_t n) {
    for (size_t k = 0; k < n; ++k) {
        uint128_t old = dst[k].lo;
        dst[k].lo = old + src[k].lo;
        uint64_t carry128 = (dst[k].lo < old);
        // add hi parts + carry128; staying in 64 bits is fine for "millions of limbs"
        dst[k].hi += src[k].hi + carry128;
    }
}

// Worker arguments
typedef struct {
    const uint64_t* a; size_t na;
    const uint64_t* b; size_t nb;
    size_t j_begin, j_end;          // range of rows (limbs of b) to process
    acc128p_t* local_acc;           // per-thread accumulator array (length = ncols)
    size_t ncols;                   // ncols = na + nb
} rc_worker_args_t;

static void* worker_rowblock(void* vp) {
    rc_worker_args_t* w = (rc_worker_args_t*)vp;
    const uint64_t* a = w->a;
    const uint64_t* b = w->b;
    const size_t na = w->na;
    const size_t ncols = w->ncols;

    // zero local accumulator
    memset(w->local_acc, 0, ncols * sizeof(acc128p_t));

    for (size_t j = w->j_begin; j < w->j_end; ++j) {
        uint64_t bj = b[j];
        if (!bj) continue; // skip zero rows fast

        for (size_t i = 0; i < na; ++i) {
            uint128_t p  = (uint128_t)a[i] * (uint128_t)bj;    // 64x64->128
            uint64_t lo = (uint64_t)p;
            uint64_t hi = (uint64_t)(p >> 64);
            size_t k = i + j;                   // column for lo
            // Bounds: k in [0 .. na-1 + j] <= na-1 + (nb-1) < na+nb = ncols
            acc_add_u64(&w->local_acc[k],     lo);
            acc_add_u64(&w->local_acc[k + 1], hi);
        }
    }
    return NULL;
}

// Expand acc[k] = hi*2^128 + lo into base-2^64 lanes T[] and do one global carry pass.
// out must have length >= ncols + 3.
static void expand_and_carry(const acc128p_t* acc, size_t ncols, uint64_t* out)
{
    // 1) clear T (we reuse 'out' as T)
    memset(out, 0, (ncols + 3) * sizeof(uint64_t));

    // 2) expand each column into up to 3 neighboring columns, locally handling tiny carries
    for (size_t k = 0; k < ncols; ++k) {
        uint128_t lo = acc[k].lo;
        uint64_t lo0 = (uint64_t)lo;
        uint64_t lo1 = (uint64_t)(lo >> 64);
        uint64_t hi0 = acc[k].hi; // each unit here is exactly 2^128 = 2 limbs

        // T[k] += lo0
        uint128_t t = (uint128_t)out[k] + lo0;
        out[k] = (uint64_t)t;
        uint64_t c0 = (uint64_t)(t >> 64);

        // T[k+1] += lo1 + c0
        t = (uint128_t)out[k+1] + lo1 + c0;
        out[k+1] = (uint64_t)t;
        uint64_t c1 = (uint64_t)(t >> 64);

        // T[k+2] += hi0 + c1
        t = (uint128_t)out[k+2] + hi0 + c1;
        out[k+2] = (uint64_t)t;
        uint64_t c2 = (uint64_t)(t >> 64);

        // Any leftover tiny carry bubbles one more limb; global pass will finish everything.
        out[k+3] += c2;
    }

    // 3) single left-to-right carry sweep
    uint64_t carry = 0;
    for (size_t k = 0; k < ncols + 3; ++k) {
        uint128_t t = (uint128_t)out[k] + carry;
        out[k]  = (uint64_t)t;
        carry   = (uint64_t)(t >> 64);
    }
    if (carry) {
        // ensure the caller reserved at least ncols+4 if they want to keep this
        out[ncols + 3] = carry;
    }
}

// Top-level: b *= m using row-parallel accumulation + single carry pass.
// threads<=0 → auto (online CPUs). threads is capped to nb.
void cfx_big_mul_rows_pthreads(cfx_big_t* b, const cfx_big_t* m, int threads)
{
    const size_t na = b->n;
    const size_t nb = m->n;

    if (!na || !nb) { cfx_big_set_val(b, 0); return; }
    if (nb == 1 && m->limb[0] == 1) { return; } // b *= 1

    // Copy multiplicand 'a' because we'll overwrite b.
    uint64_t* a_copy = (uint64_t*)malloc(na * sizeof(uint64_t));
    if (!a_copy) { /* handle OOM */ abort(); }
    memcpy(a_copy, b->limb, na * sizeof(uint64_t));

    // Plan threads
    long hw = sysconf(_SC_NPROCESSORS_ONLN);
    // CFX_PRINT_DBG("hw threads: %ld\n", hw);
    if (threads <= 0) threads = (hw > 0) ? (int)hw : 1;
    if (threads > (int)nb) threads = (int)nb;
    if (threads < 1) threads = 1;

    const size_t ncols = na + nb;          // columns 0..(na+nb-1)
    const size_t acc_len = ncols;          // acc arrays length
    const size_t out_len = ncols + 4;      // room for expansion + final carry

    // Allocate per-thread local accumulators
    acc128p_t** locals = (acc128p_t**)malloc(threads * sizeof(acc128p_t*));
    rc_worker_args_t* args = (rc_worker_args_t*)malloc(threads * sizeof(rc_worker_args_t));
    pthread_t* tids = (pthread_t*)malloc(threads * sizeof(pthread_t));
    if (!locals || !args || !tids) { abort(); }

    for (int t = 0; t < threads; ++t) {
        // 64B-align to reduce 0 sharing during reduction
        void* p = NULL;
#if defined(_ISOC11_SOURCE)
        p = aligned_alloc(64, ((acc_len * sizeof(acc128p_t) + 63) / 64) * 64);
        if (!p) { abort(); }
#else
        if (posix_memalign(&p, 64, acc_len * sizeof(acc128p_t)) != 0) { abort(); }
#endif
        locals[t] = (acc128p_t*)p;

        size_t chunk = (nb + threads - 1) / threads;
        size_t j0 = (size_t)t * chunk;
        size_t j1 = j0 + chunk; if (j1 > nb) j1 = nb;

        args[t].a = a_copy; args[t].na = na;
        args[t].b = m->limb; args[t].nb = nb;
        args[t].j_begin = j0; args[t].j_end = j1;
        args[t].local_acc = locals[t];
        args[t].ncols = ncols;

        // spawn
        int rc = pthread_create(&tids[t], NULL, worker_rowblock, &args[t]);
        if (rc != 0) { abort(); }
    }

    // Join all workers
    for (int t = 0; t < threads; ++t) {
        pthread_join(tids[t], NULL);
    }

    // Reduce locals -> global accumulator
    acc128p_t* acc = NULL;
#if defined(_ISOC11_SOURCE)
    acc = (acc128p_t*)aligned_alloc(64, ((acc_len * sizeof(acc128p_t) + 63) / 64) * 64);
    if (!acc) { abort(); }
#else
    void* accp = NULL;
    if (posix_memalign(&accp, 64, acc_len * sizeof(acc128p_t)) != 0) { abort(); }
    acc = (acc128p_t*)accp;
#endif
    memset(acc, 0, acc_len * sizeof(acc128p_t));

    for (int t = 0; t < threads; ++t) {
        acc_vec_add(acc, locals[t], acc_len);
    }

    // Expand to base-2^64 lanes and do one global carry pass
    uint64_t* out = (uint64_t*)calloc(out_len, sizeof(uint64_t));
    if (!out) { abort(); }
    expand_and_carry(acc, acc_len, out);

    // Normalize: find actual limb length (trim leading zeros)
    size_t rn = out_len;
    while (rn > 0 && out[rn - 1] == 0) { --rn; }
    if (rn == 0) { cfx_big_set_val(b, 0); }
    else {
        cfx_big_reserve(b, rn);
        // (Assumes cfx_big_reserve zeros new space; if not, it's fine—we overwrite.)
        memcpy(b->limb, out, rn * sizeof(uint64_t));
        b->n = rn;
    }

    // Cleanup
    free(out);
    free(acc);
    for (int t = 0; t < threads; ++t) free(locals[t]);
    free(locals);
    free(args);
    free(tids);
    free(a_copy);
}

void cfx_big_mul_auto(cfx_big_t* b, const cfx_big_t* m) {
    const size_t na = b->n;
    const size_t nb = m->n;

    // trivial cases
    if (!na || !nb) { cfx_big_set_val(b, 0); return; }
    if (nb == 1 && m->limb[0] == 1) return;           // b *= 1
    if (na == 1 && b->limb[0] == 1) {                  // 1 * m
        cfx_big_reserve(b, nb);
        memcpy(b->limb, m->limb, nb * sizeof(uint64_t));
        b->n = nb;
        return;
    }

    const size_t mn = (na < nb) ? na : nb;

    if (mn < 256) {
        // small/skinny -> classic single-thread
        cfx_big_mul(b, m);
        return;
    }

    // Large enough: threaded row/col/carry.
    // Make sure the *row* operand is the larger one to expose threads.
    if (nb >= na) {
        // m already larger (or equal): rows = m
        cfx_big_mul_rows_pthreads(b, m, -1);  // -1 = auto threads in your impl
    } else {
        // m is smaller: compute (m *= b) into a temp, then move into b
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, m);                // tmp = m
        cfx_big_mul_rows_pthreads(&tmp, b, -1);
        cfx_big_swap(b, &tmp);
        cfx_big_free(&tmp);
    }
}

#ifndef CFX_DEC_CHUNK_DIG
#define CFX_DEC_CHUNK_DIG 18u /* 10^18 fits in uint64_t */
#endif

#ifndef CFX_HEX_CHUNK_DIG
#define CFX_HEX_CHUNK_DIG 15u /* 16^15 = 2^60 fits in uint64_t */
#endif


/* Precomputed 10^k for k=0..18 */
static const uint64_t POW10U64[CFX_DEC_CHUNK_DIG + 1] = {
    1ULL,
    10ULL,
    100ULL,
    1000ULL,
    10000ULL,
    100000ULL,
    1000000ULL,
    10000000ULL,
    100000000ULL,
    1000000000ULL,
    10000000000ULL,
    100000000000ULL,
    1000000000000ULL,
    10000000000000ULL,
    100000000000000ULL,
    1000000000000000ULL,
    10000000000000000ULL,
    100000000000000000ULL,
    1000000000000000000ULL,
};

static void flush_chunk(cfx_big_t* out,
                        unsigned base,
                        uint64_t* chunk_val,
                        unsigned* chunk_len)
{
    if (*chunk_len == 0) return;

    if (base == 10) {
        extern const uint64_t POW10U64[]; /* from earlier */
        cfx_big_mul_sm(out, POW10U64[*chunk_len]);
        cfx_big_add_sm(out, *chunk_val);
    } else if (base == 16) {
        unsigned bits = 4u * (*chunk_len);
        cfx_big_shl_bits_eq(out, bits);
        cfx_big_add_sm(out, *chunk_val);
    } else { /* base 2 */
        unsigned bits = *chunk_len;
        cfx_big_shl_bits_eq(out, bits);
        cfx_big_add_sm(out, *chunk_val);
    }

    *chunk_val = 0;
    *chunk_len = 0;
}

/* Return 0 on success, nonzero on parse error. */
int cfx_big_from_file(cfx_big_t* out, FILE* fp, int base) {
    cfx_big_set_val(out, 0);

    enum { BASE_DEC=10, BASE_HEX=16, BASE_BIN=2 };
    int detected_base; 
    if ((base != BASE_DEC) && (base != BASE_HEX) && (base != BASE_BIN))  {
        detected_base = BASE_DEC; /* default */
    } else {
        detected_base = base;
    }
    
    int saw_digit = 0;
    int negative = 0;
    int in_prefix = 1; /* we’re skipping leading ws, sign, 0x */

    /* chunk accumulators */
    uint64_t chunk_val = 0;
    unsigned chunk_len = 0; /* digits in current chunk */

    unsigned dec_chunk_max = CFX_DEC_CHUNK_DIG;
    unsigned hex_chunk_max = CFX_HEX_CHUNK_DIG;

    unsigned char buf[64 * 1024];
    size_t nread;
    uint64_t cnt = 0;
    while ((nread = fread(buf, 1, sizeof(buf), fp)) > 0) {
        cnt += nread;
        printf("read %llu characters..... \n", cnt);
        for (size_t i = 0; i < nread; ++i) {
            unsigned char c = buf[i];

            /* Allow underscores, quotes, spaces, newlines, tabs as visual separators */
            if ( (c == '\n') ||(c == '_') || (c == '"') || isspace(c) || (c == '\t') || (c == '\r') )  continue;

            if (in_prefix) {
                if (isspace(c)) continue;
                if (c == '+') { in_prefix = 0; continue; }
                if (c == '-') { negative = 1; in_prefix = 0; continue; }

                /* Base detection: 0x / 0X (hex), 0b / 0B (bin) */
                if (c == '0') {
                    /* Peek ahead safely: if at buffer end, we’ll detect on next loop */
                    unsigned char next = (i + 1 < nread) ? buf[i + 1] : 0;
                    if (next == 'x' || next == 'X') {
                        if (detected_base != BASE_HEX) {
                            CFX_PRINT_ERR("hex '0x' prefix in file,"
                                " but different base (%d) specified!", detected_base);
                            
                        }
                        detected_base = BASE_HEX;
                        i++;
                        in_prefix = 0;
                        continue;
                    }
                    if (next == 'b' || next == 'B') {
                        if (detected_base != BASE_BIN) {
                            CFX_PRINT_ERR("binary '0b' prefix in file,"
                                " but different base (%d) specified!", detected_base);
                                return -1;
                        }
                        detected_base = BASE_BIN;
                        i++;
                        in_prefix = 0;
                        continue;
                    }
                    /* Otherwise, treat as decimal 0 and fall through as a digit */
                }
                /* We’ve seen a non-space, non-sign */
                in_prefix = 0;
            }

            /* Digit handling by base */
            if (detected_base == BASE_DEC) {
                if (c >= '0' && c <= '9') {
                    saw_digit = 1;
                    if (chunk_len == dec_chunk_max) flush_chunk(out, detected_base, &chunk_val, &chunk_len);
                    chunk_val = chunk_val * 10u + (uint64_t)(c - '0');
                    chunk_len++;
                    continue;
                }
            } else if (detected_base == BASE_HEX) {
                int v = hex_table[c];

                if (v != -1) {
                    saw_digit = 1;
                    if (chunk_len == hex_chunk_max) flush_chunk(out, detected_base, &chunk_val, &chunk_len);
                    chunk_val = (chunk_val << 4) | (uint64_t)v;
                    chunk_len++;
                    continue;
                }
            } else { /* BASE_BIN */
                if (c == '0' || c == '1') {
                    saw_digit = 1;
                    if (chunk_len == 60u) flush_chunk(out, detected_base, &chunk_val, &chunk_len); /* keep under 64 bits */
                    chunk_val = (chunk_val << 1) | (uint64_t)(c - '0');
                    chunk_len++;
                    continue;
                }
            }

            /* Any other non-space char terminates the number (or is an error). */
            if (isspace(c)) {
                /* End of the number */
                goto done_reading;
            } else {
                /* Invalid character in the number */
                // errno = EINVAL;
                CFX_PRINT_ERR("Invalid character found: '%c' (0x%x)!", c, c);
                return -1;
            }
        }
    }

done_reading:
    /* Flush any pending chunk */
    flush_chunk(out, base, &chunk_val, &chunk_len);

    if (!saw_digit) {
        CFX_PRINT_ERR("didn't find any digit!");
        return -1;
    }

    /* TODO: track sign */
    if (negative) {
        /* reject negatives for now */
        // errno = ERANGE; /* or implement signed handling */
        CFX_PRINT_ERR("negative number found");
        return -1;
    }

    return 0;
}

/* ==================== Montgomery ==================== */

/* if x >= n then x -= n */
static inline void cfx_big_cond_sub_n_(cfx_big_t* x, const cfx_big_t* n) {
    if (cfx_big_cmp(x, n) >= 0) cfx_big_sub(x, n);
}

/* n0^{-1} mod 2^64 via Newton; n0 must be odd */
static inline uint64_t inv64_mod2_64_(uint64_t n0) {
    uint64_t x = 1;
    for (int i = 0; i < 6; ++i) x *= (2 - x * n0);
    return x;
}
static inline uint64_t mont_n0inv_(uint64_t n0) {
    return (uint64_t)(0 - inv64_mod2_64_(n0)); /* -n^{-1} mod 2^64 */
}

/* rr = R^2 mod n, with R=2^(64*k). Build rr = 2^(128*k) mod n by repeated doubling. */
static int compute_rr_(cfx_big_t* rr, const cfx_big_t* n, size_t k) {
    cfx_big_set_val(rr, 1);
    for (size_t bit = 0; bit < 128ull * k; ++bit) {
        /* rr = (rr + rr) mod n (rr < n ⇒ rr+rr < 2n ⇒ one cond. subtract is enough) */
        cfx_big_reserve(rr, rr->n + 1);
        __uint128_t carry = 0;
        size_t i = 0;
        for (; i < rr->n; ++i) {
            __uint128_t s = (__uint128_t)rr->limb[i] * 2 + carry;
            rr->limb[i] = (uint64_t)s;
            carry = s >> 64;
        }
        if (carry) rr->limb[rr->n++] = (uint64_t)carry;
        cfx_big_cond_sub_n_(rr, n);
    }
    return 1;
}

int cfx_big_mont_ctx_init(cfx_big_mont_ctx_t* ctx, const cfx_big_t* n_in) {
    if (!ctx || !n_in || n_in->n == 0) return 0;
    if ((n_in->limb[0] & 1ull) == 0) return 0; /* modulus must be odd */

    cfx_big_init(&ctx->n);
    cfx_big_init(&ctx->rr);

    cfx_big_assign(&ctx->n, n_in);
    cfx_big_trim(&ctx->n);

    ctx->k = ctx->n.n;
    if (ctx->k == 0) { cfx_big_mont_ctx_free(ctx); return 0; }

    ctx->n0inv = mont_n0inv_(ctx->n.limb[0]);

    if (!compute_rr_(&ctx->rr, &ctx->n, ctx->k)) {
        cfx_big_mont_ctx_free(ctx);
        return 0;
    }
    return 1;
}

void cfx_big_mont_ctx_free(cfx_big_mont_ctx_t* ctx) {
    if (!ctx) return;
    cfx_big_free(&ctx->n);
    cfx_big_free(&ctx->rr);
    memset(ctx, 0, sizeof(*ctx));
}

/* out = a*b*R^{-1} mod n ; a,b in Montgomery domain. Alias-safe. */
int cfx_big_mont_mul(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, const cfx_big_mont_ctx_t* ctx) {
    if (!out || !a || !b || !ctx) return 0;

    const cfx_big_t* n = &ctx->n;
    const size_t k = ctx->k;
    const uint64_t n0inv = ctx->n0inv;

    cfx_big_t T;
    cfx_big_init(&T);
    cfx_big_reserve(&T, k + 1);
    /* >>>> reserve zeros the newly added limbs, so this is unnecessary >>>> */
    memset(T.limb, 0, (k + 1) * sizeof(uint64_t));
    T.n = k + 1;

    const uint64_t* an = a->limb;
    const uint64_t* bn = b->limb;
    const uint64_t* nn = n->limb;

    for (size_t i = 0; i < k; ++i) {
        const uint64_t bi = (i < b->n) ? bn[i] : 0;

        /* T += a*bi */
        __uint128_t carry = 0;
        for (size_t j = 0; j < k; ++j) {
            const uint64_t aj = (j < a->n) ? an[j] : 0;
            __uint128_t sum = (__uint128_t)T.limb[j] + (__uint128_t)aj * bi + carry;
            T.limb[j] = (uint64_t)sum;
            carry = sum >> 64;
        }
        __uint128_t top = (__uint128_t)T.limb[k] + carry;
        T.limb[k] = (uint64_t)top;

        /* m = T[0]*n0inv mod 2^64 */
        const uint64_t m = T.limb[0] * n0inv;

        /* T += m*n */
        __uint128_t carry2 = 0;
        for (size_t j = 0; j < k; ++j) {
            __uint128_t sum = (__uint128_t)T.limb[j] + (__uint128_t)m * nn[j] + carry2;
            T.limb[j] = (uint64_t)sum;
            carry2 = sum >> 64;
        }
        __uint128_t top2 = (__uint128_t)T.limb[k] + carry2;
        T.limb[k] = (uint64_t)top2;

        /* T = (T + m*n) / R : drop limb 0 */
        memmove(&T.limb[0], &T.limb[1], k * sizeof(uint64_t));
        T.limb[k] = 0; /* scratch for next iter */
    }

    T.n = k;
    cfx_big_trim(&T);
    if (cfx_big_cmp(&T, n) >= 0) cfx_big_sub(&T, n);

    cfx_big_assign(out, &T);
    cfx_big_free(&T);
    return 1;
}

/* Convert to Montgomery: aR mod n = MontMul(a, R^2) */
int cfx_big_mont_to(cfx_big_t* out, const cfx_big_t* a, const cfx_big_mont_ctx_t* ctx) {
    if (!out || !a || !ctx) return 0;
    /* If your inputs may be >> n, call your real mod here. For now, quick reduce by repeated sub when a < few*n. */
    cfx_big_t ared;
    cfx_big_init(&ared);
    cfx_big_assign(&ared, a);
    while (cfx_big_cmp(&ared, &ctx->n) >= 0) cfx_big_sub(&ared, &ctx->n);

    int ok = cfx_big_mont_mul(out, &ared, &ctx->rr, ctx);
    cfx_big_free(&ared);
    return ok;
}

/* Convert from Montgomery: aR * R^{-1} = a mod n = MontMul(aR, 1) */
int cfx_big_mont_from(cfx_big_t* out, const cfx_big_t* aR, const cfx_big_mont_ctx_t* ctx) {
    if (!out || !aR || !ctx) return 0;
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_set_val(&one, 1);
    int ok = cfx_big_mont_mul(out, aR, &one, ctx);
    cfx_big_free(&one);
    return ok;
}

/* ----------------- One-liners that hide the ctx ----------------- */
int cfx_big_mul_mod(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, const cfx_big_t* n) {
    cfx_big_mont_ctx_t C;
    if (!cfx_big_mont_ctx_init(&C, n)) return 0;
    cfx_big_t aR, bR, r;
    cfx_big_init(&aR);
    cfx_big_init(&bR);
    cfx_big_init(&r);

    int ok = cfx_big_mont_to(&aR, a, &C)
           && cfx_big_mont_to(&bR, b, &C)
           && cfx_big_mont_mul(&r, &aR, &bR, &C)
           && cfx_big_mont_from(out, &r, &C);

    cfx_big_free(&aR);
    cfx_big_free(&bR);
    cfx_big_free(&r);
    cfx_big_mont_ctx_free(&C);
    return ok;
}

int cfx_big_sqr_mod(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* n) {
    cfx_big_mont_ctx_t C; if (!cfx_big_mont_ctx_init(&C, n)) return 0;
    cfx_big_t aR, r;
    cfx_big_init(&aR); cfx_big_init(&r);
    int ok = cfx_big_mont_to(&aR, a, &C)
           && cfx_big_mont_sqr(&r, &aR, &C)
           && cfx_big_mont_from(out, &r, &C);
    cfx_big_free(&aR); cfx_big_free(&r);
    cfx_big_mont_ctx_free(&C);
    return ok;
}

int cfx_big_modexp(cfx_big_t* out, const cfx_big_t* base, const cfx_big_t* exp, const cfx_big_t* n) {
    cfx_big_mont_ctx_t C; if (!cfx_big_mont_ctx_init(&C, n)) return 0;

    cfx_big_t A, X;
    cfx_big_init(&A); cfx_big_init(&X);

    int ok = cfx_big_mont_to(&A, base, &C);
    if (ok) {
        cfx_big_t one; cfx_big_init(&one); cfx_big_set_val(&one, 1);
        ok = cfx_big_mont_to(&X, &one, &C);
        cfx_big_free(&one);
    }
    if (ok) {
        for (size_t i = 0; i < exp->n; ++i) {
            uint64_t w = exp->limb[i];
            for (int b = 0; b < 64; ++b) {
                if (w & 1u) { ok = cfx_big_mont_mul(&X, &X, &A, &C); if (!ok) break; }
                w >>= 1;
                ok = cfx_big_mont_mul(&A, &A, &A, &C); if (!ok) break; /* square */
            }
            if (!ok) break;
        }
        if (ok) ok = cfx_big_mont_from(out, &X, &C);
    }

    cfx_big_free(&A); cfx_big_free(&X);
    cfx_big_mont_ctx_free(&C);
    return ok;
}
