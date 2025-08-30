#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

 /* a single prime factor: */
typedef struct {
    uint32_t p;  // prime
    uint32_t e;  // exponent
} PF;

/**
 * factorization: a struct to hold a prime factorization of a number
 */
typedef struct {
    PF* data;
    size_t len;
    size_t cap;
} factorization;

static void fac_print(factorization* f) {
    printf("f: %p, cap: %zu, len: %zu\n", f, f->cap, f->len);
    for (size_t k = 0; k < f->len-1; ++k) {
        PF* pf = &f->data[k];
        printf("%u^%u + ", pf->p, pf->e);
    }
    PF* pf = &f->data[f->len-1];
    printf("%u^%u\n", pf->p, pf->e);
}

static void fac_init(factorization* f) {
    f->data = NULL;
    f->len = 0;
    f->cap = 0;
}

static void fac_free(factorization* f) {
    free(f->data);
    fac_init(f);
}

static void fac_reserve(factorization* f, size_t req_cap) {
    if (req_cap <= f->cap) {
        return;
    }
    size_t new_cap = f->cap ? 2 * f->cap : 16;
    if (new_cap < req_cap) {
        new_cap = req_cap;
    }
    // printf("fac_reserve %p requested cap: %zu, new cap: %zu\n", f, req_cap, new_cap);
    f->data = (PF*)realloc(f->data, new_cap * sizeof(PF));
    f->cap = new_cap;
}

static void fac_push(factorization* f, uint32_t p, uint32_t e) {
    if (e == 0) return;
    fac_reserve(f, f->len + 1);
    ++f->len;
    PF* pf = &f->data[f->len - 1];
    pf->p = p;
    pf->e = e;
    
}

/* Deep copy: dst becomes a copy of src, using fac_push for each element. */
static void fac_copy(factorization *dst, const factorization *src) {
    if (dst == src) return;  
    fac_free(dst);
    fac_init(dst);
    fac_reserve(dst, src->len);

    for (size_t i = 0; i < src->len; ++i) {
        fac_push(dst, src->data[i].p, src->data[i].e);
    }
}

// dst += src
static void fac_add(factorization* dst, factorization* src) {
    factorization out;
    fac_init(&out);
    fac_reserve(&out, src->len + dst->len);

    // Assumes data[] is sorted by p
    size_t i = 0;
    size_t j = 0;

    while (i < dst->len || j < src->len) {
        PF* pf1 = &dst->data[i];
        PF* pf2 = &src->data[j];

        if (j == src->len || (i < dst->len && pf1->p < pf2->p)) {
            // p is in dst, but not in src
            fac_push(&out, pf1->p, pf1->e);
            ++i;
        } else if (i == dst->len || (j < src->len && pf1->p > pf2->p)) {
            // p is in src, but not in dst
            fac_push(&out, pf2->p, pf2->e);
            ++j;
        } else {
            // p is in src and dst
            uint32_t p = pf1->p;
            uint32_t e = pf1->e + pf2->e;
            if (e) fac_push(&out, p, e);
            ++i; 
            ++j;
        }
    }
    free(dst->data);
    *dst = out;
}

/* dst -= src */
static void fac_sub(factorization* dst, factorization* src) {
    factorization out;
    
    fac_init(&out);
    fac_reserve(&out, dst->len);

    size_t i = 0;
    size_t j = 0;

    while (i < dst->len || j < src->len) {
        PF* pf1 = &dst->data[i];
        PF* pf2 = &src->data[j];

        if (j == src->len || (i < dst->len && pf1->p < pf2->p)) {
            // p is in dst but not src - no change
            fac_push(&out, pf1->p, pf1->e);
            ++i;
        } else if (i == dst->len || (j < src->len && pf1->p > pf2->p)) {
            // p is in src but not dst - dst does not divide src!!!!
            assert(0 && "fac_sub underflow");
            j++;
        } else {
            // p is in both src and dst:
            uint32_t p = pf1->p;
            uint32_t e;
            if (pf1->e > pf2->e) { // dst divides src
                e = pf1->e - pf2->e;
                fac_push(&out, p, e);
            } else if (pf2->e > pf1->e) { // dst does not divide src!
                assert(0 && "fac_sub underflow");
            } 
            // else, e == 0, p is cancelled out!ç
            ++i;
            ++j;
        }
    }
    free(dst->data);
    *dst = out;
}

/*
---------------------------------------------------------------------------
*/
typedef struct {
    uint32_t *v;
    size_t n;
} u32_vec;

static void u32vec_free(u32_vec* v) {
    free(v->v);
    v->n = 0;
}

static void print_vec(u32_vec* v) {
    printf("v: %p, size: %zu\n", v, v->n);
    for (size_t i = 0; i < v->n; ++i) {
        printf("%u ", v->v[i]);
    }
    printf("\n");
}

/* find all primes up to n using sieve of eratosthenes */
static u32_vec sieve_primes(uint32_t n) {
    u32_vec primes = {0};
    if (n < 2) return primes;

    /* mark all composites, multiples of i until sqrt(n) */
    /* mark[x] == 0: x prime, mark[x] == 1: x composite */
    uint8_t* mark = (uint8_t*)calloc(n + 1, 1);

    for (uint32_t i = 2; (i * i) <= n; ++i) {
        if (!mark[i]) {
            for (uint32_t j = i*i; j <= n; j += i) {
                // printf("marking %u composite\n", j);
                mark[j] = 1;
            }
        }
    }
    
    size_t p_cnt = 0;
    for (uint32_t i = 2; i < n + 1; ++i) {
        if (!mark[i]) ++p_cnt;
    }

    primes.v = (uint32_t*)malloc(p_cnt * sizeof(uint32_t));
    primes.n = p_cnt;

    size_t k = 0;
    for (uint32_t i = 2; i <= n; ++i) {
        if (!mark[i]) {
            primes.v[k] = i;
            ++k;
        }
    }

    free(mark);
    return primes;
}

/* uses legendre's formula to give th power of p that divides n! */
static uint32_t legendre(uint32_t n, uint32_t p) {
    uint32_t e = 0;
    while (n) {
        n /= p;
        e += n;
    }
    return e;
}

/* calculate the factorization of n.
we pass in a list of primes to not have to calculate it on every call of this function */
static factorization fac_factorial(uint32_t n, const u32_vec *primes) {
    factorization f;
    fac_init(&f);
    for (size_t i = 0; i < primes->n; ++i) {
        uint32_t p = primes->v[i];
        if (p > n) break;
        uint32_t e = legendre(n, p);
        if (e) fac_push(&f, p, e);
    }
    return f;
}

/* Factorization of C(n,k) = n! / (k! (n-k)!) */
static factorization fac_binom(uint32_t n, uint32_t k){
    if (k > n) { factorization z; fac_init(&z); return z; }
    if (k > n-k) k = n-k;
    u32_vec primes = sieve_primes(n);
    factorization fn = fac_factorial(n, &primes);
    factorization fk = fac_factorial(k, &primes);
    factorization fnk = fac_factorial(n-k, &primes);
    fac_sub(&fn, &fk);
    fac_sub(&fn, &fnk);
    u32vec_free(&primes);
    return fn; /* sorted, coalesced */
}

/*------------------------------------------------------*/
/* minimal big int struct */
typedef struct {
    uint32_t* limb; /* the 'digits' of the number with base BIG_BASE*/
    size_t n;
    size_t cap;
} big;

/* 1e9 */
#define BIG_BASE 1000000000u

static void big_init(big* b) {
    b->limb = NULL;
    b->n = 0;
    b->cap = 0;
}

static void big_free(big* b) {
    free(b->limb);
    big_init(b);
}

static void big_reserve(big* b, size_t need) {
    if (need < b->cap) return;
    size_t new_cap = b->cap ? 2*b->cap : 16;
    if (new_cap < need) {
        new_cap = need;
    }
    b->limb = (uint32_t*)realloc(b->limb, new_cap*sizeof(uint32_t));
    b->cap = new_cap;
}

static void big_set_u32(big* b, uint32_t v) {
    big_reserve(b, 1);
    b->n = v ? 1:0;
    if (v) b->limb[0] = v;
}

static void big_mul_u32(big* b, uint32_t m) {
    if (m == 1) return;
    if (b->n == 0) {
        big_set_u32(b, 0);
        return;
    }
    uint64_t carry = 0;
    for (size_t i = 0; i < b->n; ++i) {
        uint64_t t = b->limb[i] * m + carry;
        b->limb[i] = (uint32_t)(t % BIG_BASE);
        carry = t / BIG_BASE;
    }
    while (carry) {
        // have to add a digit 
        big_reserve(b, b->n+1);
        b->limb[b->n] = (uint32_t)(carry % BIG_BASE);
        b->n++;
        carry /= BIG_BASE;
    }

}

/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
static void big_powmul_prime(big *b, uint32_t p, uint32_t e){
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

    for (uint32_t i = 0; i < q; i++) big_mul_u32(b, pow_t);

    /* multiply remainder by binary exponentiation (still fits u32) */
    uint32_t rempow = 1;
    acc = p;
    uint32_t rr = r;

    while (rr) {
        if (rr & 1u) rempow = (uint32_t)(rempow * acc);
        rr >>= 1u;
        if (rr) acc = acc*acc;
    }
    if (rempow != 1) big_mul_u32(b, rempow);
}

/* Materialize factorization into big */
static void big_from_factorization(big *b, const factorization *f){
    big_set_u32(b, 1);
    for (size_t i=0;i<f->len;i++){
        big_powmul_prime(b, f->data[i].p, f->data[i].e);
    }
}

/* Convert big to decimal string */
static char* big_to_string(const big *b){
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

/* ---------- Demo: compute C(n,k) ---------- */

static void print_binom(uint32_t n, uint32_t k){
    factorization f = fac_binom(n,k);
    big B; big_init(&B);
    big_from_factorization(&B, &f);
    char *s = big_to_string(&B);
    printf("C(%u, %u) = %s\n", n, k, s);
    free(s);
    big_free(&B);
    fac_free(&f);
}

int main(int argc, char* argv[]) {
    uint32_t n, k;
    if (argc == 2) {
        n = (uint32_t)strtol(argv[1], NULL, 10);
        for (int i = 1; i <= n; ++i) {
            for (int j = 0; j <= i; ++j) {
                print_binom(i, j);
            }
            printf(".....\n");
        }
    } else if (argc == 3) {
        n = (uint32_t)strtol(argv[1], NULL, 10);
        k = (uint32_t)strtol(argv[2], NULL, 10);
        print_binom(n, k);
    }
/*
    factorization f1, f2;

    factorization* p1 = &f1;
    factorization* p2 = &f2;
    
    fac_init(p1);
    fac_init(p2);
    
    fac_push(p1, 2, 4);
    fac_push(p1, 3, 2);
    
    fac_push(p2, 2, 5);
    fac_push(p2, 5, 2);

    printf("p1: \n");
    fac_print(p1);
    printf("p2: \n");
    fac_print(p2);

    printf("p1 += p2 -------\n");
    fac_add(p1, p2);
    fac_print(p1);

    printf("p1 -= p2 -------\n");
    fac_sub(p1, p2);
    fac_print(p1);

    factorization f3;
    fac_copy(&f3, &f1);
    printf("f3 copy of p1 -------\n");
    fac_print(&f3);

    printf("primes until %u:\n", n);
    u32_vec primes = sieve_primes(n);
    print_vec(&primes);

    printf("-------------- \n");
    factorization f4 = fac_factorial(n, &primes);
    fac_print(&f4);

    print_binom(n, k);
*/
    

    return 0;
}