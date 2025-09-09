#include "cfx/algo.h"
#include "cfx/types.h"
#include "cfx/macros.h"
#include "cfx/primes.h"

#include <stdlib.h>

/* bit index i -> x = 2i + 3 */
#define BIT_GET(A,i)  (((A)[(i)>>3] >> ((i)&7)) & 1u)
#define BIT_SET(A,i)  ((A)[(i)>>3] |= (uint8_t)(1u << ((i)&7)))

static inline uint64_t isqrt_u64(uint64_t n) {
    // Integer sqrt via double for speed, then correct by adjustment.
    // No UB and avoids i*i overflow checks in loops.
    double d = (double)n;
    uint64_t x = (uint64_t)(d > 0 ? __builtin_floor(__builtin_sqrt(d)) : 0);
    while ((x+1) > 0 && (x+1) <= n / (x+1)) ++x;
    while (x > 0 && x > n / x) --x;
    return x;
}

/* Sieve of Eratosthenes (odd-only, bitset), primes up to and including n */
cfx_vec_t cfx_sieve_primes(uint64_t n) {
    cfx_vec_t primes = {0};
    if (n < 2) return primes;

    // Map odds only: index k (0..m-1) corresponds to value p = 2*k + 3
    uint64_t m = (n >= 3) ? ((n - 3) >> 1) + 1 : 0;      // #odd candidates in [3..n]
    size_t bytes = (size_t)((m + 7) >> 3);
    uint8_t *mark = (uint8_t*)calloc(bytes ? bytes : 1, 1);
    if (!mark) return primes;

    uint64_t r = isqrt_u64(n);

    // Sieve: i runs over odd primes; index ki maps i = 2*ki+3
    for (uint64_t i = 3; i <= r; i += 2) {
        uint64_t ki = (i - 3) >> 1;
        if (!BIT_GET(mark, ki)) {
            // Start at j = i*i (odd), step by 2*i (skip evens)
            uint128_t ii = ( uint128_t)i * ( uint128_t)i; // avoid overflow on i*i
            uint64_t j = (uint64_t)ii;
            uint64_t step = i << 1; // 2*i

            // Convert j to bit index kj: j = 2*kj + 3  => kj = (j-3)/2
            uint64_t kj = (j - 3) >> 1;
            for (; kj < m; kj += step >> 1) { // step/2 in bit-index space
                BIT_SET(mark, kj);
                j += step; // (not strictly needed since we advance kj directly)
            }
        }
    }

    // Count primes
    size_t count = 0;
    count += 1; // prime 2
    for (uint64_t k = 0; k < m; ++k) if (!BIT_GET(mark, k)) ++count;

    // Emit primes
    primes.data = (uint64_t*)malloc(count * sizeof(uint64_t));
    if (!primes.data) { free(mark); return (cfx_vec_t){0}; }
    size_t pos = 0;
    primes.data[pos++] = 2;
    for (uint64_t k = 0; k < m; ++k) {
        if (!BIT_GET(mark, k)) primes.data[pos++] = (k << 1) + 3;
    }
    primes.size = pos;

    free(mark);
    return primes;
}


/* uses legendre's formula to give th power of p that divides n! */
uint64_t cfx_legendre(uint64_t n, uint64_t p) {
    uint64_t e = 0;
    while (n) {
        n /= p;
        e += n;
    }
    return e;
}

static inline uint64_t _cfx_mulmod_u64(uint64_t a, uint64_t b, uint64_t m) {
    uint128_t p = (uint128_t)a * b;
    return (uint64_t) (p % m);
}

/* binary exponentiation mod m */
static inline uint64_t _cfx_powmod_u64(uint64_t a, uint64_t e, uint64_t m) {
    if (m == 0) return 0;
    if (m == 1) return 0;
    a %= m;
    uint64_t r = 1;
    while (e) {
        if (e & 1) r = _cfx_mulmod_u64(r, a, m);
        a = _cfx_mulmod_u64(a, a, m);
        e >>= 1;
    }
    return r;
}

uint64_t cfx_mulmod_u64(uint64_t a, uint64_t b, uint64_t m) {
    return _cfx_mulmod_u64(a, b, m); 
}

/* modular exponentiation (a^e) mod m */
uint64_t cfx_powmod_u64(uint64_t a, uint64_t e, uint64_t m) {
    return _cfx_powmod_u64(a, e, m);
}

/** Uses Miller–Rabin 'probable primality' test 
 * - take n-1 (which is even), express it as d*2^s.
 * - find a which is coprime to n
 * - then n is "probably" prime if either:
 *      - a^d == 1 mod n, or 
 *      - a^(2r*d) == -1 mod n for some 0 <= r < s;
*/
int cfx_is_prime_u64(uint64_t n) {
   if (n < 2) return 0;

    // small primes precheck
    static const uint32_t small[] = {2,3,5,7,11,13,17,19,23,29,31,37,0};
    for (int i = 0; small[i]; ++i) {
        uint32_t p = small[i];
        if (n == p) return 1;
        if (n % p == 0) return 0;
    }

    // write n-1 = d * 2^s with d odd
    uint64_t d = n - 1, s = 0;
    while ((d & 1) == 0) { d >>= 1; ++s; }

    /* The "Sinclair 7" bases, determinstic for 64 bit */
    static const uint64_t bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};

    for (int i = 0; bases[i]; ++i) {
        uint64_t a = bases[i];
        if (a % n == 0) continue;             // skip if a == n
        uint64_t x = _cfx_powmod_u64(a, d, n);
        if (x == 1 || x == n - 1) continue;   // this base passes

        int witness = 1;
        for (uint64_t r = 1; r < s; ++r) {
            x = _cfx_mulmod_u64(x, x, n);
            if (x == n - 1) { witness = 0; break; }
        }
        if (witness) return 0; // composite
    }
    return 1; // prime
}


uint64_t cfx_gcd_u64(uint64_t a, uint64_t b) {
    while (b) { uint64_t t = a % b; a = b; b = t; }
    return a;
}

uint64_t cfx_rho_brent(uint64_t n) {
    if ((n & 1) == 0) return 2;
    uint64_t y = (uint64_t)rand() % (n-1) + 1;
    uint64_t c = (uint64_t)rand() % (n-1) + 1;
    uint64_t m = 128;

    uint64_t g = 1, r = 1, q = 1, x, ys;
    while (g == 1) {
        x = y;
        for (uint64_t i=0; i<r; ++i)
            y = (_cfx_mulmod_u64(y,y,n) + c) % n;

        uint64_t k = 0;
        while (k < r && g == 1) {
            ys = y;
            uint64_t lim = (m < (r-k)) ? m : (r-k);
            for (uint64_t i=0; i<lim; ++i) {
                y = (_cfx_mulmod_u64(y,y,n) + c) % n;
                uint64_t diff = x>y ? x-y : y-x;
                q = _cfx_mulmod_u64(q, diff, n);
            }
            g = cfx_gcd_u64(q, n);
            k += lim;
        }
        r <<= 1;
    }

    if (g == n) {
        do {
            ys = (_cfx_mulmod_u64(ys,ys,n) + c) % n;
            uint64_t diff = x>ys ? x-ys : ys-x;
            g = cfx_gcd_u64(diff, n);
        } while (g == 1);
    }
    return g;
}

/* comparator fn for qsort */
static int cmp_u64(const void* a, const void* b) {
    uint64_t x = *(const uint64_t*)a, y = *(const uint64_t*)b;
    return (x < y) ? -1 : (x > y);
}


int cfx_factor_u64(cfx_vec_t* primes, cfx_vec_t* exps, uint64_t n) {
   
    if (n == 0) return -1;
    if (n == 1) return 0;
    
    cfx_vec_t v; /* local vec for all primes */
    cfx_vec_init(&v);

    cfx_vec_init(primes);
    cfx_vec_init(exps);
    
    /* 1) Strip tiny (<255) primes fast */
    for (size_t i = 0; i < 54; ++i) {
        uint64_t p = cfx_primes[i];
        while ((n % p) == 0) {
            cfx_vec_push(&v, p);
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
                    cfx_vec_push(&v, p); m /= p;
                }
            }
            if (m > 1) cfx_vec_push(&v, m);
            continue;
        }

        if (cfx_is_prime_u64(m)) {
            cfx_vec_push(&v, m);
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
                cfx_vec_push(&v, m);
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
    qsort(v.data, v.size, sizeof(uint64_t), cmp_u64);

    for (size_t i = 0; i < v.size; ) {
        uint64_t p = v.data[i];
        uint32_t e = 1;
        size_t j = i + 1;
        while (j < v.size && v.data[j] == p) { ++e; ++j; }
        cfx_vec_push(primes, p);
        cfx_vec_push(exps, e);
        i = j;
    }
    cfx_vec_free(&v);
    return 0;
}
