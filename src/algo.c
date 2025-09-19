#include "cfx/algo.h"
#include "cfx/types.h"
#include "cfx/macros.h"
#include "cfx/primes.h"

#include <stdlib.h>
#include <stdint.h>
#include <string.h>   // memset
#include <assert.h>

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
    static const uint64_t bases[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022, 0};

    for (int i = 0; i < sizeof(bases)/sizeof(bases[0]); ++i) {
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



/* Add a 128-bit value (add_hi:add_lo) into one accumulator cell */
static inline void csa_add128_into(csa128_t* acc, uint64_t add_lo, uint64_t add_hi) {
    // 128-bit add: (acc_hi:acc_lo) += (add_hi:add_lo)   without leaking carries to neighbors
    uint128_t t  = (uint128_t)acc->lo + add_lo;
    acc->lo        = (uint64_t)t;
    acc->hi       += add_hi + (uint64_t)(t >> 64);  // keep the carry in the local 'hi' word
}

// Zero only the first `nout` entries that the multiplication will touch
void cfx_mul_scratch_zero(cfx_mul_scratch_t* s, size_t nout) {
    assert(nout <= s->cap);
    memset(s->acc,   0, nout * sizeof(*s->acc));
    memset(s->spill, 0, (nout+1) * sizeof(uint64_t));
}

void cfx_mul_scratch_alloc(cfx_mul_scratch_t* s, size_t needed) {
    if (s->cap >= needed) return;  // already big enough

    // free old
    free(s->acc);
    free(s->spill);

    /* or use typeof(s->acc) :D gcc extension*/
    s->acc   = (csa128_t*)calloc(needed, sizeof(*s->acc));
    s->spill = (uint64_t*)calloc(needed+1, sizeof(uint64_t)); // we write spill[k+1] when k==nout-1)
    s->cap   = needed;

    assert(s->acc && s->spill);
}

void cfx_mul_scratch_free(cfx_mul_scratch_t* s) {
    free(s->acc);
    free(s->spill);
    s->acc = NULL;
    s->spill = NULL;
    s->cap = 0;
}

// Assumes scratch->acc[0..nout-1] and scratch->spill[0..nout-1] are zeroed by the caller.
void cfx_mul_csa_portable_fast(const uint64_t* A, size_t na,
                               const uint64_t* B, size_t nb,
                               uint64_t* R,
                               cfx_mul_scratch_t* scratch)
{
    const size_t nout = na + nb;
    csa128_t* acc = scratch->acc;
    uint64_t* spill = scratch->spill;

    // --- Accumulate all partial products (no ripple spills) ---
    for (size_t i = 0; i < na; ++i) {
        uint64_t ai = A[i];
        for (size_t j = 0; j < nb; ++j) {
            uint128_t p  = (uint128_t)ai * (uint128_t)B[j];
            uint64_t add_lo = (uint64_t)p;
            uint64_t add_hi = (uint64_t)(p >> 64);
            size_t k = i + j;

            uint128_t t = (uint128_t)acc[k].lo + add_lo;
            uint64_t c0   = (uint64_t)(t >> 64);
            acc[k].lo     = (uint64_t)t;

            uint128_t u = (uint128_t)acc[k].hi + add_hi + c0;
            acc[k].hi     = (uint64_t)u;

            // Defer the overflow of acc[k].hi into a 1-bit counter for next diagonal.
            spill[k + 1] += (uint64_t)(u >> 64);  // 0 or 1
        }
    }

    // --- Fold deferred spills into acc[].hi with a single linear propagate ---
    // spill[i] can become >1, so we propagate once left-to-right; this is cheap & predictable.
    uint64_t c = 0;
    for (size_t k = 0; k < nout; ++k) {
        uint128_t h = (uint128_t)acc[k].hi + c;
        acc[k].hi = (uint64_t)h;
        c = (uint64_t)(h >> 64);
        // carry from adding c merges with next diagonal’s spill counter
        spill[k + 1] += c;   // safe: spill has nout entries; spill[nout] is ignored later
        c = spill[k + 1];
    }

    // --- Final normalization (single carry-propagating pass across limbs) ---
    uint128_t carry = 0;
    for (size_t k = 0; k < nout; ++k) {
        uint128_t wide = (((uint128_t)acc[k].hi) << 64) | acc[k].lo;
        wide += carry;
        R[k]  = (uint64_t)wide;
        carry = wide >> 64;
    }
    // By construction carry==0 for na+nb limbs.
}


void cfx_mul_csa_portable(const uint64_t* A, size_t na,
                          const uint64_t* B, size_t nb, uint64_t* R) {

    assert(na > 0 && nb > 0);
    const size_t nout = na + nb;              // result limbs
    csa128_t* acc = (csa128_t*)alloca(nout * sizeof(csa128_t));
    for (size_t k = 0; k < nout; ++k) { acc[k].lo = 0; acc[k].hi = 0; }

    // Accumulate all partial products per diagonal, carry-save style
    for (size_t i = 0; i < na; ++i) {
        uint64_t ai = A[i];
        for (size_t j = 0; j < nb; ++j) {
            uint128_t p  = (uint128_t)ai * (uint128_t)B[j];
            uint64_t add_lo = (uint64_t)p;
            uint64_t add_hi = (uint64_t)(p >> 64);
            size_t k = i + j;

            // (acc_hi:acc_lo) += (add_hi:add_lo)
            uint128_t t = (uint128_t)acc[k].lo + add_lo;
            uint64_t new_lo = (uint64_t)t;
            uint64_t c0     = (uint64_t)(t >> 64);

            uint128_t u = (uint128_t)acc[k].hi + add_hi + c0;
            acc[k].lo = new_lo;
            acc[k].hi = (uint64_t)u;

            // *** spill overflow of the diagonal's hi into the NEXT diagonal's hi ***
            uint64_t spill = (uint64_t)(u >> 64);   // 0 or 1
            size_t kk = k + 1;
            while (spill && kk < nout) {
                uint128_t w = (uint128_t)acc[kk].hi + spill;
                acc[kk].hi = (uint64_t)w;
                spill = (uint64_t)(w >> 64);        // might ripple further
                ++kk;
            }
            // If spill is still nonzero here, it would be beyond the m+n limb bound.
            // For 64-bit limb products of sizes na, nb, true products fit in exactly na+nb limbs.
            assert(spill == 0);
        }
    }

    // Final single carry-propagation across limbs
    uint128_t carry = 0;
    for (size_t k = 0; k < nout; ++k) {
        uint128_t wide = (((uint128_t)acc[k].hi) << 64) | acc[k].lo;
        wide += carry;
        R[k]  = (uint64_t)wide;
        carry = wide >> 64;
    }

    // By construction, carry should be 0 here for na+nb limbs.
    assert(carry == 0);

    // CFX_PRINT_DBG("A:\t");
    // PRINT_ARR(A, na);
    // CFX_PRINT_DBG("B:\t");
    // PRINT_ARR(B, nb);
    // CFX_PRINT_DBG("R:\t");
    // PRINT_ARR(R, nout);
}

// Accumulate rows i in [ia, ia+na_rows) of A against full B into acc/spill (size nout).
void cfx_mul_csa_portable_fast_rows(const uint64_t* A, size_t ia, size_t na_rows,
                                    const uint64_t* B, size_t nb,
                                    csa128_t* acc, uint64_t* spill)
{
    for (size_t ii = 0; ii < na_rows; ++ii) {
        size_t i = ia + ii;
        uint64_t ai = A[i];
        for (size_t j = 0; j < nb; ++j) {
            __uint128_t p  = (__uint128_t)ai * (__uint128_t)B[j];
            uint64_t add_lo = (uint64_t)p;
            uint64_t add_hi = (uint64_t)(p >> 64);
            size_t k = i + j;

            __uint128_t t = (__uint128_t)acc[k].lo + add_lo;
            uint64_t c0   = (uint64_t)(t >> 64);
            acc[k].lo     = (uint64_t)t;

            __uint128_t u = (__uint128_t)acc[k].hi + add_hi + c0;
            acc[k].hi     = (uint64_t)u;

            spill[k + 1] += (uint64_t)(u >> 64);
        }
    }
}

void cfx_mul_csa_fold_and_normalize(csa128_t* acc, uint64_t* spill,
                                    size_t nout, uint64_t* R)
{
    // Fold spills left->right (predictable)
    uint64_t pending = spill[0];
    for (size_t k = 0; k < nout; ++k) {
        __uint128_t h = (__uint128_t)acc[k].hi + pending;
        acc[k].hi = (uint64_t)h;
        pending = (uint64_t)(h >> 64);
        pending += spill[k + 1]; // spill sized nout+1
    }
    // Final normalization
    __uint128_t carry = 0;
    for (size_t k = 0; k < nout; ++k) {
        __uint128_t w = (((__uint128_t)acc[k].hi) << 64) | acc[k].lo;
        w += carry;
        R[k]  = (uint64_t)w;
        carry = (uint64_t)(w >> 64);
    }
}
