/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/* ntt.c - Number-Theoretic Transform implementation */

#include "cfx/ntt.h"
#include "cfx/arch.h"
#include <stdlib.h>
#include <string.h>

/*
 * Modular Arithmetic
 *
 * For NTT we need fast modular multiplication of 64-bit values mod a 64-bit prime.
 * The product can be up to 128 bits, so we need either:
 *   - Native __uint128_t (GCC/Clang on 64-bit)
 *   - Portable 64x64->128 via 32-bit halves
 */

uint64_t cfx_ntt_mod_add(uint64_t a, uint64_t b, uint64_t p) {
    uint64_t sum = a + b;
    /* if overflow or sum >= p, subtract p */
    if (sum < a || sum >= p) {
        sum -= p;
    }
    return sum;
}

uint64_t cfx_ntt_mod_sub(uint64_t a, uint64_t b, uint64_t p) {
    if (a >= b) {
        return a - b;
    } else {
        return p - (b - a);
    }
}

/*
 * 64x64 -> 128-bit multiplication, then mod p
 *
 * We compute (a * b) mod p where a, b, p are all 64-bit.
 */
#if CFX_HAS_UINT128

uint64_t cfx_ntt_mod_mul(uint64_t a, uint64_t b, uint64_t p) {
    __uint128_t prod = (__uint128_t)a * (__uint128_t)b;
    return (uint64_t)(prod % p);
}

#else

/* portable 64x64 -> 128-bit product, then division */
static void mul128(uint64_t a, uint64_t b, uint64_t* hi, uint64_t* lo) {
    uint64_t a_lo = (uint32_t)a;
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = (uint32_t)b;
    uint64_t b_hi = b >> 32;

    uint64_t p0 = a_lo * b_lo;
    uint64_t p1 = a_lo * b_hi;
    uint64_t p2 = a_hi * b_lo;
    uint64_t p3 = a_hi * b_hi;

    uint64_t mid = (p0 >> 32) + (uint32_t)p1 + (uint32_t)p2;
    *lo = (p0 & 0xFFFFFFFF) | (mid << 32);
    *hi = p3 + (p1 >> 32) + (p2 >> 32) + (mid >> 32);
}

/* 128-bit mod 64-bit using binary long division
 * For NTT primes near 2^64, we need to handle overflow when shifting rem.
 * We track a conceptual 65-bit value as (carry, rem).
 */
static uint64_t mod128(uint64_t hi, uint64_t lo, uint64_t p) {
    if (hi == 0) {
        return lo % p;
    }

    /* precompute 2^64 mod p for overflow handling
     * 2^64 mod p = 2^64 - p (since p < 2^64 and 2^64 = 1*p + (2^64 - p))
     * in unsigned arithmetic: 0 - p = 2^64 - p
     */
    uint64_t pow2_64_mod_p = (uint64_t)0 - p;

    uint64_t rem = 0;

    for (int i = 63; i >= 0; --i) {
        int carry = (rem >> 63) & 1;
        rem = (rem << 1) | ((hi >> i) & 1);
        if (carry) {
            /* conceptual value was 2^64 + rem, now add 2^64 mod p */
            rem += pow2_64_mod_p;
            if (rem >= p) rem -= p;
        }
        if (rem >= p) rem -= p;
    }
    for (int i = 63; i >= 0; --i) {
        int carry = (rem >> 63) & 1;
        rem = (rem << 1) | ((lo >> i) & 1);
        if (carry) {
            rem += pow2_64_mod_p;
            if (rem >= p) rem -= p;
        }
        if (rem >= p) rem -= p;
    }
    return rem;
}

uint64_t cfx_ntt_mod_mul(uint64_t a, uint64_t b, uint64_t p) {
    uint64_t hi, lo;
    mul128(a, b, &hi, &lo);
    return mod128(hi, lo, p);
}

#endif /* CFX_HAS_UINT128 */

uint64_t cfx_ntt_mod_pow(uint64_t base, uint64_t exp, uint64_t p) {
    if (p == 1) return 0;

    uint64_t result = 1;
    base = base % p;

    while (exp > 0) {
        if (exp & 1) {
            result = cfx_ntt_mod_mul(result, base, p);
        }
        exp >>= 1;
        base = cfx_ntt_mod_mul(base, base, p);
    }
    return result;
}

uint64_t cfx_ntt_mod_inv(uint64_t a, uint64_t p) {
    /* by Fermat's little theorem: a^(-1) = a^(p-2) mod p */
    return cfx_ntt_mod_pow(a, p - 2, p);
}

/*
 * Primitive Root Utilities
 */

uint64_t cfx_ntt_root_of_unity(uint64_t g, uint64_t p, size_t n) {
    /* omega_n = g^((p-1)/n) mod p */
    uint64_t exp = (p - 1) / n;
    return cfx_ntt_mod_pow(g, exp, p);
}

/*
 * Twiddle Factor Precomputation
 */

int cfx_ntt_twiddles_init(cfx_ntt_twiddles_t* tw, size_t n, uint64_t p, uint64_t g) {
    /* n must be power of 2 */
    if (n == 0 || (n & (n - 1)) != 0) {
        return -1;
    }

    tw->forward = (uint64_t*)malloc(n * sizeof(uint64_t));
    tw->inverse = (uint64_t*)malloc(n * sizeof(uint64_t));
    if (!tw->forward || !tw->inverse) {
        free(tw->forward);
        free(tw->inverse);
        tw->forward = NULL;
        tw->inverse = NULL;
        return -1;
    }

    tw->n = n;
    tw->p = p;
    tw->n_inv = cfx_ntt_mod_inv(n, p);

    /* forward twiddles: omega^k for k = 0..n-1 where omega = g^((p-1)/n) */
    uint64_t omega = cfx_ntt_root_of_unity(g, p, n);
    uint64_t omega_inv = cfx_ntt_mod_inv(omega, p);

    tw->forward[0] = 1;
    tw->inverse[0] = 1;
    for (size_t k = 1; k < n; ++k) {
        tw->forward[k] = cfx_ntt_mod_mul(tw->forward[k - 1], omega, p);
        tw->inverse[k] = cfx_ntt_mod_mul(tw->inverse[k - 1], omega_inv, p);
    }

    return 0;
}

void cfx_ntt_twiddles_free(cfx_ntt_twiddles_t* tw) {
    free(tw->forward);
    free(tw->inverse);
    tw->forward = NULL;
    tw->inverse = NULL;
    tw->n = 0;
}

/*
 * Bit-Reversal Permutation
 */

static size_t reverse_bits(size_t x, int bits) {
    size_t result = 0;
    for (int i = 0; i < bits; ++i) {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

void cfx_ntt_bit_reverse(uint64_t* a, size_t n) {
    if (n <= 1) return;

    int bits = 0;
    size_t tmp = n;
    while (tmp > 1) {
        tmp >>= 1;
        ++bits;
    }

    for (size_t i = 0; i < n; ++i) {
        size_t j = reverse_bits(i, bits);
        if (i < j) {
            uint64_t t = a[i];
            a[i] = a[j];
            a[j] = t;
        }
    }
}

/*
 * Core NTT Transforms (Cooley-Tukey, radix-2, decimation-in-time)
 *
 * Input must be bit-reversed before calling.
 * Twiddles array contains omega^k for k = 0..n-1.
 */

void cfx_ntt_forward(uint64_t* a, size_t n, uint64_t p, const uint64_t* twiddles) {
    if (n <= 1) return;

    /* iterative Cooley-Tukey */
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len >> 1;
        size_t step = n / len;

        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < half; ++j) {
                uint64_t w = twiddles[j * step];
                uint64_t u = a[i + j];
                uint64_t v = cfx_ntt_mod_mul(a[i + j + half], w, p);

                a[i + j] = cfx_ntt_mod_add(u, v, p);
                a[i + j + half] = cfx_ntt_mod_sub(u, v, p);
            }
        }
    }
}

void cfx_ntt_inverse(uint64_t* a, size_t n, uint64_t p,
                     const uint64_t* twiddles, uint64_t n_inv) {
    if (n <= 1) return;

    /* same as forward but with inverse twiddles */
    for (size_t len = 2; len <= n; len <<= 1) {
        size_t half = len >> 1;
        size_t step = n / len;

        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < half; ++j) {
                uint64_t w = twiddles[j * step];
                uint64_t u = a[i + j];
                uint64_t v = cfx_ntt_mod_mul(a[i + j + half], w, p);

                a[i + j] = cfx_ntt_mod_add(u, v, p);
                a[i + j + half] = cfx_ntt_mod_sub(u, v, p);
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        a[i] = cfx_ntt_mod_mul(a[i], n_inv, p);
    }
}

/*
 * Polynomial Convolution via NTT
 *
 * c = a * b (polynomial multiplication) using NTT
 */
int cfx_ntt_convolve(uint64_t* c, const uint64_t* a, size_t len_a,
                     const uint64_t* b, size_t len_b,
                     size_t n, uint64_t p, uint64_t g) {
    if (n == 0 || (n & (n - 1)) != 0) {
        return -1;
    }
    if (len_a == 0 || len_b == 0) {
        for (size_t i = 0; i < n; ++i) {
            c[i] = 0;
        }
        return 0;
    }
    if (n < len_a + len_b - 1) {
        return -1;
    }

    uint64_t* a_ntt = (uint64_t*)malloc(n * sizeof(uint64_t));
    uint64_t* b_ntt = (uint64_t*)malloc(n * sizeof(uint64_t));
    if (!a_ntt || !b_ntt) {
        free(a_ntt);
        free(b_ntt);
        return -1;
    }

    for (size_t i = 0; i < len_a; ++i) {
        a_ntt[i] = a[i] % p;
    }
    for (size_t i = len_a; i < n; ++i) {
        a_ntt[i] = 0;
    }

    for (size_t i = 0; i < len_b; ++i) {
        b_ntt[i] = b[i] % p;
    }
    for (size_t i = len_b; i < n; ++i) {
        b_ntt[i] = 0;
    }

    cfx_ntt_twiddles_t tw;
    if (cfx_ntt_twiddles_init(&tw, n, p, g) != 0) {
        free(a_ntt);
        free(b_ntt);
        return -1;
    }

    cfx_ntt_bit_reverse(a_ntt, n);
    cfx_ntt_forward(a_ntt, n, p, tw.forward);

    cfx_ntt_bit_reverse(b_ntt, n);
    cfx_ntt_forward(b_ntt, n, p, tw.forward);

    for (size_t i = 0; i < n; ++i) {
        c[i] = cfx_ntt_mod_mul(a_ntt[i], b_ntt[i], p);
    }

    cfx_ntt_bit_reverse(c, n);
    cfx_ntt_inverse(c, n, p, tw.inverse, tw.n_inv);

    cfx_ntt_twiddles_free(&tw);
    free(a_ntt);
    free(b_ntt);

    return 0;
}

/* number of 16-bit chunks per limb */
#if CFX_LIMB_BITS == 64
#define CHUNKS_PER_LIMB 4
#elif CFX_LIMB_BITS == 32
#define CHUNKS_PER_LIMB 2
#else
#error "Unsupported CFX_LIMB_BITS"
#endif

size_t cfx_ntt_limbs_to_chunks(uint64_t* chunks, size_t max_chunks,
                                const cfx_limb_t* limbs, size_t n_limbs) {
    size_t n_chunks = n_limbs * CHUNKS_PER_LIMB;
    if (n_chunks > max_chunks) {
        n_chunks = max_chunks;
    }

    size_t idx = 0;
    for (size_t i = 0; i < n_limbs && idx < max_chunks; ++i) {
        cfx_limb_t limb = limbs[i];
        for (int c = 0; c < CHUNKS_PER_LIMB && idx < max_chunks; ++c) {
            chunks[idx++] = limb & CFX_NTT_CHUNK_MASK;
            limb >>= CFX_NTT_CHUNK_BITS;
        }
    }
    return idx;
}

size_t cfx_ntt_chunks_to_limbs(cfx_limb_t* limbs, size_t max_limbs,
                                const uint64_t* chunks, size_t n_chunks) {
    size_t n_limbs = (n_chunks + CHUNKS_PER_LIMB - 1) / CHUNKS_PER_LIMB;
    if (n_limbs > max_limbs) n_limbs = max_limbs;

    for (size_t i = 0; i < n_limbs; ++i) limbs[i] = 0;

    uint64_t carry = 0;
    size_t limb_idx = 0;
    int chunk_in_limb = 0;

    for (size_t i = 0; i < n_chunks; ++i) {
        uint64_t val = chunks[i] + carry;
        uint64_t chunk_val = val & CFX_NTT_CHUNK_MASK;
        carry = val >> CFX_NTT_CHUNK_BITS;

        if (limb_idx < n_limbs) {
            limbs[limb_idx] |= (cfx_limb_t)(chunk_val << (chunk_in_limb * CFX_NTT_CHUNK_BITS));
        }

        chunk_in_limb++;
        if (chunk_in_limb == CHUNKS_PER_LIMB) {
            chunk_in_limb = 0;
            limb_idx++;
        }
    }

    while (carry && limb_idx < n_limbs) {
        if (chunk_in_limb == 0) {
            uint64_t val = carry;
            uint64_t chunk_val = val & CFX_NTT_CHUNK_MASK;
            carry = val >> CFX_NTT_CHUNK_BITS;
            limbs[limb_idx] = (cfx_limb_t)chunk_val;
            chunk_in_limb = 1;
        }
        while (carry && chunk_in_limb < CHUNKS_PER_LIMB && limb_idx < n_limbs) {
            uint64_t val = carry;
            uint64_t chunk_val = val & CFX_NTT_CHUNK_MASK;
            carry = val >> CFX_NTT_CHUNK_BITS;
            limbs[limb_idx] |= (cfx_limb_t)(chunk_val << (chunk_in_limb * CFX_NTT_CHUNK_BITS));
            chunk_in_limb++;
        }
        if (chunk_in_limb == CHUNKS_PER_LIMB) {
            chunk_in_limb = 0;
            limb_idx++;
        }
    }

    while (n_limbs > 0 && limbs[n_limbs - 1] == 0) n_limbs--;

    return n_limbs;
}

static size_t next_pow2(size_t n) {
    if (n == 0) return 1;
    n--;
    n |= n >> 1;
    n |= n >> 2;
    n |= n >> 4;
    n |= n >> 8;
    n |= n >> 16;
#if SIZE_MAX > 0xFFFFFFFF
    n |= n >> 32;
#endif
    return n + 1;
}

size_t cfx_ntt_mul_limbs(cfx_limb_t* out, size_t out_cap,
                          const cfx_limb_t* a, size_t n_a,
                          const cfx_limb_t* b, size_t n_b) {
    if (n_a == 0 || n_b == 0) {
        return 0;
    }

    size_t chunks_a = n_a * CHUNKS_PER_LIMB;
    size_t chunks_b = n_b * CHUNKS_PER_LIMB;
    size_t chunks_result = chunks_a + chunks_b - 1;
    size_t ntt_size = next_pow2(chunks_result);

    uint64_t* ca = (uint64_t*)malloc(ntt_size * sizeof(uint64_t));
    uint64_t* cb = (uint64_t*)malloc(ntt_size * sizeof(uint64_t));
    uint64_t* cc = (uint64_t*)malloc(ntt_size * sizeof(uint64_t));
    if (!ca || !cb || !cc) {
        free(ca);
        free(cb);
        free(cc);
        return 0;
    }

    cfx_ntt_limbs_to_chunks(ca, chunks_a, a, n_a);
    for (size_t i = chunks_a; i < ntt_size; ++i) ca[i] = 0;

    cfx_ntt_limbs_to_chunks(cb, chunks_b, b, n_b);
    for (size_t i = chunks_b; i < ntt_size; ++i) cb[i] = 0;

    int rc = cfx_ntt_convolve(cc, ca, chunks_a, cb, chunks_b,
                               ntt_size, CFX_NTT_P1, CFX_NTT_G1);
    free(ca);
    free(cb);

    if (rc != 0) {
        free(cc);
        return 0;
    }

    size_t result_limbs = cfx_ntt_chunks_to_limbs(out, out_cap, cc, chunks_result);
    free(cc);

    return result_limbs;
}
