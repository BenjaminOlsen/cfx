/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/big.h"
#include "cfx/fac.h"
#include "cfx/algo.h"
#include "cfx/arith.h"
#include "cfx/macros.h"
#include "cfx/primes.h"
#include "cfx/base64.h"
#include "cfx/ntt.h"
#include "cfx/compat.h"

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <stdio.h>

/* Platform-specific includes */
#if defined(_WIN32) || defined(_WIN64)
    #include <intrin.h>  /* MSVC intrinsics for _BitScanReverse, etc. */
#else
    #include <unistd.h>
    #ifdef _POSIX_THREADS
        #define CFX_HAS_PTHREAD 1
        #include <pthread.h>
    #endif
#endif
#include <assert.h>
#include <inttypes.h>



static inline void cfx_big_trim(cfx_big_t* b) {
    while (b->n && b->limb[b->n-1] == 0) --b->n;
    if (b->n == 0 && b->cap){ b->limb[0] = 0; } /* dont trim last zero */
}


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
}

int cfx_big_copy(cfx_big_t* dst, const cfx_big_t* src) {
    if (dst == src) return 0;

    if (src->n == 0) {
        dst->n = 0;
        if (dst->cap) dst->limb[0] = 0;
        return 0;
    }

    /* fast path - dst already has enough capacity */
    if (dst->cap >= src->n) {
        memcpy(dst->limb, src->limb, src->n * sizeof(*src->limb));
        dst->n = src->n;
        if (dst->n < dst->cap) dst->limb[dst->n] = 0;
        return 0;
    }

    /*
     * Slow path: allocate new storage without clobbering dst on OOM.
     * Build in tmp, then swap (commit), then free old dst through tmp.
     */
    cfx_big_t tmp;
    cfx_big_init(&tmp);

    cfx_big_reserve(&tmp, src->n);
    memcpy(tmp.limb, src->limb, src->n * sizeof(*src->limb));
    tmp.n = src->n;

    cfx_big_swap(&tmp, dst);
    cfx_big_free(&tmp);
    return 0;
}

void cfx_big_assign(cfx_big_t* dst, const cfx_big_t* src) {
    if (dst == src) return;
    if (src->n == 0) {
        dst->n = 0;
        if (dst->cap) memset(dst->limb, 0, dst->cap * sizeof(cfx_limb_t));
    }
    cfx_big_reserve(dst, src->n);
    memcpy(dst->limb, src->limb, src->n * sizeof(cfx_limb_t));
    dst->n = src->n;
    cfx_big_trim(dst);
}

void cfx_big_assign_sm(cfx_big_t* dst, const cfx_limb_t src) {
    if (src == 0) {
        dst->n = 0;
        if (dst->cap) memset(dst->limb, 0, dst->cap * sizeof(cfx_limb_t));
    }

    cfx_big_reserve(dst, 1); /* todo: error if != 0 */

    memset(dst->limb, 0, dst->cap * sizeof(cfx_limb_t));
    dst->n = 1;
    dst->limb[0] = src;
}

void cfx_big_assign_zero(cfx_big_t* b) {
    b->n = 0;
    if (b->cap == 0) cfx_big_reserve(b, 1);
    b->limb[0] = 0;
}

void cfx_big_assign_one(cfx_big_t* b) {
    if (b->cap == 0) cfx_big_reserve(b, 1);
    b->n = 1;
    b->limb[0] = 1;
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

int cfx_big_is_one(const cfx_big_t* b) {
    return (b->n == 1 && b->limb[0] == 1);
}

int cfx_big_is_even(const cfx_big_t* b) {
    /* zero is even (0 = 2*0), otherwise check least significant bit */
    return (b->n == 0) || !(b->limb[0] & 0x1);
}

int cfx_big_eq_u64(const cfx_big_t* b, cfx_limb_t n) {
    return (n == 0 && cfx_big_is_zero(b)) || (b->n == 1 && b->limb[0] == n);
}

int cfx_big_eq(const cfx_big_t* a, const cfx_big_t* b) {
    if (a->n != b->n) return 0;
    cfx_limb_t diff = 0;
    for (size_t i = 0; i < a->n; ++i) {
        diff |= (a->limb[i] ^ b->limb[i]);
    }
    return diff == 0;
}

/*
static inline size_t _nz_len(const cfx_big_t* x) {
    size_t n = x->n;
    while (n && x->limb[n-1] == 0) --n;
    return n;
}
*/

/** compare two bigs -
 * returns:
 * -1 if a < b
 * 0 if a == b
 * 1 if a > b
 *
 * neither a nor b can have leading zeros!
 * */
int cfx_big_cmp(const cfx_big_t* a, const cfx_big_t* b) {
    if (a->n != b->n) return (a->n < b->n) ? -1 : 1;
    for (size_t i = a->n; i-- > 0; ) {
        if (a->limb[i] != b->limb[i]) return (a->limb[i] < b->limb[i]) ? -1 : 1;
    }
    return 0;
}

/** compare a big and a small: -
 * returns:
 * -1 if a < b
 * 0 if a == b
 * 1 if a > b
 */
int cfx_big_cmp_sm(const cfx_big_t* a, cfx_limb_t b) {
    if (a->n == 0) /* a == 0 */ return b == 0 ? 0 : -1;
    if (a->n > 1) return 1;
    if (a->limb[0] != b) return (a->limb[0] < b) ? -1 : 1;
    return 0;
}

void cfx_big_swap(cfx_big_t* a, cfx_big_t* b) {
    if (a == b) return;
    cfx_big_t tmp = *a;
    *a = *b;
    *b = tmp;
}

/* ---- Bitwise operations ---- */

typedef enum { BITOP_AND, BITOP_OR, BITOP_XOR } bitop_t;

static inline cfx_limb_t bitop_apply(cfx_limb_t x, cfx_limb_t y, bitop_t op) {
    switch (op) {
        case BITOP_AND: return x & y;
        case BITOP_OR:  return x | y;
        case BITOP_XOR: return x ^ y;
    }
    return 0;
}

/*
 * Shared helper for AND, OR, XOR.
 * - AND: result size = min(a.n, b.n), high limbs vanish
 * - OR/XOR: result size = max(a.n, b.n), high limbs copied from larger
 */
static void big_bitop(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, bitop_t op) {
    /* XOR with self is zero */
    if (op == BITOP_XOR && a == b) {
        cfx_big_assign_zero(out);
        return;
    }

    const size_t an = a->n, bn = b->n;
    const size_t min_n = (an < bn) ? an : bn;
    const size_t max_n = (an > bn) ? an : bn;
    const cfx_big_t* lg = (an >= bn) ? a : b;

    /* AND shrinks to min, OR/XOR grow to max */
    const size_t out_n = (op == BITOP_AND) ? min_n : max_n;

    /* Handle aliasing: copy inputs if out overlaps */
    cfx_big_t tmp_a, tmp_b;
    int alias_a = (out == a), alias_b = (out == b);
    if (alias_a) { cfx_big_init(&tmp_a); cfx_big_copy(&tmp_a, a); a = &tmp_a; }
    if (alias_b) { cfx_big_init(&tmp_b); cfx_big_copy(&tmp_b, b); b = &tmp_b; }

    cfx_big_reserve(out, out_n);
    out->n = out_n;

    /* Apply op to overlapping limbs */
    for (size_t i = 0; i < min_n; ++i)
        out->limb[i] = bitop_apply(a->limb[i], b->limb[i], op);

    /* For OR/XOR, copy remaining limbs from the larger operand */
    if (op != BITOP_AND) {
        for (size_t i = min_n; i < max_n; ++i)
            out->limb[i] = lg->limb[i];
    }

    cfx_big_trim(out);

    if (alias_a) cfx_big_free(&tmp_a);
    if (alias_b) cfx_big_free(&tmp_b);
}

void cfx_big_and(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    big_bitop(out, a, b, BITOP_AND);
}

void cfx_big_or(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    big_bitop(out, a, b, BITOP_OR);
}

void cfx_big_xor(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    big_bitop(out, a, b, BITOP_XOR);
}



/* In-place variants */
void cfx_big_and_eq(cfx_big_t* a, const cfx_big_t* b) {
    big_bitop(a, a, b, BITOP_AND);
}

void cfx_big_or_eq(cfx_big_t* a, const cfx_big_t* b) {
    big_bitop(a, a, b, BITOP_OR);
}

void cfx_big_xor_eq(cfx_big_t* a, const cfx_big_t* b) {
    big_bitop(a, a, b, BITOP_XOR);
}

void cfx_big_rotl_w(cfx_big_t* out, const cfx_big_t* b, unsigned r, unsigned w) {
    if (w == 0) {
        cfx_big_from_limb(out, 0);
        return;
    }

    unsigned rot = r % w;
    if (rot == 0) {
        cfx_big_copy(out, b);
        cfx_big_mask_bits(out, w);
        return;
    }

    cfx_big_t hi, lo;
    cfx_big_init(&hi);
    cfx_big_init(&lo);
    cfx_big_shl_bits(&hi, b, rot);
    cfx_big_shr_bits(&lo, b, w - rot);
    cfx_big_or(out, &hi, &lo);
    cfx_big_mask_bits(out, w);
    cfx_big_free(&hi);
    cfx_big_free(&lo);
}

void cfx_big_rotl(cfx_big_t* out, const cfx_big_t* b, unsigned r) {
    cfx_big_rotl_w(out, b, r, cfx_big_bitlen(b));
}

void cfx_big_rotr_w(cfx_big_t* out, const cfx_big_t* b, unsigned r, unsigned w) {
    if (w == 0) {
        cfx_big_from_limb(out, 0);
        return;
    }

    unsigned rot = r % w;
    if (rot == 0) {
        cfx_big_copy(out, b);
        cfx_big_mask_bits(out, w);
        return;
    }

    cfx_big_t hi, lo;
    cfx_big_init(&hi);
    cfx_big_init(&lo);
    cfx_big_shr_bits(&lo, b, rot);
    cfx_big_shl_bits(&hi, b, w - rot);
    cfx_big_or(out, &hi, &lo);
    cfx_big_mask_bits(out, w);
    cfx_big_free(&hi);
    cfx_big_free(&lo);
}

void cfx_big_rotr(cfx_big_t* out, const cfx_big_t* b, unsigned r) {
    cfx_big_rotr_w(out, b, r, cfx_big_bitlen(b));
}

void cfx_big_mask_bits(cfx_big_t* a, unsigned nbits) {
    if (nbits == 0) {
        cfx_big_assign_zero(a);
        return;
    }
    unsigned full_limbs = nbits / CFX_LIMB_BITS;
    unsigned rem_bits   = nbits % CFX_LIMB_BITS;

    /* all bits fit in existing limbs */
    if (full_limbs >= a->n) return;

    /* truncate to full_limbs (+1 if partial) */
    if (rem_bits == 0) {
        a->n = full_limbs;
    } else {
        a->n = full_limbs + 1;
        cfx_limb_t mask = ((cfx_limb_t)1 << rem_bits) - 1;
        a->limb[full_limbs] &= mask;
    }
    cfx_big_trim(a);
}

int cfx_big_endswith_u64(const cfx_big_t* x, uint64_t value) {
    if (!x || x->n == 0) {
        return value == 0;
    }

#if CFX_LIMB_BITS >= 64
    /* Single limb covers 64 bits */
    return x->limb[0] == (cfx_limb_t)value;

#else
    /* 32-bit limbs: need up to two limbs */
    uint64_t lo = 0;

    if (x->n >= 1) {
        lo |= (uint64_t)x->limb[0];
    }
    if (x->n >= 2) {
        lo |= ((uint64_t)x->limb[1] << 32);
    }

    return lo == value;
#endif
}

int cfx_big_bit_is_set(const cfx_big_t* x, size_t bit) {
    if (!x || x->n == 0) return 0;

    const size_t limb_idx = bit / CFX_LIMB_BITS;
    const size_t bit_idx  = bit % CFX_LIMB_BITS;

    if (limb_idx >= x->n) return 0;

    return (x->limb[limb_idx] >> bit_idx) & 1;
}

void cfx_big_bit_set(cfx_big_t* x, size_t bit) {
    const size_t limb_idx = bit / CFX_LIMB_BITS;
    const size_t bit_idx  = bit % CFX_LIMB_BITS;

    /* Grow if needed */
    if (limb_idx >= x->n) {
        cfx_big_reserve(x, limb_idx + 1);
        for (size_t i = x->n; i <= limb_idx; ++i)
            x->limb[i] = 0;
        x->n = limb_idx + 1;
    }

    x->limb[limb_idx] |= (cfx_limb_t)1 << bit_idx;
}

void cfx_big_bit_clear(cfx_big_t* x, size_t bit) {
    const size_t limb_idx = bit / CFX_LIMB_BITS;
    const size_t bit_idx  = bit % CFX_LIMB_BITS;

    if (limb_idx >= x->n) return; /* already zero */

    x->limb[limb_idx] &= ~((cfx_limb_t)1 << bit_idx);
    cfx_big_trim(x);
}

void cfx_big_bit_flip(cfx_big_t* x, size_t bit) {
    const size_t limb_idx = bit / CFX_LIMB_BITS;
    const size_t bit_idx  = bit % CFX_LIMB_BITS;

    /* Grow if needed (flipping a zero bit to one) */
    if (limb_idx >= x->n) {
        cfx_big_reserve(x, limb_idx + 1);
        for (size_t i = x->n; i <= limb_idx; ++i)
            x->limb[i] = 0;
        x->n = limb_idx + 1;
    }

    x->limb[limb_idx] ^= (cfx_limb_t)1 << bit_idx;
    cfx_big_trim(x);
}

size_t cfx_big_popcount(const cfx_big_t* x) {
    if (!x || x->n == 0) return 0;

    size_t count = 0;
    for (size_t i = 0; i < x->n; ++i) {
        cfx_limb_t v = x->limb[i];
        /* Portable popcount via bit manipulation */
        while (v) {
            count += v & 1;
            v >>= 1;
        }
    }
    return count;
}

/*
 * Constant-time (maybe) conditional swap: if condition != 0, swap a and b.
 */
void cfx_big_cswap(cfx_big_t* a, cfx_big_t* b, int condition) {
    if (a == b) return;

    /* create all ones or all-zeros mask (without triggering MSVC warning) */
    const cfx_limb_t mask = (cfx_limb_t)0 - (cfx_limb_t)(condition != 0);
    const size_t max_n = (a->n > b->n) ? a->n : b->n;

    cfx_big_reserve(a, max_n);
    cfx_big_reserve(b, max_n);

    for (size_t i = a->n; i < max_n; ++i) a->limb[i] = 0;
    for (size_t i = b->n; i < max_n; ++i) b->limb[i] = 0;

    for (size_t i = 0; i < max_n; ++i) {
        cfx_limb_t diff = mask & (a->limb[i] ^ b->limb[i]);
        a->limb[i] ^= diff;
        b->limb[i] ^= diff;
    }

    size_t n_diff = (size_t)(mask & (cfx_limb_t)(a->n ^ b->n));
    a->n ^= n_diff;
    b->n ^= n_diff;

    cfx_big_trim(a);
    cfx_big_trim(b);
}

size_t cfx_big_bitlen(const cfx_big_t* b) {
    /* Returns 0 for zero; otherwise assumes no leading zero limbs */
    if (b->n == 0) return 0;
    const size_t limb_bits = 8u * sizeof(cfx_limb_t);
    cfx_limb_t v = b->limb[b->n - 1];
    size_t lz = cfx_clz(v);
    return b->n * limb_bits - lz;
}

void cfx_big_reserve(cfx_big_t* b, size_t need) {
    if (need <= b->cap) return;

    size_t old_cap = b->cap;
    size_t new_cap = old_cap ? old_cap : 32;

    while (new_cap < need) {
        if (new_cap > SIZE_MAX / 2) {  /* doubling would overflow */
            new_cap = need;            /* fall back to exactly need */
            break;
        }
        new_cap *= 2;
    }

    if (new_cap > SIZE_MAX / sizeof(cfx_limb_t)) {
        fprintf(stderr, "cfx_big_reserve: allocation size overflow\n");
        abort();
    }

    size_t new_bytes = new_cap * sizeof(cfx_limb_t);

    void* tmp = realloc(b->limb, new_bytes);
    if (!tmp) {
        fprintf(stderr, "cfx_big_reserve: out of memory (requested %zu bytes)\n", new_bytes);
        abort();
    }

    b->limb = (cfx_limb_t*)tmp;

    /* Zero the newly added region only (realloc doesn't do this for you) */
    if (new_cap > old_cap) {
        size_t add = new_cap - old_cap;
        memset(b->limb + old_cap, 0, add * sizeof(cfx_limb_t));
    }

    b->cap = new_cap;
}


int cfx_big_from_u64(cfx_big_t* b, uint64_t v) {
    cfx_limb_t limbs[2];

#if (CFX_LIMB_BITS == 64)

    limbs[0] = (cfx_limb_t)v;
    cfx_big_from_limbs(b, limbs, 1);

#elif (CFX_LIMB_BITS == 32)

    limbs[0] = (cfx_limb_t)(uint32_t)(v & 0xffffffffu);
    limbs[1] = (cfx_limb_t)(uint32_t)(v >> 32);

    if (limbs[1] != 0) {
        cfx_big_from_limbs(b, limbs, 2);
    } else {
        cfx_big_from_limbs(b, limbs, 1);
    }

#else
#   error "Unsupported CFX_LIMB_BITS"
#endif
    return 0;
}

/* assumes b is already initted */
int cfx_big_from_limb(cfx_big_t* b, cfx_limb_t v) {
    if (v == 0) {
        cfx_big_assign_zero(b);
        return 0;
    }

    if (b->cap < 1) {
        cfx_big_reserve(b, 1);
    }

    b->limb[0] = v;

    /* If the old value had more limbs, wipe them for hygiene */
    if (b->n > 1) {
        memset(b->limb + 1, 0, (b->n - 1) * sizeof(*b->limb));
    }

    b->n = 1;

    /* cheap invariant */
    /* CFX_ASSERT(b->n == 1 && b->limb[0] != 0); */
    return 0;
}

/* out == NULL is allowed, as a size query;
 returns 0 on success, -1 buf too small*/
int cfx_big_to_bytes_be(uint8_t* out, size_t* outlen, const cfx_big_t* b) {

    if (b->n == 0) {
        if (outlen) *outlen = 0;
        return 0;
    }

    cfx_limb_t limb = b->limb[b->n-1];
    size_t bits = cfx_big_bitlen(b);
    size_t bytes = (bits + 7) / 8;

    if (!out) {
        *outlen = bytes;
        return 0;
    }

    size_t o = bytes;
    if (*outlen < bytes) return -1;
    for (size_t i = 0; i < b->n; ++i) {
        limb = b->limb[i];
        for (size_t k = 0; k < sizeof(cfx_limb_t); ++k) {
            if (o == 0) break;
            out[--o] = (uint8_t)(limb & 0xFF);
            limb >>= 8;
        }
    }

    *outlen = bytes;
    return 0;
}

/* copy big endian bytes into cfx_big_t. */
/* - be[0] is the most-significant byte */
/* - trims any leading zero bytes */
/* Returns 0 on success, nonzero on allocation/argument error. */
int cfx_big_from_bytes_be(cfx_big_t* out, const uint8_t* be, size_t len) {
    if (!out || (!be && len)) return -1;

    size_t off = 0;
    while (off < len && be[off] == 0) ++off;

    if (off == len) {
        /* value == 0 */
        cfx_big_reserve(out, 1);
        out->n = 0;
        out->limb[0] = 0;
        return 0;
    }

    const size_t nbytes = len - off;
    const size_t lb     = sizeof(cfx_limb_t);
    const size_t nlimbs = (nbytes + lb - 1) / lb;

    cfx_big_reserve(out, nlimbs);

    /* pack bytes into (little-endian) limbs */
    /* For limb i (0 = least significant), we read up to lb bytes from the end. */
    size_t src_end = len; /* one-past-the-last valid index */
    for (size_t i = 0; i < nlimbs; ++i) {
        cfx_limb_t limb = 0;
        size_t     take = (src_end > off) ? ((src_end - off) < lb ? (src_end - off) : lb) : 0;

        /* copy 'take' bytes into this limb, least significant byte last in BE stream */
        for (size_t j = 0; j < take; ++j) {
            uint8_t b = be[src_end - 1 - j];
            limb |= (cfx_limb_t)b << (8u * j);
        }

        out->limb[i] = limb;
        src_end  -= take;
    }

    out->n = nlimbs;
    cfx_big_trim(out);
    return 0;
}


static inline void _mul_sm_fast(cfx_big_t* b, cfx_limb_t m) {
    size_t n = b->n;
    cfx_big_reserve(b, n + 1);
    cfx_limb_t* p = b->limb;

#if (CFX_USE_X86_INTRINSICS == 1) && defined(__BMI2__) && (CFX_LIMB_BITS == 64)
    cfx_limb_t carry = 0;
    for (size_t i = 0; i < n; ++i) {
        unsigned long long lo, hi;

        /* lo = (p[i] * m) low 64, hi = high 64 — flags NOT clobbered */
        lo = _mulx_u64(p[i], m, &hi);

        /* add previous carry into lo, get carry-out in c1 */
        unsigned char c1 = _addcarry_u64(0, lo, carry, &lo);
        hi += c1;
        p[i] = (cfx_limb_t)lo;
        carry = (cfx_limb_t)hi;
    }

    p[n] = carry;
    b->n = n + (carry != 0);
#else /* Portable fallback (cfx_acc_t) */
    cfx_limb_t carry = 0;
    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t t;
        cfx_acc_mul(&t, p[i], m);
        cfx_acc_add_lo(&t, carry);
        p[i] = cfx_acc_lo(t);
        carry = cfx_acc_hi(t);
    }
    p[n] = carry;
    b->n = n + (carry != 0);
#endif
}

/* Multiply by p^e by repeated squaring using small chunks to avoid u32 overflow */
void cfx_big_expmul_prime(cfx_big_t* b, cfx_limb_t p, cfx_limb_t e) {

    /* Find largest t so that p^t fits in 32 bits -> p^2t fits in 64 */
    cfx_limb_t t = 1;
    cfx_acc_t acc = cfx_acc_from_lo(p);
    const cfx_acc_t lim = cfx_acc_from_lo(CFX_SQRT_ACC_MAX);

    cfx_acc_t quot, rem;
    cfx_acc_divrem(&quot, &rem, &lim, &acc);
    while (cfx_acc_le(&acc, &quot)) {
        cfx_acc_mul_eq(&acc, &acc);
        t *= 2u;
        cfx_acc_divrem(&quot, &rem, &lim, &acc);
    }
    /* now t is the largest s.t. p^t <= lim */
    CFX_PRINT_DBG("multiplying by " CFX_PRIuLIMB "^" CFX_PRIuLIMB " by breaking it into (" CFX_PRIuLIMB "^" CFX_PRIuLIMB ")^(" CFX_PRIuLIMB "/" CFX_PRIuLIMB ") \n", p, e, p, t, e, t);

    /* Now multiply by (p^t)^(e/t) and then the remainder */
    /* Compute p^t */
    /* compute power safely */
    cfx_limb_t pow_t = p;

    /* fast power to get pow_t = p^t by using binary expansion of t:
    p^t = p^(b_i*2^i) * p^(b_(i-1)*2^(i-1) * ... */
    pow_t = 1;
    acc = p;
    cfx_limb_t tt = t;
    while (tt) {
        if (tt & 1u) pow_t = (cfx_limb_t)(pow_t * acc);
        tt >>= 1u;
        acc = acc*acc;
    }
    cfx_limb_t q = e / t;
    cfx_limb_t r = e % t;

    for (cfx_limb_t i = 0; i < q; i++) _mul_sm_fast(b, pow_t);

    /* multiply remainder by binary exponentiation (still fits u32) */
    cfx_limb_t rempow = 1;
    acc = p;
    cfx_limb_t rr = r;

    while (rr) {
        if (rr & 1u) rempow = (cfx_limb_t)(rempow * acc);
        rr >>= 1u;
        acc = acc*acc;
    }
    if (rempow != 1) _mul_sm_fast(b, rempow);
}

static inline void cfx_big_mul_small_inplace(cfx_big_t* out, cfx_limb_t m) {
    if (m == 0) { cfx_big_from_limb(out, 0); return; }
    if (m == 1) return;
    _mul_sm_fast(out, m);
}

/* out = p^e (p is a limb prime) */
void cfx_big_pow_sm(cfx_big_t* out, cfx_limb_t p, cfx_limb_t e) {
    cfx_big_from_limb(out, 1);
    if (e == 0) return;


    if (p == 2) {
        /* out = 1 << e */
        cfx_big_t one;
        cfx_big_init(&one);
        cfx_big_from_limb(&one, 1);
        cfx_big_shl_bits(out, &one, (unsigned)e);
        cfx_big_free(&one);
        return;
    }

    if (e <= 3) {
        for (cfx_limb_t i = 0; i < e; i++) cfx_big_mul_small_inplace(out, p);
        return;
    }

    /* Binary exponentiation */
    cfx_big_t base;
    cfx_big_init(&base);
    cfx_big_from_limb(&base, p);

    cfx_limb_t exp = e;
    while (exp) {
        if (exp & 1u) {
            cfx_big_mul_auto(out, &base); /* out *= base */
        }
        exp >>= 1u;
        if (exp) {
            /* base = base^2 */
            cfx_big_mul_auto(&base, &base);
        }
    }

    cfx_big_free(&base);
}

void cfx_big_exp(cfx_big_t* out, const cfx_big_t* n, const cfx_big_t* p) {
    if (cfx_big_is_zero(p))  { cfx_big_from_limb(out, 1); return; }
    if (cfx_big_is_zero(n))  { PRINT_BIG(">>>>>>>>>> n is zero!", n); cfx_big_from_limb(out, 0); return; }
    if (cfx_big_eq_u64(n, 1)) { cfx_big_from_limb(out, 1); return; }

    cfx_big_t acc, pp, np; /* accumulator, p copy, n^p*/
    cfx_big_init(&acc);
    cfx_big_init(&pp);
    cfx_big_init(&np);

    cfx_big_from_limb(&np, 1);
    cfx_big_copy(&pp, p);
    cfx_big_copy(&acc, n);

    while (!cfx_big_is_zero(&pp)) {
        if (pp.n && (pp.limb[0] & 1)) {
            cfx_big_mul_auto(&np, &acc);
        }
        cfx_big_shr_bits(&pp, &pp, 1);
        if (!cfx_big_is_zero(&pp)) cfx_big_mul_auto(&acc, &acc);
    }
    cfx_big_move(out, &np);
    cfx_big_free(&np);
    cfx_big_free(&pp);
    cfx_big_free(&acc);
}

void cfx_big_exp_u64(cfx_big_t* out, const cfx_big_t* n, cfx_limb_t p) {
    if (p == 0)              { cfx_big_from_limb(out, 1); return; }
    if (cfx_big_is_zero(n))  { cfx_big_from_limb(out, 0); return; }
    if (cfx_big_eq_u64(n, 1)) { cfx_big_from_limb(out, 1); return; }

    cfx_big_t acc, np; /* accumulator, p copy, n^p*/
    cfx_big_init(&acc);
    cfx_big_init(&np);
    cfx_big_from_limb(&np, 1);
    cfx_big_copy(&acc, n);
    while (p) {
        if (p & 1) {
            cfx_big_mul_auto(&np, &acc);
        }
        p >>= 1;
        if (p) cfx_big_mul_auto(&acc, &acc);
    }
    cfx_big_move(out, &np);
    cfx_big_free(&np);
    cfx_big_free(&acc);
}

/* if p is secret, this leaks info */
void cfx_big_mod_exp(cfx_big_t* out, const cfx_big_t* n, const cfx_big_t* p, const cfx_big_t* m) {
    if (cfx_big_is_zero(p))  { cfx_big_from_limb(out, 1); return; }
    if (cfx_big_is_zero(n))  { cfx_big_from_limb(out, 0); return; }
    if (cfx_big_eq_u64(n, 1)) { cfx_big_from_limb(out, 1); return; }

    cfx_big_t acc, pp, np; /* accumulator, p copy, n^p*/
    cfx_big_init(&acc);
    cfx_big_init(&pp);
    cfx_big_init(&np);
    cfx_big_from_limb(&np, 1);
    cfx_big_copy(&pp, p);
    cfx_big_copy(&acc, n);
    while (!cfx_big_is_zero(&pp)) {
        if (pp.n && (pp.limb[0] & 1)) {
            cfx_big_mul_auto(&np, &acc);
            cfx_big_mod(&np, &np, m);
        }
        cfx_big_shr_bits(&pp, &pp, 1);
        if (!cfx_big_is_zero(&pp)) cfx_big_mul_auto(&acc, &acc);
    }
    cfx_big_move(out, &np);
    cfx_big_free(&np);
    cfx_big_free(&pp);
    cfx_big_free(&acc);
}


/* out = (a^e) mod m */
static void powmod_u64_base(cfx_big_t* out, cfx_limb_t a, const cfx_big_t* e, const cfx_big_t* mod) {
    /* out = a^e mod mod, with small base a and big exponent e */
    cfx_big_t base, res, ee;
    cfx_big_init(&base);
    cfx_big_init(&res);
    cfx_big_init(&ee);
    cfx_big_assign_sm(&base, a);
    cfx_big_mod(&base, &base, mod);

    cfx_big_assign_sm(&res, 1);
    cfx_big_copy(&ee, e);

    /* square & multiply */
    while (!cfx_big_is_zero(&ee)) {
        if (!cfx_big_is_even(&ee)) {  /* if (ee.limb[0] & 1) */
            cfx_big_mulmod(&res, &res, &base, mod);
        }
        cfx_big_mulmod(&base, &base, &base, mod);
        cfx_big_shr_bits_eq(&ee, 1);
    }
    cfx_big_copy(out, &res);
    cfx_big_free(&ee);
    cfx_big_free(&res);
    cfx_big_free(&base);

}

/* Single Miller–Rabin round.
 * n must be odd > 2.  n-1 = d * 2^s, with d odd.
 * a is a small base (2,3,5,7,...) < 2^64.
 */
static int miller_rabin_once(const cfx_big_t* n, cfx_limb_t a,
                             const cfx_big_t* d, cfx_limb_t s) {
    cfx_big_t x, nm1;
    int rc = 0;

    cfx_big_init(&x);
    cfx_big_init(&nm1);

    /* x = a^d mod n */
    powmod_u64_base(&x, a, d, n);

    /* nm1 = n-1 */
    cfx_big_copy(&nm1, n);
    cfx_big_sub_sm_eq(&nm1, 1);

    if (cfx_big_is_one(&x)) {
        rc = 1;
        goto cleanup;
    }
    if (cfx_big_cmp(&x, &nm1) == 0) {
        rc = 1;
        goto cleanup;
    }

    for (cfx_limb_t r = 1; r < s; ++r) {
        cfx_big_mulmod(&x, &x, &x, n);
        if (cfx_big_cmp(&x, &nm1) == 0) {
            rc = 1;
            goto cleanup;
        }
    }

    rc = 0;     /* witness: composite */

cleanup:
    cfx_big_free(&x);
    cfx_big_free(&nm1);
    return rc;
}



/* Does miller - rabin's primality test on a big */
int cfx_big_is_prime(const cfx_big_t* n) {

    if (cfx_big_cmp_sm(n, 2) < 0)
        return 0;

    static const uint32_t small_primes[] = {
        2,3,5,7,11,13,17,19,23,29,31,37,
        41,43,47,53,59,61,67,71,73,79,83,89,97,0
    };

    /* Small trial division first */
    for (int i = 0; small_primes[i]; ++i) {
        uint32_t p = small_primes[i];
        int cmp = cfx_big_cmp_sm(n, p);
        if (cmp == 0)
            return 1;
        if (cmp > 0 && cfx_big_mod_sm(n, p) == 0)
            return 0;
    }

    if (cfx_big_is_even(n))
        return 0;

    /* Write n-1 = d * 2^s with d odd */
    cfx_big_t d;
    cfx_big_init(&d);
    cfx_limb_t s = 0;

    cfx_big_copy(&d, n);
    cfx_big_sub_sm_eq(&d, 1);          /* d = n-1 */

    while (!cfx_big_is_zero(&d) && cfx_big_is_even(&d)) {
        cfx_big_shr_bits_eq(&d, 1);
        ++s;
    }

    if (cfx_big_is_zero(&d)) {
        cfx_big_free(&d);
        return 0;
    }

    /* for numbers up to 64 bits, use the "Sinclair 7" bases which are
       deterministic for all 64-bit inputs. For larger numbers, use
       additional small prime bases for higher confidence. */
    size_t bitlen = cfx_big_bitlen(n);

    if (bitlen <= 64) {
        static const uint64_t sinclair7[] = {
            2ULL, 325ULL, 9375ULL, 28178ULL, 450775ULL, 9780504ULL, 1795265022ULL, 0
        };
        for (size_t i = 0; sinclair7[i]; ++i) {
            if (!miller_rabin_once(n, (cfx_limb_t)sinclair7[i], &d, s)) {
                cfx_big_free(&d);
                return 0;
            }
        }
    } else {
        /* For larger numbers, use first 12 primes as witnesses.
           This is probabilistic but with very low false positive rate. */
        static const cfx_limb_t bases[] = {
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 0
        };
        for (size_t i = 0; bases[i]; ++i) {
            if (!miller_rabin_once(n, bases[i], &d, s)) {
                cfx_big_free(&d);
                return 0;
            }
        }
    }

    cfx_big_free(&d);
    return 1;
}

/* Binary GCD algorithm for big integers */
void cfx_big_gcd(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    if (cfx_big_is_zero(a)) {
        cfx_big_copy(out, b);
        return;
    }
    if (cfx_big_is_zero(b)) {
        cfx_big_copy(out, a);
        return;
    }

    cfx_big_t u, v;
    cfx_big_init(&u);
    cfx_big_init(&v);
    cfx_big_copy(&u, a);
    cfx_big_copy(&v, b);

    /* cnt common factors of 2 */
    size_t shift = 0;
    while (cfx_big_is_even(&u) && cfx_big_is_even(&v)) {
        cfx_big_shr_bits_eq(&u, 1);
        cfx_big_shr_bits_eq(&v, 1);
        shift++;
    }

    /* remove remaining factors of 2 from u */
    while (cfx_big_is_even(&u)) {
        cfx_big_shr_bits_eq(&u, 1);
    }

    do {
        /* remove factors of 2 from v */
        while (cfx_big_is_even(&v)) {
            cfx_big_shr_bits_eq(&v, 1);
        }

        /* make sure u <= v */
        if (cfx_big_cmp(&u, &v) > 0) {
            cfx_big_swap(&u, &v);
        }

        cfx_big_sub_eq(&v, &u);
        cfx_big_trim(&v);

    } while (!cfx_big_is_zero(&v));

    /* restore common factors of 2 */
    cfx_big_shl_bits_eq(&u, shift);

    cfx_big_swap(out, &u);
    cfx_big_free(&u);
    cfx_big_free(&v);
}

/* Pollard-Rho factorization using Montgomery multiplication (Brent's improvement).
 * Returns a non-trivial factor in 'factor', or copies n if n is prime/unfactorable. */
void cfx_big_pollard_rho(cfx_big_t* factor, const cfx_big_t* n) {
    if (cfx_big_cmp_sm(n, 2) < 0) {
        cfx_big_copy(factor, n);
        return;
    }
    if (cfx_big_is_even(n)) {
        cfx_big_from_limb(factor, 2);
        return;
    }
    if (cfx_big_cmp_sm(n, 3) == 0) {
        cfx_big_from_limb(factor, 3);
        return;
    }

    /* monty context for modular arithmetic mod n */
    cfx_big_mont_ctx_t ctx;
    if (cfx_big_mont_ctx_init(&ctx, n) == 0) {
        cfx_big_copy(factor, n);
        return;
    }

    cfx_big_t y, c, x, ys, q, g, one, diff, temp;
    cfx_big_init(&y);
    cfx_big_init(&c);
    cfx_big_init(&x);
    cfx_big_init(&ys);
    cfx_big_init(&q);
    cfx_big_init(&g);
    cfx_big_init(&one);
    cfx_big_init(&diff);
    cfx_big_init(&temp);

    /* init one in Montgomery form */
    cfx_big_from_limb(&one, 1);
    cfx_big_mont_to(&one, &one, &ctx);

    /* Try different c values if we fail with the first one */
    static const cfx_limb_t c_vals[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23};
    const size_t num_c_vals = sizeof(c_vals) / sizeof(c_vals[0]);

    for (size_t c_idx = 0; c_idx < num_c_vals; ++c_idx) {
        /* Initialize with new c value */
        cfx_big_from_limb(&y, 2);
        cfx_big_mont_to(&y, &y, &ctx);
        cfx_big_from_limb(&c, c_vals[c_idx]);
        cfx_big_mont_to(&c, &c, &ctx);

        cfx_big_copy(&x, &y);
        cfx_big_copy(&ys, &y);
        cfx_big_copy(&q, &one);
        cfx_big_from_limb(&g, 1);

        size_t nbits = n->n * CFX_LIMB_BITS;
        size_t max_iters = 1ULL << (nbits / 4 + 4);  // e.g., 128-bit -> 2^36
        if (max_iters > (1ULL << 32)) max_iters = 1ULL << 32;  // cap at 4B

        size_t r = 1;
        const size_t batch = 128;
        size_t iters = 0;

        while (cfx_big_is_one(&g)) {
            cfx_big_copy(&x, &y);

            /* advance y by r steps: y = y^2 + c mod n */
            for (size_t i = 0; i < r; ++i) {
                cfx_big_mont_sqr(&y, &y, &ctx);
                cfx_big_add_eq(&y, &c);
                if (cfx_big_cmp(&y, n) >= 0) {
                    cfx_big_sub_eq(&y, n);
                }
            }

            /* GCD computation batch*/
            for (size_t k = 0; k < r && cfx_big_is_one(&g); k += batch) {
                cfx_big_copy(&ys, &y);
                size_t lim = (batch < r - k) ? batch : (r - k);

                for (size_t i = 0; i < lim; ++i) {
                    /* y = y^2 + c mod n */
                    cfx_big_mont_sqr(&y, &y, &ctx);
                    cfx_big_add_eq(&y, &c);
                    if (cfx_big_cmp(&y, n) >= 0) {
                        cfx_big_sub_eq(&y, n);
                    }

                    /* diff = |x - y| */
                    if (cfx_big_cmp(&x, &y) >= 0) {
                        cfx_big_copy(&diff, &x);
                        cfx_big_sub_eq(&diff, &y);
                    } else {
                        cfx_big_copy(&diff, &y);
                        cfx_big_sub_eq(&diff, &x);
                    }

                    /* q = q * diff mod n (in Montgomery form) */
                    if (!cfx_big_is_zero(&diff)) {
                        cfx_big_mont_mul(&q, &q, &diff, &ctx);
                    }
                }

                /* g = gcd(q, n) */
                cfx_big_mont_from(&temp, &q, &ctx);
                cfx_big_gcd(&g, &temp, n);
            }

            r <<= 1;
            iters += r;
            if (iters > max_iters) {
                /* Try next c value */
                break;
            }
        }

        /* If g == n, backtrack step-by-step */
        if (cfx_big_cmp(&g, n) == 0) {
            cfx_big_from_limb(&g, 1);
            cfx_big_mont_from(&x, &x, &ctx);  /* Convert x back from Montgomery */

            for (size_t i = 0; i < batch && cfx_big_is_one(&g); ++i) {
                /* ys = ys^2 + c mod n */
                cfx_big_mont_sqr(&ys, &ys, &ctx);
                cfx_big_add_eq(&ys, &c);
                if (cfx_big_cmp(&ys, n) >= 0) {
                    cfx_big_sub_eq(&ys, n);
                }

                /* Convert ys from Montgomery form for GCD */
                cfx_big_mont_from(&temp, &ys, &ctx);

                /* diff = |x - ys| */
                if (cfx_big_cmp(&x, &temp) >= 0) {
                    cfx_big_copy(&diff, &x);
                    cfx_big_sub_eq(&diff, &temp);
                } else {
                    cfx_big_copy(&diff, &temp);
                    cfx_big_sub_eq(&diff, &x);
                }

                cfx_big_gcd(&g, &diff, n);
            }
        }

        /* Check if we found a non-trivial factor */
        if (cfx_big_cmp(&g, n) != 0 && !cfx_big_is_one(&g)) {
            cfx_big_copy(factor, &g);
            goto cleanup;
        }
    }

    /* All c values exhausted - return n (failure) */
    cfx_big_copy(factor, n);

cleanup:
    cfx_big_free(&y);
    cfx_big_free(&c);
    cfx_big_free(&x);
    cfx_big_free(&ys);
    cfx_big_free(&q);
    cfx_big_free(&g);
    cfx_big_free(&one);
    cfx_big_free(&diff);
    cfx_big_free(&temp);
    cfx_big_mont_ctx_free(&ctx);
}


void cfx_big_sq_eq(cfx_big_t* b) {

    const size_t n = b->n;
    cfx_big_t ret;
    cfx_big_init(&ret);

    if (n == 0) {
        return;
    }

    cfx_big_reserve(&ret, 2*n);
    memset(ret.limb, 0, 2*n * sizeof(cfx_limb_t));
    ret.n = 2*n;

    /* 1) Cross terms: for i < j, add 2*b[i]*b[j] into ret[i+j] */
    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t carry = 0;
        for (size_t j = i + 1; j < n; ++j) {
            cfx_acc_t p = (cfx_acc_t)b->limb[i] * b->limb[j];

            /* add p once */
            cfx_acc_t t = (cfx_acc_t)ret.limb[i + j]
                            + (cfx_limb_t)p
                            + (cfx_limb_t)carry;
            ret.limb[i + j] = (cfx_limb_t)t;
            carry = (carry >> CFX_LIMB_BITS) + (t >> CFX_LIMB_BITS) + (p >> CFX_LIMB_BITS);

            /* add p again (to double) -- same carry rule */
            t = (cfx_acc_t)ret.limb[i + j]
                + (cfx_limb_t)p
                + (cfx_limb_t)carry;
            ret.limb[i + j] = (cfx_limb_t)t;
            carry = (carry >> CFX_LIMB_BITS) + (t >> CFX_LIMB_BITS) + (p >> CFX_LIMB_BITS);
        }

        /* propagate whatever is left in 'carry' */
        size_t k = i + n; /* next column after the last updated (i + (n-1)) */
        while (carry) {
            cfx_acc_t t = (cfx_acc_t)ret.limb[k] + (cfx_limb_t)carry;
            ret.limb[k] = (cfx_limb_t)t;
            carry = (carry >> CFX_LIMB_BITS) + (t >> CFX_LIMB_BITS);
            ++k;
        }
    }

    /* 2) diagonals: add b[i]^2 once at ret[2*i] */
    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t sq = (cfx_acc_t)b->limb[i] * b->limb[i];

        cfx_acc_t t = (cfx_acc_t)ret.limb[2*i] + (cfx_limb_t)sq;
        ret.limb[2*i] = (cfx_limb_t)t;
        cfx_acc_t c = (t >> CFX_LIMB_BITS) + (sq >> CFX_LIMB_BITS);

        size_t k = 2*i + 1;
        while (c) {
            t = (cfx_acc_t)ret.limb[k] + (cfx_limb_t)c;
            ret.limb[k] = (cfx_limb_t)t;
            c = (c >> CFX_LIMB_BITS) + (t >> CFX_LIMB_BITS);
            ++k;
        }
    }

    /* trim */
    while (ret.n && ret.limb[ret.n - 1] == 0) --ret.n;
    cfx_big_swap(&ret, b);
    cfx_big_free(&ret);
}

void cfx_big_mul_fft(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    if (cfx_big_is_zero(a) || cfx_big_is_zero(b)) {
        cfx_big_from_limb(out, 0);
        return;
    }

    const size_t na = a->n;
    const size_t nb = b->n;
    const size_t nout = na + nb;

    cfx_big_reserve(out, nout);

    size_t n = cfx_ntt_mul_limbs(out->limb, nout, a->limb, na, b->limb, nb);
    if (n == 0) {
        cfx_big_mul(out, a, b);
        return;
    }
    out->n = n;
}

void cfx_big_mul_eq_fft(cfx_big_t* b, const cfx_big_t* m) {
    cfx_big_t result;
    cfx_big_init(&result);
    cfx_big_mul_fft(&result, b, m);
    cfx_big_swap(&result, b);
    cfx_big_free(&result);
}

void cfx_big_mul_eq_csa(cfx_big_t* b, const cfx_big_t* m) {
    if (cfx_big_is_zero(b) || cfx_big_is_zero(m)) {
        cfx_big_from_limb(b, 0);
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
    /* cfx_mul_scratch_alloc(scratch, nout); */
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

/*  out = a * b */
void cfx_big_mul(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    if (cfx_big_is_zero(a) || cfx_big_is_zero(b)) {
        cfx_big_from_limb(out, 0);
        return;
    }

    size_t na = a->n;
    size_t nb = b->n;

    cfx_big_reserve(out, na + nb);
    out->n = na + nb;
    memset(out->limb, 0, out->n * sizeof(cfx_limb_t));

    for (size_t i = 0; i < nb; ++i) {
        cfx_acc_t carry = 0;
        cfx_acc_t bi = (cfx_acc_t)b->limb[i];

        size_t k = i;
        for (size_t j = 0; j < na; ++j, ++k) {
            cfx_acc_t s = (cfx_acc_t)out->limb[k]
                        + bi * (cfx_acc_t)a->limb[j]
                        + carry;
            out->limb[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
        }

        while (carry) {
            cfx_acc_t s = (cfx_acc_t)out->limb[k] + carry;
            out->limb[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
            ++k;
        }
    }

    cfx_big_trim(out);
}

/* In-place multiplication: b *= m */
void cfx_big_mul_eq(cfx_big_t* b, const cfx_big_t* m) {

    if (cfx_big_is_zero(b) || cfx_big_is_zero(m)) {
        cfx_big_from_limb(b, 0);
        return;
    }
    /* not sure this is worth it:*/
    /* if (cfx_big_eq(b, m)) {
        cfx_big_sq_eq(b);
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
        cfx_acc_t carry = 0;
        cfx_acc_t mi = (cfx_acc_t)m->limb[i];

        size_t k = i; /* index into tmp */
        for (size_t j = 0; j < nb; ++j, ++k) {
            cfx_acc_t s = (cfx_acc_t)tmp.limb[k]
                        + mi * (cfx_acc_t)b->limb[j]
                        + carry;
            tmp.limb[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
        }

        /* propagate any remaining carry */
        while (carry) {
            cfx_acc_t s = (cfx_acc_t)tmp.limb[k] + carry;
            tmp.limb[k] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
            ++k; /* safe: nb+nm is enough for two nb/nm-digit numbers */
        }
    }

    cfx_big_trim(&tmp);
    cfx_big_swap(&tmp, b);
    cfx_big_free(&tmp);
}

/* out = a + b */
void cfx_big_add(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    /* handle aliasing: if out aliases a or b, work on a temp */
    if (out == a) {
        cfx_big_add_eq(out, b);
        return;
    }
    if (out == b) {
        cfx_big_add_eq(out, a);
        return;
    }

    /* out is distinct from a and b */
    cfx_big_assign(out, a);
    cfx_big_add_eq(out, b);
}

/* b += a */
void cfx_big_add_eq(cfx_big_t* b, const cfx_big_t* a) {
    cfx_limb_t carry = 0;
    size_t i = 0;

    while (i < a->n || carry) {
        cfx_big_reserve(b, i + 1);

        cfx_acc_t bi = (i < b->n) ? b->limb[i] : 0;
        cfx_acc_t ai = (i < a->n) ? a->limb[i] : 0;
        cfx_acc_t s  = bi + ai + carry;

        if (i >= b->n) b->n = i + 1;      /* we’re extending b */
        b->limb[i] = (cfx_limb_t)s;
        carry      = (cfx_limb_t)(s >> CFX_LIMB_BITS);
        ++i;
    }

    /* If a->n > b->n and carry ended at 0, we may still need to bump b->n */
    if (i > b->n) b->n = i;
}


void cfx_big_add_sm_eq(cfx_big_t* b, cfx_limb_t n) {
    if (n == 0) return;
    if (b->n == 0) {
        cfx_big_from_limb(b, n);
        return;
    }

    size_t i = 0;
    while (i < b->n) {
        cfx_acc_t s = (cfx_acc_t)b->limb[i] + n;
        b->limb[i] = (cfx_limb_t)s;
        n = (cfx_limb_t)(s >> CFX_LIMB_BITS);
        if (n == 0) {
            return;
        }
        ++i;
    }

    /* carry left; append one limb */
    cfx_big_reserve(b, b->n + 1);
    b->limb[b->n++] = n;      /* n < 2^64 here */
}

/* out = a - b (assumes a >= b) */
void cfx_big_sub(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b) {
    /* handle aliasing */
    if (out == a) {
        cfx_big_sub_eq(out, b);
        return;
    }
    if (out == b) {
        /* out = a - out, need temp */
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_assign(&tmp, a);
        cfx_big_sub_eq(&tmp, b);
        cfx_big_swap(out, &tmp);
        cfx_big_free(&tmp);
        return;
    }

    /* out is distinct from a and b */
    cfx_big_assign(out, a);
    cfx_big_sub_eq(out, b);
}

void cfx_big_sub_eq(cfx_big_t* a, const cfx_big_t* b) {
    /* assumes a >= b; subtract b from a */
    cfx_acc_t borrow = 0;
    size_t i = 0, n = b->n;
    for (; i < n; ++i) {
        cfx_acc_t ai = a->limb[i];
        cfx_acc_t s  = (cfx_acc_t)b->limb[i] + borrow;
        cfx_acc_t r  = ai - s;
        a->limb[i] = (cfx_limb_t)r;
        borrow = (ai < s);
    }
    while (borrow && i < a->n) {
        cfx_limb_t ai = a->limb[i];
        cfx_limb_t r  = ai - 1u;
        a->limb[i++] = r;
        borrow = (ai == 0u);
    }
    cfx_big_trim(a);
}

void cfx_big_sub_sm_eq(cfx_big_t* b, cfx_limb_t n) {
    if (n == 0 || b->n == 0) return;

    /* Subtract n from limb[0], propagate borrow if needed */
    cfx_limb_t borrow = (b->limb[0] < n) ? 1 : 0;
    b->limb[0] -= n;

    /* Propagate borrow through remaining limbs */
    size_t i = 1;
    while (borrow && i < b->n) {
        cfx_limb_t old = b->limb[i];
        b->limb[i] -= 1;
        borrow = (old == 0);
        ++i;
    }

    /* If borrow remains, we had underflow (b < n). This is a caller bug.
     * The result is now garbage (wrapped), but we don't crash. */
    assert(!borrow && "cfx_big_sub_sm_eq: underflow (b < n)");

    cfx_big_trim(b);
}


void cfx_big_mul_sm_eq(cfx_big_t* b, cfx_limb_t m) {
    if (m == 1) return;
    if (m == 0 || b->n == 0) {
        cfx_big_init(b);
        cfx_big_from_limb(b, 0);
        return;
    }
    _mul_sm_fast(b, m);
}

void cfx_big_from_limbs(cfx_big_t* b, const cfx_limb_t* limbs, size_t n) {
    cfx_big_free(b);
    cfx_big_reserve(b, n);
    memcpy(b->limb, limbs, n * sizeof(cfx_limb_t));
    b->n = n;
    cfx_big_trim(b);
}

/* Materialize factorization into cfx_big_t */
void cfx_big_from_fac(cfx_big_t* b, const cfx_fac_t *f) {
    cfx_big_from_limb(b, 1);
    for (size_t i = 0; i < f->len; i++){
        cfx_big_expmul_prime(b, f->data[i].p, f->data[i].e);
    }
}

/* ******************************************************************************* */
/* bucket-carry strategy for making bigs out of fac : */
#define CFX_FAC_BUCKETS 128 /* plenty: log2(bitlen(n!)) */

static void bucket_init(cfx_big_t buckets[CFX_FAC_BUCKETS], uint8_t used[CFX_FAC_BUCKETS]) {
    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        cfx_big_init(&buckets[i]);
        used[i] = 0;
    }
}

static void bucket_free(cfx_big_t buckets[CFX_FAC_BUCKETS], uint8_t used[CFX_FAC_BUCKETS]) {
    (void)used;
    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        cfx_big_free(&buckets[i]);
    }
}

/* Inserts x into buckets by multiplying upward while occupied. x is consumed (moved from). */
static void bucket_insert(cfx_big_t buckets[CFX_FAC_BUCKETS], uint8_t used[CFX_FAC_BUCKETS], cfx_big_t* x) {
    size_t level = 0;

    if (cfx_big_is_one(x)) {
        cfx_big_free(x);
        cfx_big_init(x);
        return;
    }

    while (level < CFX_FAC_BUCKETS) {
        if (!used[level]) {
            cfx_big_move(&buckets[level], x);  /* buckets[level] takes ownership */
            used[level] = 1;
            return;
        }

        /* x *= buckets[level]; clear bucket; carry upward */
        cfx_big_mul_auto(x, &buckets[level]);    /* x = x * buckets[level] */
        cfx_big_free(&buckets[level]);
        cfx_big_init(&buckets[level]);
        used[level] = 0;
        level++;
    }

    /* ---------> todo: increase CFX_FAC_BUCKETS. */
}

void cfx_big_from_fac_fast(cfx_big_t* out, const cfx_fac_t* f) {
    cfx_big_t buckets[CFX_FAC_BUCKETS];
    uint8_t   used[CFX_FAC_BUCKETS];

    bucket_init(buckets, used);


    for (size_t i = 0; i < f->len; i++) {
        const cfx_limb_t p = f->data[i].p;
        const cfx_limb_t e = f->data[i].e;
        if (!e) continue;

        cfx_big_t x;
        cfx_big_init(&x);
        cfx_big_from_limb(&x, 1);
        cfx_big_expmul_prime(&x, p, e);  /* x = p^e */

        bucket_insert(buckets, used, &x);

        cfx_big_free(&x);   /* x was moved into a bucket or freed; ensure no leak if insert didn't consume */
    }


    cfx_big_from_limb(out, 1);
    for (size_t level = 0; level < CFX_FAC_BUCKETS; level++) {
        if (used[level]) {
            cfx_big_mul_auto(out, &buckets[level]); /* out *= buckets[level] */
        }
    }

    bucket_free(buckets, used);
}

/* x > 0 */
static size_t floor_log2_size_t(size_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return (size_t)(8u * sizeof(size_t) - 1u - (unsigned)__builtin_clzl(x));
#else
    size_t r = 0;
    while (x >>= 1) r++;
    return r;
#endif
}

static void bucket_clear(cfx_big_t* b) { cfx_big_free(b); cfx_big_init(b); }


void cfx_big_from_fac_faster(cfx_big_t* out, const cfx_fac_t* f) {
    cfx_big_t buckets[CFX_FAC_BUCKETS];
    uint8_t used[CFX_FAC_BUCKETS];

    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        cfx_big_init(&buckets[i]);
        used[i] = 0;
    }

    for (size_t i = 0; i < f->len; i++) {
        cfx_limb_t p = f->data[i].p;
        cfx_limb_t e = f->data[i].e;
        if (!e) continue;

        cfx_big_t x;
        cfx_big_init(&x);
        cfx_big_pow_sm(&x, p, e);

        if (cfx_big_is_one(&x)) {
            cfx_big_free(&x);
            continue;
        }

        size_t bl = cfx_big_bitlen(&x);
        size_t level = (bl > 0) ? floor_log2_size_t(bl) : 0;
        if (level >= CFX_FAC_BUCKETS) level = CFX_FAC_BUCKETS - 1;

        for (;;) {
            if (!used[level]) {
                cfx_big_move(&buckets[level], &x);
                used[level] = 1;
                break;
            }
            /* x *= buckets[level], clear bucket, carry up */
            cfx_big_mul_auto(&x, &buckets[level]);
            bucket_clear(&buckets[level]);
            used[level] = 0;

            if (level + 1 < CFX_FAC_BUCKETS) {
                ++level;
            } else {
                /* extremely unlikely: fold into last bucket */
                if (!used[level]) {
                    cfx_big_move(&buckets[level], &x);
                    used[level] = 1;
                } else {
                    cfx_big_mul_auto(&buckets[level], &x);
                }
                break;
            }
        }

        cfx_big_free(&x);
    }

    cfx_big_from_limb(out, 1);
    for (size_t i = 0; i < CFX_FAC_BUCKETS; i++) {
        if (used[i]) cfx_big_mul_auto(out, &buckets[i]);
        cfx_big_free(&buckets[i]);
    }
}


/* Helper: check if a big integer fits in 64 bits */
static int big_fits_in_64(const cfx_big_t* b) {
#if CFX_LIMB_BITS == 64
    return b->n <= 1;
#elif CFX_LIMB_BITS == 32
    return b->n <= 2;
#endif
}

/* Helper: convert big integer to uint64 (assumes it fits) */
static uint64_t big_to_u64(const cfx_big_t* b) {
    if (b->n == 0) return 0;
#if CFX_LIMB_BITS == 64
    return b->limb[0];
#elif CFX_LIMB_BITS == 32
    uint64_t val = b->limb[0];
    if (b->n >= 2) val |= ((uint64_t)b->limb[1] << 32);
    return val;
#endif
}

/*
 * Factorize a big integer into prime factors.
 * Strategy:
 *   1. Trial division by small primes from cfx_primes[]
 *   2. If remainder fits in 64 bits, delegate to cfx_fac_from_u64
 *   3. Otherwise, use Pollard-Rho to find factors recursively
 *   4. If we encounter a prime > 64 bits, return incomplete
 *
 * Returns:
 *   0  - complete factorization
 *   1  - incomplete (remainder holds unfactored prime > 64 bits)
 *  -1  - error
 */
int cfx_big_to_fac(cfx_fac_t* f, const cfx_big_t* b, cfx_big_t* remainder) {
    cfx_fac_init(f);

    if (remainder) {
        cfx_big_from_limb(remainder, 1);
    }

    if (cfx_big_is_zero(b)) return 0;
    if (b->n == 1 && b->limb[0] == 1) return 0;

    /* work on a copy */
    cfx_big_t work;
    cfx_big_init(&work);
    cfx_big_copy(&work, b);

    /* Trial division by small primes */
    for (size_t i = 0; i < cfx_primes_len && !cfx_big_is_zero(&work); ++i) {
        uint32_t p = cfx_primes[i];

        /* If work < p, we're done with trial division */
        if (work.n == 1 && work.limb[0] < p) break;

        /* Count exponent */
        uint32_t e = 0;
        while (cfx_big_mod_sm(&work, p) == 0) {
            cfx_big_div_sm_eq(&work, p);
            cfx_big_trim(&work);
            e++;
        }
        if (e > 0) {
            cfx_fac_push(f, (uint64_t)p, e);
        }

        if (cfx_big_is_one(&work)) break;
    }

    /* If work == 1, factorization is complete */
    if (cfx_big_is_zero(&work) || cfx_big_is_one(&work)) {
        cfx_big_free(&work);
        return 0;
    }

    /* use a stack for iterative factorization */
    cfx_big_t stack[64];  /* Enough for any reasonable factorization */
    size_t stack_top = 0;
    stack_top = 1;
    cfx_big_init(&stack[0]);
    cfx_big_swap(&stack[0], &work);

    int result = 0;  /* 0 = complete */

    while (stack_top > 0) {
        cfx_big_t* cur = &stack[--stack_top];

        if (cfx_big_is_one(cur) || cfx_big_is_zero(cur)) {
            cfx_big_free(cur);
            continue;
        }

        /* If it fits in 64 bits, use fast 64-bit factorization */
        if (big_fits_in_64(cur)) {
            uint64_t val = big_to_u64(cur);
            cfx_fac_t fac64;
            if (cfx_fac_from_u64(&fac64, val) == 0) {
                for (size_t i = 0; i < fac64.len; ++i) {
                    cfx_fac_push(f, fac64.data[i].p, fac64.data[i].e);
                }
            }
            cfx_fac_free(&fac64);
            cfx_big_free(cur);
            continue;
        }

        /* Check if it's prime */
        if (cfx_big_is_prime(cur)) {
            /* Prime > 64 bits - can't store in fac_t */
            if (remainder) {
                /* Multiply remainder by this prime (remainder *= cur) */
                cfx_big_mul_eq(remainder, cur);
            }
            result = 1;  /* incomplete */
            cfx_big_free(cur);
            continue;
        }

        /* Use Pollard-Rho to find a factor */
        cfx_big_t factor;
        cfx_big_init(&factor);
        cfx_big_pollard_rho(&factor, cur);

        /* If Pollard-Rho failed (returned n), try a different approach */
        if (cfx_big_cmp(&factor, cur) == 0) {
            /* Last resort: the number might be a prime we couldn't detect */
            if (remainder) {
                /* Multiply remainder by this composite */
                cfx_big_mul_eq(remainder, cur);
            }
            result = 1;  /* incomplete */
            cfx_big_free(&factor);
            cfx_big_free(cur);
            continue;
        }

        /* Push both factors onto the stack */
        cfx_big_t quotient;
        cfx_big_init(&quotient);
        cfx_big_copy(&quotient, cur);
        cfx_big_divrem_eq(&quotient, &factor, NULL);
        cfx_big_trim(&quotient);

        /* Free cur BEFORE pushing, since cur points into the stack array
         * and we're about to reuse that slot. This avoids double-free. */
        cfx_big_free(cur);

        /* Push factor and quotient */
        if (stack_top < 62) {
            cfx_big_init(&stack[stack_top]);
            cfx_big_swap(&stack[stack_top], &factor);
            stack_top++;

            cfx_big_init(&stack[stack_top]);
            cfx_big_swap(&stack[stack_top], &quotient);
            stack_top++;
        }

        cfx_big_free(&factor);
        cfx_big_free(&quotient);
    }

    cfx_big_free(&work);

    /* Sort and coalesce the factorization */
    /* (cfx_fac_push maintains sorted order, but we might have duplicates) */
    if (f->len > 1) {
        /* Coalesce duplicate primes */
        cfx_fac_t coalesced;
        cfx_fac_init(&coalesced);

        /* Sort by prime (bubble sort - factors are usually few) */
        for (size_t i = 0; i < f->len; ++i) {
            for (size_t j = i + 1; j < f->len; ++j) {
                if (f->data[i].p > f->data[j].p) {
                    cfx_pf_t tmp = f->data[i];
                    f->data[i] = f->data[j];
                    f->data[j] = tmp;
                }
            }
        }

        /* Coalesce adjacent equal primes */
        uint64_t cur_p = f->data[0].p;
        uint32_t cur_e = f->data[0].e;
        for (size_t i = 1; i < f->len; ++i) {
            if (f->data[i].p == cur_p) {
                cur_e += f->data[i].e;
            } else {
                cfx_fac_push(&coalesced, cur_p, cur_e);
                cur_p = f->data[i].p;
                cur_e = f->data[i].e;
            }
        }
        cfx_fac_push(&coalesced, cur_p, cur_e);

        /* Swap */
        cfx_fac_free(f);
        *f = coalesced;
    }

    return result;
}

cfx_limb_t cfx_big_div_sm_eq(cfx_big_t* b, cfx_limb_t d) {
    cfx_acc_t rem = 0;
    for (size_t i = b->n; i--;) {
        cfx_acc_t cur = ((cfx_acc_t)rem << CFX_LIMB_BITS) | b->limb[i];
        cfx_limb_t q = (cfx_limb_t)(cur / d);
        rem = (cfx_limb_t)(cur % d);
        b->limb[i] = q;
    }
    return rem;
}

/* Divides x (base 2^64) by uint32_t d, returns remainder. */
uint32_t cfx_big_div_sm_u32_eq(cfx_big_t* b, uint32_t d) {
    cfx_limb_t rem = 0;
    for (size_t i = b->n; i--;) {
        cfx_acc_t cur = ((cfx_acc_t)rem << CFX_LIMB_BITS) | b->limb[i];
        cfx_limb_t q = (cfx_limb_t)(cur / d);   /* 128/32 -> 128/64 promoted; OK */
        rem = (cfx_limb_t)(cur % d);
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
cfx_limb_t cfx_big_mod_sm(const cfx_big_t* b, cfx_limb_t m) {
    if (b->n == 0) return 0;
    cfx_acc_t acc = 0;
    for (size_t i = b->n; i--;) {
        /* acc << 64 is acc * 2^64, which is (a_n * B); limb[i] is a_(n-1)*/
        acc = ((acc << CFX_LIMB_BITS) + b->limb[i]) % m;
    }
    return (cfx_limb_t)acc;
}

int cfx_fac_from_big(cfx_fac_t* fac, const cfx_big_t* in) {
    (void)fac;
    (void)in;
    /* TODO: implement big integer factorization into cfx_fac_t
     * See cfx_big_to_fac() for a working implementation that uses
     * trial division + Pollard-Rho. This function would need helper
     * functions like cfx_big_strip_small, cfx_big_fits_u64, etc. */
    return -1;
}

static inline size_t hex_digits_limb(cfx_limb_t v) {
    if (!v) return 1;
#if defined(__GNUC__) || defined(__clang__)
    #if CFX_LIMB_BITS == 64
        unsigned lead = (unsigned)__builtin_clzll(v);
    #else
        unsigned lead = (unsigned)__builtin_clz(v);
    #endif
    unsigned bits = CFX_LIMB_BITS - lead;
    return (bits + 3u) / 4u; /* ceil(bits/4) */
#else
    size_t d = 0;
    while (v) { v >>= 4; ++d; }
    return d;
#endif
}
static size_t bits_in_limb(cfx_limb_t x) {
    if (!x) return 0;
#if defined(__GNUC__) || defined(__clang__)
    #if CFX_LIMB_BITS == 64
        return CFX_LIMB_BITS - (size_t)__builtin_clzll(x);
    #else
        return CFX_LIMB_BITS - (size_t)__builtin_clz(x);
    #endif
#else
    size_t n = 0;
    while (x) { ++n; x >>= 1; }
    return n;
#endif
}

char* cfx_big_b64_alloc(const cfx_big_t* src, size_t* sz_out) {
    if (!src || src->n == 0) {
        char* s = (char*)malloc(2);
        if (!s) return NULL;
        s[0] = '0';
        s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    size_t bytecnt = 0;
    cfx_big_to_bytes_be(NULL, &bytecnt, src);
    uint8_t* bytes = (uint8_t*)malloc(bytecnt);
    cfx_big_to_bytes_be(bytes, &bytecnt, src);
    size_t charcnt = 0;
    cfx_base64_encode(NULL, &charcnt, bytes, bytecnt);
    char* s = (char*)malloc(charcnt);
    cfx_base64_encode(s, &charcnt, bytes, bytecnt);
    free(bytes);
    return s;
}

char* cfx_big_bin_alloc(const cfx_big_t* src, size_t* sz_out) {
    if (!src || src->n == 0) {
        char* s = (char*)malloc(2);
        if (!s) return NULL;
        s[0] = '0';
        s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    const size_t ms_idx  = src->n - 1;              /* most-significant limb index */
    const cfx_limb_t msval  = src->limb[ms_idx];
    const size_t ms_bits  = bits_in_limb(msval);      /* 1..CFX_LIMB_BITS */
    const size_t total_len = ms_bits + (size_t)CFX_LIMB_BITS * ms_idx; /* total bits as characters */

    char* s = (char*)malloc(total_len + 1);
    if (!s) return NULL;

    size_t pos = 0;
    const char bch[2] = {'0', '1'};

    for (size_t b = ms_bits; b-- > 0; ) {
        s[pos++] = bch[(size_t)((msval >> b) & 1u)];
    }


    for (size_t i = ms_idx; i-- > 0; ) {
        cfx_limb_t limb = src->limb[i];
        for (size_t b = CFX_LIMB_BITS; b-- > 0; ) {
            s[pos++] = bch[(size_t)((limb >> b) & 1u)];
        }
    }

    s[total_len] = '\0';
    if (sz_out) *sz_out = total_len;
    return s;
}

char* cfx_big_hex_alloc(const cfx_big_t* src, size_t* sz_out) {
    /* Treat empty/zero as "0" */
    if (!src || src->n == 0) {
        char* s = (char*)malloc(2);
        if (!s) return NULL;
        s[0] = '0'; s[1] = '\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    /* trim leading zeros */
    size_t ms = src->n;

    const cfx_limb_t ms_val = src->limb[ms - 1];
    const size_t ms_digits = hex_digits_limb(ms_val);
    const size_t hex_per_limb = CFX_LIMB_BITS / 4;  /* 8 for 32-bit, 16 for 64-bit */
    const size_t total_len = ms_digits + (ms - 1) * hex_per_limb;

    char* s = (char*)malloc(total_len + 1); /* +1 for NUL */
    if (!s) return NULL;

    char* p = s;
    size_t rem = total_len + 1;

    /* Most-significant limb without leading zeros */
    int written = snprintf(p, rem, CFX_PRIxLIMB, (cfx_limb_t)ms_val);
    assert(written > 0 && (size_t)written == ms_digits);
    p   += written;
    rem -= (size_t)written;

    /* Remaining limbs, zero-padded to hex_per_limb hex chars each */
    #ifdef CFX_DEBUG
    size_t k = 0;
    size_t cnt = 0;
    #endif
    for (size_t i = ms-1; i--;) {
        written = snprintf(p, rem, "" CFX_PRI0xLIMB, (cfx_limb_t)src->limb[i]);
        assert(written == (int)hex_per_limb);
        p   += written;
        #ifdef CFX_DEBUG
        cnt += written;
        #endif
        rem -= (size_t)written;

        #ifdef CFX_DEBUG
        if (!(cnt % 7)) {  /* some number.. */
            const char spinner[] = "|/-\\";
            ++k;
            CFX_PRINT_DBG("%zu hex digits done... %zu/%zu limbs remain... %c        \r",
                cnt, rem, total_len, spinner[k % 4]);
            fflush(stdout);
        }
        #endif
    }
    CFX_PRINT_DBG("\n");
    /* `snprintf` already wrote the final '\0' on the last call */
    if (sz_out) *sz_out = total_len;
    return s;
}

/* Convert cfx_big_t to decimal string */
char* cfx_big_dec_alloc(const cfx_big_t* src, size_t *sz_out) {
    if (src->n == 0) {
        char* s = (char*)malloc(2);
        s[0]='0';
        s[1]='\0';
        if (sz_out) *sz_out = 1;
        return s;
    }

    cfx_big_t tmp = *src;
    tmp.cap = tmp.n;
    tmp.limb = (cfx_limb_t*)malloc(tmp.n * sizeof(cfx_limb_t));
    memcpy(tmp.limb, src->limb, tmp.n * sizeof(cfx_limb_t));

    enum { CHUNK_BASE = 1000000000u, CHUNK_DIGS = 9 };
    size_t maxdig = src->n * 20; /* log10(2^64) == 19.2659... */
    uint32_t *chunks = (uint32_t*)malloc(maxdig);
    size_t k = 0;

    /* printf("[cfx_big_dec_alloc] building base %u representation... max digits: %zu\n",CHUNK_BASE, maxdig); */
    #ifdef CFX_DEBUG
    int cnt = 0;
    const size_t n0 = tmp.n;
    #endif
    while (tmp.n) {

        #ifdef CFX_DEBUG
        if (!(tmp.n % 100)) {
            const char spinner[] = "|/-\\";
            CFX_PRINT_DBG("%zu decimal digits done... %zu/%zu limbs remain... %c        \r",
                k*9, tmp.n, n0, spinner[cnt++ % 4]);
            fflush(stdout);
        }
        #endif
        chunks[k++] = cfx_big_div_sm_u32_eq(&tmp, CHUNK_BASE);
    }
    /* printf("\n"); */
    cfx_big_free(&tmp);

    /* build string */
    /* first chunk has no zero-padding, others are 9-digit padded */
    size_t len = 0;
    {
        uint32_t first = chunks[k-1];
        char buf[16];
        int l = snprintf(buf, sizeof buf, "%" PRIu32, first);
        len += l + (k-1) * CHUNK_DIGS;
    }
    char* s = (char*)malloc(len+1);
    char* p = s;

    /* write first */
    p += snprintf(p, len+1, "%" PRIu32, chunks[k-1]);
    /* write rest padded (indices k-2 -> 0) */
    for (size_t i = k-1; i--;) {
        p += snprintf(p, len+1, "%09" PRIu32, chunks[i]);
    }
    s[len] = '\0';
    if (sz_out) *sz_out = len;
    free(chunks);
    return s;
}

size_t cfx_big_snprint_dec(const cfx_big_t* b, char* out, size_t outlen) {
    size_t len;
    char* s = cfx_big_dec_alloc(b, &len);
    if (!s) {
        if (out && outlen > 0) out[0] = '\0';
        return 0;
    }
    if (out && outlen > 0) {
        size_t copy = (len < outlen) ? len : outlen - 1;
        memcpy(out, s, copy);
        out[copy] = '\0';
    }
    free(s);
    return len;
}

size_t cfx_big_snprint_hex(const cfx_big_t* b, char* out, size_t outlen) {
    size_t len;
    char* s = cfx_big_hex_alloc(b, &len);
    if (!s) {
        if (out && outlen > 0) out[0] = '\0';
        return 0;
    }
    if (out && outlen > 0) {
        size_t copy = (len < outlen) ? len : outlen - 1;
        memcpy(out, s, copy);
        out[copy] = '\0';
    }
    free(s);
    return len;
}

size_t cfx_big_snprint_bin(const cfx_big_t* b, char* out, size_t outlen) {
    size_t len;
    char* s = cfx_big_bin_alloc(b, &len);
    if (!s) {
        if (out && outlen > 0) out[0] = '\0';
        return 0;
    }
    if (out && outlen > 0) {
        size_t copy = (len < outlen) ? len : outlen - 1;
        memcpy(out, s, copy);
        out[copy] = '\0';
    }
    free(s);
    return len;
}

/* Scan a single numeric literal token at in[0..in_len).
   Accepts optional prefixes: 0x/0X, 0b/0B, 0o/0O.
   No internal whitespace, no sign.
   On success: parses into out, sets *consumed to total chars consumed, returns 0.
   On failure (not a number at start): sets *consumed=0, returns -1. */
int cfx_big_scan_num_n(cfx_big_t* out, const uint8_t* in, size_t in_len, size_t* consumed) {
    if (consumed) *consumed = 0;
    if (!out || (!in && in_len)) return -1;
    if (in_len == 0) return -1;

    size_t prefix = 0;
    int ret = 0;
    size_t digits = 0;

    /* base prefix detection */
    if (in_len >= 2 && in[0] == (uint8_t)'0') {
        uint8_t p = in[1];
        if (p == (uint8_t)'x' || p == (uint8_t)'X' ||
            p == (uint8_t)'b' || p == (uint8_t)'B' ||
            p == (uint8_t)'o' || p == (uint8_t)'O') {
            prefix = 2;
        }
    } else if (in_len >= 4 && strncmp((const char*)in, "b64:", 4) == 0) {
        prefix = 4;
    }

    // printf("SCAN_NUM prefix=0%c, calling bin on '%c' '%c' '%c'...\n",
    //     (char)in[1],
    //     (char)in[prefix + 0],
    //     (char)in[prefix + 1],
    //     (char)in[prefix + 2]);

    if (prefix) {
        if (prefix == 2) {
            uint8_t p = in[1];
            if (p == (uint8_t)'x' || p == (uint8_t)'X') {
                ret = cfx_big_scan_hex_n(out, in + prefix, in_len - prefix, &digits);
            } else if (p == (uint8_t)'b' || p == (uint8_t)'B') {
                ret = cfx_big_scan_bin_n(out, in + prefix, in_len - prefix, &digits);
            } else { /* o/O */
                ret = cfx_big_from_oct_n(out, in + prefix, in_len - prefix, &digits);
            }
        } else if (prefix == 4) {
            ret = cfx_big_scan_b64_n(out, in + prefix, in_len - prefix, &digits);
        }

        if (ret != 0) return -1;          /* no digits after prefix */
        if (consumed) *consumed = prefix + digits;
        return 0;
    }

    /* no prefix => decimal */
    ret = cfx_big_scan_dec_n(out, in, in_len, &digits);
    if (ret != 0) return -1;
    if (consumed) *consumed = digits;
    return 0;
}

int cfx_big_from_str(cfx_big_t* out, const char* s) {
    if (!out || !s) return -1;

    size_t len = strlen(s);
    if (len == 0) return -1;

    while (len > 0 && isspace(s[0])) {
        ++s;
        --len;
    }

    size_t consumed = 0;  /* digits consumed */
    size_t prefix_len = 0;
    int ret = 0;

    if (len >= 2 && (s[0]=='0' && (s[1]=='x' || s[1]=='X'))) {
        prefix_len = 2;
        ret = cfx_big_scan_hex_n(out, (const uint8_t*)(s+prefix_len), len-prefix_len, &consumed);
    } else if (len >=2 && (s[0]=='0' && (s[1]=='b' || s[1]=='B'))) {
        prefix_len = 2;
        ret = cfx_big_scan_bin_n(out, (const uint8_t*)(s+prefix_len), len-prefix_len, &consumed);
    } else if (len >= 4 && (s[0]=='b' || s[0]=='B') && s[1]=='6' && s[2]=='4' && s[3]==':') {
        prefix_len = 4;
        ret = cfx_big_scan_b64_n(out, (const uint8_t*)(s + prefix_len), len - prefix_len, &consumed);
    } else {
        prefix_len = 0;
        ret = cfx_big_scan_dec_n(out, (const uint8_t*)s, len, &consumed);
    }

    s += (prefix_len + consumed);

    while (*s && isspace((unsigned char)*s)) ++s;   /* allow trailing spaces only */
    if (ret != 0 || *s != '\0') return -1;          /* no digits or junk trailing */

    return 0;
}

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

int cfx_big_scan_hex_n(cfx_big_t* out, const uint8_t* in, size_t in_len, size_t* consumed) {
    cfx_big_assign_zero(out);
    if (consumed) *consumed = 0;
    size_t pos = 0;
    for (; pos < in_len; ++pos) {
        const uint8_t c = in[pos];
        int d = hex_table[c];
        if (d < 0) break;
        cfx_big_shl_bits(out, out, 4);
        cfx_big_add_sm_eq(out, (cfx_limb_t)d);
    }
    if (pos == 0) return -1;
    if (consumed) *consumed = pos;
    return 0;
}

int cfx_big_from_hex(cfx_big_t* out, const char* s) {
    while (*s && isspace((unsigned char)*s)) ++s;

    if (s[0]=='0' && (s[1]=='x' || s[1]=='X')) s += 2;
    size_t len = strlen(s);
    size_t digits = 0;
    int ret = cfx_big_scan_hex_n(out, (const uint8_t*)s, len, &digits);
    s += digits;
    while (isspace((unsigned char)*s)) ++s;
    if (ret != 0 || *s != '\0') return -1;
    return 0;
}

int cfx_big_scan_bin_n(cfx_big_t* out, const uint8_t* in, size_t in_len, size_t* consumed) {
    size_t ndig = 0;

    if (consumed) *consumed = 0;
    if (!out || (!in && in_len)) return -1;

    while (ndig < in_len) {
        uint8_t c = in[ndig];
        if (c != (uint8_t)'0' && c != (uint8_t)'1') break;
        ++ndig;
    }

    if (consumed) *consumed = ndig;
    if (ndig == 0) { cfx_big_assign_zero(out); return -1; }

    cfx_big_assign_zero(out);

    size_t limb_cnt = (ndig + CFX_LIMB_BITS - 1) / CFX_LIMB_BITS;
    cfx_big_reserve(out, limb_cnt);
    out->n = limb_cnt;

    for (size_t j = 0; j < ndig; ++j) {
        if (in[ndig - 1 - j] == (uint8_t)'1') {
            out->limb[j / CFX_LIMB_BITS] |= (cfx_limb_t)1 << (j % CFX_LIMB_BITS);
        }
    }

    while (out->n > 0 && out->limb[out->n - 1] == 0) --out->n;
    return 0;
}



int cfx_big_from_oct_n(cfx_big_t* out, const uint8_t* in, size_t in_len, size_t* consumed) {
    (void)out;
    (void)in;
    (void)in_len;
    (void)consumed;
    /* TODO: implement octal parsing */
    return -1;
}

int cfx_big_from_bin(cfx_big_t* b, const char* s) {
    size_t len = strlen(s);
    if ((len > 2) && (s[0] == '0' && ((s[1] == 'b') || (s[1] == 'B')))) {
        s += 2;
        len -= 2;
    }
    size_t digits = 0;
    int ret = cfx_big_scan_bin_n(b, (const uint8_t*)s, len, &digits);
    return ret;
}

int cfx_big_scan_dec_n(cfx_big_t* out, const uint8_t* in, size_t in_len, size_t* consumed) {
    cfx_big_assign_zero(out);
    if (consumed) *consumed = 0;
    size_t pos = 0;
    while (pos < in_len) {
        if (!isdigit((unsigned char)in[pos])) break;
        uint32_t digit = (uint32_t)(in[pos] - '0');
        cfx_big_mul_sm_eq(out, 10); /* out = out * 10 */
        cfx_big_add_sm_eq(out, digit); /* out = out + digit */
        ++pos;
    }
    if (pos == 0) return -1;
    if (consumed) *consumed = pos;
    return 0;
}

int cfx_big_from_dec(cfx_big_t* out, const char* s) {
    while (isspace((unsigned char)*s)) s++;
    size_t len = strlen(s);
    size_t digits = 0;
    int ret = cfx_big_scan_dec_n(out, (const uint8_t*)s, len, &digits);
    return ret;
}

static int is_b64_char(unsigned char c) {
    if (c == '+' || c == '/' || c == '=') return 1;
    if (c >= 'A' && c <= 'Z') return 1;
    if (c >= 'a' && c <= 'z') return 1;
    if (c >= '0' && c <= '9') return 1;
    /*if (isspace(c)) return 1;*/
    return 0;
}

/* Decodes 'in' as base 64 into 'out', up until the end of token: either '=' padding or first invalid character.
   Writes consumed characters to *consumed. */
int cfx_big_scan_b64_n(cfx_big_t* out, const uint8_t* in, size_t in_len, size_t* consumed) {
    if (!out || !in || !consumed) return -1;

    size_t n = 0;        /* bytes consumed from input */
    size_t nonws = 0;    /* count of non-ws chars consumed */
    size_t eq = 0;
    int seen_eq = 0;

    /* consume one base64 token */
    while (n < in_len) {
        unsigned char c = (unsigned char)in[n];

        if (isspace(c)) { n++; continue; } /* allow internal whitespace */

        if (c == '=') {
            seen_eq = 1;
            eq++;
            nonws++;
            n++;
            continue;
        }

        if (is_b64_char(c)) {
            if (seen_eq) break;      /* data after '=' ends token */
            nonws++;
            n++;
            continue;
        }

        break; /* non-b64, non-ws terminates token */
    }

    if (nonws == 0) return -1;
    if ((nonws % 4) != 0) return -1;
    if (eq > 2) return -1;

    /* now decode n bytes of base64: */
    size_t bin_len = 0;
    if (cfx_base64_decode(NULL, &bin_len, (const char *)in, n) != 0) return -1;

    uint8_t stack_buf[256];
    uint8_t *bin = stack_buf;

    if (bin_len > sizeof(stack_buf)) {
        bin = (uint8_t *)malloc(bin_len);
        if (!bin) return -1;
    }

    size_t got = bin_len;
    if (cfx_base64_decode(bin, &got, (const char *)in, n) != 0 || got != bin_len) {
        if (bin != stack_buf) free(bin);
        return -1;
    }

    if (cfx_big_from_bytes_be(out, bin, bin_len) != 0) {
        if (bin != stack_buf) free(bin);
        return -1;
    }

    if (bin != stack_buf) free(bin);

    *consumed = n;
    return 0;
}


/* ------------------------------------------------------------- */
/* NOTE: Use cfx_clz() from algo.h for limb-aware leading-zero count */

/* b <<= s (in-place, no temp allocation) */
void cfx_big_shl_bits_eq(cfx_big_t* b, unsigned s) {
    if (s == 0 || cfx_big_is_zero(b)) return;

    const unsigned W = CFX_LIMB_BITS;
    const unsigned limb_shift = s / W;
    const unsigned bit_shift = s % W;
    const size_t old_n = b->n;
    const size_t new_n = old_n + limb_shift + (bit_shift ? 1 : 0);

    cfx_big_reserve(b, new_n);

    if (bit_shift == 0) {
        /* Pure limb shift: move right-to-left to avoid overwriting */
        for (size_t i = old_n; i-- > 0; )
            b->limb[i + limb_shift] = b->limb[i];
        b->n = old_n + limb_shift;
    } else {
        /* combined limb + bit shift: work right-to-left */
        const unsigned r = W - bit_shift;
        cfx_limb_t carry = 0;
        for (size_t i = old_n; i-- > 0; ) {
            cfx_limb_t cur = b->limb[i];
            b->limb[i + limb_shift + 1] = carry | (cur >> r);
            carry = cur << bit_shift;
        }
        b->limb[limb_shift] = carry;
        b->n = new_n;
    }

    
    for (size_t i = 0; i < limb_shift; ++i) {
        b->limb[i] = 0;
    }

    cfx_big_trim(b);
}

/* b >>= s (in-place, no temp allocation) */
void cfx_big_shr_bits_eq(cfx_big_t* b, unsigned s) {
    if (s == 0 || cfx_big_is_zero(b)) return;

    const unsigned W = CFX_LIMB_BITS;
    const unsigned limb_shift = s / W;
    const unsigned bit_shift = s % W;

    if (b->n <= limb_shift) {
        cfx_big_assign_zero(b);
        return;
    }

    const size_t new_n = b->n - limb_shift;

    if (bit_shift == 0) {
        /* Pure limb shift: move left-to-right */
        for (size_t i = 0; i < new_n; ++i)
            b->limb[i] = b->limb[i + limb_shift];
    } else {
        /* Combined limb + bit shift: work left-to-right */
        const unsigned r = W - bit_shift;
        for (size_t i = 0; i < new_n; ++i) {
            cfx_limb_t lo = b->limb[i + limb_shift];
            cfx_limb_t hi = (i + limb_shift + 1 < b->n) ? b->limb[i + limb_shift + 1] : 0;
            b->limb[i] = (lo >> bit_shift) | (hi << r);
        }
    }

    b->n = new_n;
    cfx_big_trim(b);
}


/* out = a << s */
void cfx_big_shl_bits(cfx_big_t* out, const cfx_big_t* a, unsigned s) {
    if (cfx_big_is_zero(a) || s == 0) {
        cfx_big_copy(out, a);
        return;
    }

    cfx_big_t tmp;
    cfx_big_init(&tmp);

    if (a == out) {
        cfx_big_copy(&tmp, a);
        a = &tmp;          /* read from tmp */
        /* write into out */
    }

    const unsigned W = CFX_LIMB_BITS;
    const unsigned limb_shift = s / W;
    const unsigned bit_shift  = s % W;
    const unsigned r = bit_shift ? (W - bit_shift) : 0;

    size_t new_n = a->n + limb_shift + (bit_shift ? 1 : 0);
    cfx_big_reserve(out, new_n);

    /* (optional but good hygiene) clear destination region you’ll use */
    for (size_t i = 0; i < new_n; ++i) out->limb[i] = 0;

    if (bit_shift == 0) {
        for (size_t i = 0; i < a->n; ++i) {
            out->limb[i + limb_shift] = a->limb[i];
        }
        out->n = a->n + limb_shift;
    } else {
        cfx_limb_t carry = 0;
        for (size_t i = 0; i < a->n; ++i) {
            cfx_limb_t lo = a->limb[i];
            out->limb[i + limb_shift] = (lo << bit_shift) | carry;
            carry = (r ? (lo >> r) : 0);
        }
        out->limb[limb_shift + a->n] = carry;
        out->n = limb_shift + a->n + (carry ? 1 : 0);
    }

    cfx_big_trim(out);
    cfx_big_free(&tmp);
}



/* out = a >> s (0..63)  */
/* out = a >> s  (bitwise), base b = 2^64 */
void cfx_big_shr_bits(cfx_big_t* out, const cfx_big_t* a, unsigned s) {
    if (cfx_big_is_zero(a) || s == 0) {
        cfx_big_copy(out, a);
        return;
    }

    const unsigned W = CFX_LIMB_BITS;
    const unsigned limb_shift = s / W;
    const unsigned bit_shift  = s % W;

    if (a->n <= limb_shift) {
        cfx_big_from_limb(out, 0);
        return;
    }

    cfx_big_t tmp;
    cfx_big_init(&tmp);

    /* If aliasing, copy source into tmp and read from tmp. */
    if (a == out) {
        cfx_big_copy(&tmp, a);
        a = &tmp;
    }

    const size_t new_n = a->n - limb_shift;
    cfx_big_reserve(out, new_n);
    out->n = new_n;

    if (bit_shift == 0) {
        for (size_t i = 0; i < new_n; ++i) {
            out->limb[i] = a->limb[i + limb_shift];
        }
    } else {
        const unsigned r = W - bit_shift; /* 1..W-1 */
        for (size_t i = 0; i < new_n; ++i) {
            const cfx_limb_t lo = a->limb[i + limb_shift];
            const cfx_limb_t hi = (i + limb_shift + 1 < a->n) ? a->limb[i + limb_shift + 1] : 0;
            out->limb[i] = (lo >> bit_shift) | (hi << r);
        }
    }

    cfx_big_trim(out);
    cfx_big_free(&tmp);
}

/* Core: q = u / v; r = u % v; any of q or r may be NULL. Returns 0, or -1 if v==0. */
int cfx_big_divrem(cfx_big_t* q, cfx_big_t* r,
                   const cfx_big_t* u, const cfx_big_t* v) {
    #if CFX_LIMB_BITS==64
    const cfx_limb_t B_hi_bit = 1ull << 63;
    #elif CFX_LIMB_BITS==32
    const cfx_limb_t B_hi_bit = 1ul << 31;
    #endif

    /* ---- Trivial cases */
    if (cfx_big_is_zero(v)) return -1;

    if (cfx_big_is_zero(u)) {
        if (q) cfx_big_from_limb(q, 0);
        if (r) cfx_big_from_limb(r, 0);
        return 0;
    }

    if (cfx_big_cmp(u, v) < 0) {
        if (q) cfx_big_from_limb(q, 0);
        if (r) cfx_big_copy(r, u);
        return 0;
    }

    /* ---- single-limb fast path */
    if (v->n == 1) {
        /* PRINT_DBG("[FAST] single-limb divisor\n"); */
        cfx_limb_t div = v->limb[0];
        cfx_acc_t rem = 0;
        if (q) { cfx_big_reserve(q, u->n); q->n = u->n; }
        for (size_t i = u->n; i--;) {
            cfx_acc_t cur = (rem << CFX_LIMB_BITS) | u->limb[i];
            cfx_limb_t qi = (cfx_limb_t)(cur / div);
            rem = cur % div;
            if (q) q->limb[i] = qi;
            /* PRINT_DBG("  [FAST] i=%zd  cur=(" CFX_PRI0xLIMB "|" CFX_PRI0xLIMB ")  qi=" CFX_PRI0xLIMB "  rem=" CFX_PRI0xLIMB "\n", */
            /*        i, (cfx_limb_t)(cur >> CFX_LIMB_BITS), (cfx_limb_t)cur, */
            /*        qi, (cfx_limb_t)rem); */
        }
        if (q) cfx_big_trim(q);
        if (r) cfx_big_from_limb(r, (cfx_limb_t)rem);
        return 0;
    }

    /* ---- Algorithm D (Knuth) */
    /* D1. Normalize so MSB of V's top limb is set. */
    unsigned s = cfx_clz(v->limb[v->n - 1]);
    cfx_big_t V, U;
    cfx_big_init(&V);
    cfx_big_init(&U);
    cfx_big_shl_bits(&V, v, s);    /* V = v << s */
    cfx_big_shl_bits(&U, u, s);    /* U = u << s */

    PRINT_DBG("---- Normalize: s = %u (shift left by s bits)\n", s);
#ifdef CFX_PRINT_DEBUG
    PRINT_BIG("u (orig)     ", u);
    PRINT_BIG("v (orig)     ", v);
    PRINT_BIG("U = u << s   ", &U);
    PRINT_BIG("V = v << s   ", &V);
#endif
    /* Ensure the invariant V[n-1] >= B/2 */
    if (!(V.limb[V.n - 1] & B_hi_bit)) {
        PRINT_DBG("  [WARN] V not normalized: V[n-1] MSB not set!\n");
    }

    /* Append one extra zero limb to U (Knuth D1 convenience) */
    cfx_big_reserve(&U, U.n + 1);
    U.limb[U.n] = 0;
    U.n += 1;

    size_t n = V.n;          /* divisor length after normalization */
    size_t m_plus_n = U.n;   /* dividend length incl the extra top limb */
    size_t m = (m_plus_n >= n) ? (m_plus_n - n) : 0;  /* number of quotient digits (qlen) */

    /* PRINT_DBG("---- Sizes: n(divisor limbs)=%zu, U.n=%zu, qlen=%zu\n", n, U.n, m); */

    cfx_big_t Q;
    cfx_big_init(&Q);
    cfx_big_reserve(&Q, m);
    Q.n = m;
    for (size_t i = 0; i < m; ++i) Q.limb[i] = 0;

    const cfx_limb_t Vn1 = V.limb[n - 1];
    const cfx_limb_t Vn2 = V.limb[n - 2]; /* safe since n>=2 here */

    /* D2..D7 main loop: j = m-1 down to 0 */
    for (size_t j = m; j--;) {
        /* D3: compute qhat, rhat from top two limbs of U and top limb of V */
        cfx_limb_t Ujm  = U.limb[j + n];       /* U_{j+n} */
        cfx_limb_t Ujm1 = U.limb[j + n - 1];   /* U_{j+n-1} */

        cfx_limb_t qhat;
        cfx_acc_t rhat;  /* Use wider type to detect overflow */
        int skip_d3_refine = 0;

        if (Ujm == Vn1) {
            /* Ujm == Vn1 means the trial quotient would be >= b, so clamp to b-1.
             * According to Knuth, rhat = (Ujm*b + Ujm1) - qhat*Vn1
             *                         = (Vn1*b + Ujm1) - (b-1)*Vn1
             *                         = Ujm1 + Vn1
             * If this overflows (>= b), then D3 refinement is skipped entirely
             * since rhat >= b means the test qhat*Vn2 > rhat*b + Ujm2 is always false. */
            qhat = CFX_LIMB_MAX;               /* b-1 */
            rhat = (cfx_acc_t)Ujm1 + Vn1;      /* may be >= b */
            if (rhat >= ((cfx_acc_t)1 << CFX_LIMB_BITS)) {
                skip_d3_refine = 1;            /* rhat >= b, skip refinement */
            }
            /* PRINT_DBG("j=%zd  qhat CLAMP (Ujm==Vn1)  qhat=" CFX_PRI0xLIMB " rhat=0x%" PRIx64 " skip=%d\n", */
            /*        j, qhat, (uint64_t)rhat, skip_d3_refine); */
        } else {
            cfx_acc_t top = ((cfx_acc_t)Ujm << CFX_LIMB_BITS) | Ujm1;
            qhat = (cfx_limb_t)(top / Vn1);
            rhat = top % Vn1;
            /* PRINT_DBG("j=%zd  top=[" CFX_PRI0xLIMB "|" CFX_PRI0xLIMB "]/" CFX_PRI0xLIMB " -> qhat=" CFX_PRI0xLIMB " rhat=" CFX_PRI0xLIMB "\n", */
            /*        j, Ujm, Ujm1, */
            /*        Vn1, qhat, (cfx_limb_t)rhat); */
        }

        /* D3 refine: while qhat*V[n-2] > rhat*b + U[j+n-2], decrement qhat.
         * Only do this if rhat < b (i.e., fits in a limb). */
        if (!skip_d3_refine) {
            cfx_limb_t Ujm2 = U.limb[j + n - 2];
            cfx_acc_t lhs = (cfx_acc_t)qhat * Vn2;
            cfx_acc_t rhs = (rhat << CFX_LIMB_BITS) | Ujm2;

            /* First possible decrement */
            if (lhs > rhs) {
                /* PRINT_DBG("  D3: adjust qhat (lhs>rhs): qhat--, rhat+=Vn1\n"); */
                qhat--;
                rhat += Vn1; /* may now be >= b */
                /* Second possible decrement only if rhat < b */
                if (rhat < ((cfx_acc_t)1 << CFX_LIMB_BITS)) {
                    lhs = (cfx_acc_t)qhat * Vn2;
                    rhs = (rhat << CFX_LIMB_BITS) | Ujm2;
                    if (lhs > rhs) {
                        /* PRINT_DBG("  D3: second adjust qhat: qhat--, rhat+=Vn1\n"); */
                        qhat--;
                    }
                }
            }
        }

        /* D4: Multiply-subtract qhat * V from U window [j .. j+n] */
        cfx_acc_t carry = 0;
        for (size_t i = 0; i < n; ++i) {
            cfx_acc_t t = (cfx_acc_t)qhat * V.limb[i] + carry;
            cfx_limb_t t_lo = (cfx_limb_t)t;
            cfx_limb_t t_hi = (cfx_limb_t)(t >> CFX_LIMB_BITS);

            cfx_limb_t uji = U.limb[j + i];
            cfx_limb_t uji_new = uji - t_lo;
            cfx_limb_t borrow = (uji_new > uji);   /* 1 if underflow */

            U.limb[j + i] = uji_new;
            carry = (cfx_acc_t)t_hi + borrow;

            /* PRINT_DBG("  D4: j=%zd i=%zu  q*V=" CFX_PRI0xLIMB "|" CFX_PRI0xLIMB "  U=" CFX_PRI0xLIMB " -> " CFX_PRI0xLIMB "  carry=" CFX_PRI0xLIMB "(+%u)\n", */
            /*        j, i, */
            /*        t_hi, t_lo, */
            /*        uji, uji_new, */
            /*        (cfx_limb_t)carry, (unsigned)((cfx_limb_t)(carry >> CFX_LIMB_BITS))); */
        }

        /* D6: subtract final carry from U[j+n] */
        cfx_acc_t ujmw = (cfx_acc_t)U.limb[j + n];
        cfx_acc_t diff = ujmw - carry;
        cfx_limb_t ujm_new = (cfx_limb_t)diff;
        cfx_limb_t borrow_out = (ujmw < carry);  /* true iff underflow */

        /* PRINT_DBG("  D6: j=%zd  U[j+n]=" CFX_PRI0xLIMB " - carry=" CFX_PRI0xLIMB "(+%u) -> " CFX_PRI0xLIMB "  borrow_out=%u\n", */
        /*        j, */
        /*        (cfx_limb_t)ujmw, */
        /*        (cfx_limb_t)carry, */
        /*        (unsigned)((cfx_limb_t)(carry >> CFX_LIMB_BITS)), */
        /*        ujm_new, */
        /*        (unsigned)borrow_out); */

        U.limb[j + n] = ujm_new;

        /* D7: if we subtracted too much, add V back and decrement qhat */
        if (borrow_out) {
            /* PRINT_DBG("  D7: add-back (qhat too large): qhat=" CFX_PRI0xLIMB " -> " CFX_PRI0xLIMB "\n", */
            /*        qhat, (qhat - 1)); */
            qhat--;

            cfx_limb_t c = 0;
            for (size_t i = 0; i < n; ++i) {
                cfx_acc_t ssum = (cfx_acc_t)U.limb[j + i] + V.limb[i] + c;
                U.limb[j + i] = (cfx_limb_t)ssum;
                c = (cfx_limb_t)(ssum >> CFX_LIMB_BITS);
                /* PRINT_DBG("    add-back: i=%zu  U+=V+c -> U=" CFX_PRI0xLIMB "  c=" CFX_PRI0xLIMB "\n", */
                /*        i, U.limb[j + i], c); */
            }
            U.limb[j + n] += c;
            /* PRINT_DBG("    add-back: U[j+n]+=c -> " CFX_PRI0xLIMB "\n", */
            /*        U.limb[j + n]); */
        }

        Q.limb[j] = qhat;
        /* PRINT_DBG("  => Q[%zd] = " CFX_PRI0xLIMB "\n", j, qhat); */
    }

    /* D8: Unnormalize remainder: r = (U[0..n-1] >> s) */
    if (r) {
        cfx_big_t Rn;
        cfx_big_init(&Rn);
        cfx_big_reserve(&Rn, n);
        Rn.n = n;
        for (size_t i = 0; i < n; ++i) Rn.limb[i] = U.limb[i];
#ifdef CFX_PRINT_DEBUG
        PRINT_BIG("R (norm pre>>s)", &Rn);
#endif
        cfx_big_shr_bits(r, &Rn, s);
        cfx_big_free(&Rn);
        cfx_big_trim(r);
#ifdef CFX_PRINT_DEBUG
        PRINT_BIG("r (unnormalized)", r);
#endif
    }

    /* Output quotient */
    cfx_big_trim(&Q);
    if (q) {
        cfx_big_swap(&Q, q);
    }
    cfx_big_free(&Q);

    /* Cleanup */
    cfx_big_free(&U);
    cfx_big_free(&V);

#ifdef CFX_PRINT_DEBUG
    if (q) PRINT_BIG("Q (final)", q);
#endif
    return 0;
}

int cfx_big_div(cfx_big_t* q, const cfx_big_t* u, const cfx_big_t* v) {
    return cfx_big_divrem(q, NULL, u, v);
}

int cfx_big_mod(cfx_big_t* r, const cfx_big_t* u, const cfx_big_t* v) {
    return cfx_big_divrem(NULL, r, u, v);
}


int cfx_big_mulmod(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, const cfx_big_t* m) {
    cfx_big_t tmp;
    cfx_big_init(&tmp);
    cfx_big_mul(&tmp, a, b);
    int ret = cfx_big_mod(out, &tmp, m);
    cfx_big_free(&tmp);
    return ret;
}

/* In-place: u := floor(u/v); optional remainder r. Alias-safe for any combination. */
int cfx_big_divrem_eq(cfx_big_t* u, const cfx_big_t* v, cfx_big_t* r /*nullable*/) {
    cfx_big_t qtmp, rtmp;
    cfx_big_init(&qtmp);
    cfx_big_init(&rtmp);

    int rc = cfx_big_divrem(&qtmp, r ? &rtmp : NULL, u, v);
    if (rc == 0) {
        cfx_big_swap(&qtmp, u);
        if (r) cfx_big_swap(&rtmp, r);
    }

    cfx_big_free(&qtmp);
    cfx_big_free(&rtmp);
    return rc;
}


/* 128-bit + wrap counter accumulator per column. */
/* Value represented = hi * 2^128 + lo */
typedef struct {
    cfx_acc_t lo;
    cfx_limb_t hi;
    cfx_limb_t pad; /* pad to 32 bytes to reduce 0 sharing during reductions */
} acc128p_t;

/* add 64-bit x into accumulator (as a 128-bit add) */
static inline void acc_add_u64(acc128p_t* a, cfx_limb_t x) {
    cfx_acc_t old = a->lo;
    a->lo = old + (cfx_acc_t)x;
    a->hi += (a->lo < old); /* wrap over 2^128 */
}

/* add 128-bit x into accumulator */
/* static inline void acc_add_u128(acc128p_t* a, cfx_acc_t x) { */
/*     cfx_acc_t old = a->lo; */
/*     a->lo = old + x; */
/*     a->hi += (a->lo < old); */
/* } */

/* dst[k] += src[k] for k in [0..n) */
static inline void acc_vec_add(acc128p_t* dst, const acc128p_t* src, size_t n) {
    for (size_t k = 0; k < n; ++k) {
        cfx_acc_t old = dst[k].lo;
        dst[k].lo = old + src[k].lo;
        cfx_limb_t carry128 = (dst[k].lo < old);
        /* add hi parts + carry128; staying in 64 bits is fine for "millions of limbs" */
        dst[k].hi += src[k].hi + carry128;
    }
}

#if CFX_HAS_PTHREAD
/* Worker arguments */
typedef struct {
    const cfx_limb_t* a; size_t na;
    const cfx_limb_t* b; size_t nb;
    size_t j_begin, j_end;          /* range of rows (limbs of b) to process */
    acc128p_t* local_acc;           /* per-thread accumulator array (length = ncols) */
    size_t ncols;                   /* ncols = na + nb */
} rc_worker_args_t;



static void* worker_rowblock(void* vp) {
    rc_worker_args_t* w = (rc_worker_args_t*)vp;
    const cfx_limb_t* a = w->a;
    const cfx_limb_t* b = w->b;
    const size_t na = w->na;
    const size_t ncols = w->ncols;

    /* zero local accumulator */
    memset(w->local_acc, 0, ncols * sizeof(acc128p_t));

    for (size_t j = w->j_begin; j < w->j_end; ++j) {
        cfx_limb_t bj = b[j];
        if (!bj) continue; /* skip zero rows fast */

        for (size_t i = 0; i < na; ++i) {
            cfx_acc_t p  = (cfx_acc_t)a[i] * (cfx_acc_t)bj;    /* 64x64->128 */
            cfx_limb_t lo = (cfx_limb_t)p;
            cfx_limb_t hi = (cfx_limb_t)(p >> CFX_LIMB_BITS);
            size_t k = i + j;                   /* column for lo */
            /* Bounds: k in [0 .. na-1 + j] <= na-1 + (nb-1) < na+nb = ncols */
            acc_add_u64(&w->local_acc[k],     lo);
            acc_add_u64(&w->local_acc[k + 1], hi);
        }
    }
    return NULL;
}

/* Expand acc[k] = hi*2^128 + lo into base-2^64 lanes T[] and do one global carry pass. */
/* out must have length >= ncols + 3. */
static void expand_and_carry(const acc128p_t* acc, size_t ncols, cfx_limb_t* out)
{
    /* 1) clear T (we reuse 'out' as T) */
    memset(out, 0, (ncols + 3) * sizeof(cfx_limb_t));

    /* 2) expand each column into up to 3 neighboring columns, locally handling tiny carries */
    for (size_t k = 0; k < ncols; ++k) {
        cfx_acc_t lo = acc[k].lo;
        cfx_limb_t lo0 = (cfx_limb_t)lo;
        cfx_limb_t lo1 = (cfx_limb_t)(lo >> CFX_LIMB_BITS);
        cfx_limb_t hi0 = acc[k].hi; /* each unit here is exactly 2^128 = 2 limbs */

        /* T[k] += lo0 */
        cfx_acc_t t = (cfx_acc_t)out[k] + lo0;
        out[k] = (cfx_limb_t)t;
        cfx_limb_t c0 = (cfx_limb_t)(t >> CFX_LIMB_BITS);

        /* T[k+1] += lo1 + c0 */
        t = (cfx_acc_t)out[k+1] + lo1 + c0;
        out[k+1] = (cfx_limb_t)t;
        cfx_limb_t c1 = (cfx_limb_t)(t >> CFX_LIMB_BITS);

        /* T[k+2] += hi0 + c1 */
        t = (cfx_acc_t)out[k+2] + hi0 + c1;
        out[k+2] = (cfx_limb_t)t;
        cfx_limb_t c2 = (cfx_limb_t)(t >> CFX_LIMB_BITS);

        /* Any leftover tiny carry bubbles one more limb; global pass will finish everything. */
        out[k+3] += c2;
    }

    /* 3) single left-to-right carry sweep */
    cfx_limb_t carry = 0;
    for (size_t k = 0; k < ncols + 3; ++k) {
        cfx_acc_t t = (cfx_acc_t)out[k] + carry;
        out[k]  = (cfx_limb_t)t;
        carry   = (cfx_limb_t)(t >> CFX_LIMB_BITS);
    }
    if (carry) {
        /* ensure the caller reserved at least ncols+4 if they want to keep this */
        out[ncols + 3] = carry;
    }
}
#endif
/* Top-level: b *= m using row-parallel accumulation + single carry pass. */
/* threads<=0 -> auto (online CPUs). threads is capped to nb. */
void cfx_big_mul_rows_pthreads(cfx_big_t* b, const cfx_big_t* m, int threads)
{
    #if CFX_HAS_PTHREAD
    const size_t na = b->n;
    const size_t nb = m->n;

    if (!na || !nb) { cfx_big_from_limb(b, 0); return; }
    if (nb == 1 && m->limb[0] == 1) { return; } /* b *= 1 */

    /* Copy multiplicand 'a' because we'll overwrite b. */
    cfx_limb_t* a_copy = (cfx_limb_t*)malloc(na * sizeof(cfx_limb_t));
    if (!a_copy) { /* handle OOM */ abort(); }
    memcpy(a_copy, b->limb, na * sizeof(cfx_limb_t));

    /* Plan threads */
    int hw = cfx_cpu_count();
    if (threads <= 0) threads = hw;
    if (threads > (int)nb) threads = (int)nb;
    if (threads < 1) threads = 1;

    const size_t ncols = na + nb;          /* columns 0..(na+nb-1) */
    const size_t acc_len = ncols;          /* acc arrays length */
    const size_t out_len = ncols + 4;      /* room for expansion + final carry */

    /* per-thread local accumulators */
    acc128p_t** locals = (acc128p_t**)malloc(threads * sizeof(acc128p_t*));
    rc_worker_args_t* args = (rc_worker_args_t*)malloc(threads * sizeof(rc_worker_args_t));
    pthread_t* tids = (pthread_t*)malloc(threads * sizeof(pthread_t));
    if (!locals || !args || !tids) { abort(); }

    for (int t = 0; t < threads; ++t) {
        /* 64B-align to reduce 0 sharing during reduction */
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

        int rc = pthread_create(&tids[t], NULL, worker_rowblock, &args[t]);
        if (rc != 0) { abort(); }
    }

    for (int t = 0; t < threads; ++t) {
        pthread_join(tids[t], NULL);
    }

    /* Reduce locals -> global accumulator */
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

    /* Expand to base-2^64 lanes and do one global carry pass */
    cfx_limb_t* out = (cfx_limb_t*)calloc(out_len, sizeof(cfx_limb_t));
    if (!out) { abort(); }
    expand_and_carry(acc, acc_len, out);

    /* Normalize: find actual limb length (trim leading zeros) */
    size_t rn = out_len;
    while (rn > 0 && out[rn - 1] == 0) { --rn; }
    if (rn == 0) { cfx_big_from_limb(b, 0); }
    else {
        cfx_big_reserve(b, rn);
        /* (Assumes cfx_big_reserve zeros new space; if not, it's fine we overwrite.) */
        memcpy(b->limb, out, rn * sizeof(cfx_limb_t));
        b->n = rn;
    }

    free(out);
    free(acc);
    for (int t = 0; t < threads; ++t) free(locals[t]);
    free(tids);
    free(locals);
    free(args);
    free(a_copy);
    #else
    /* Fallback: use regular multiplication when pthreads not available */
    (void)threads;
    cfx_big_mul_eq(b, m);
    #endif
}

/* threshold for switching to NTT multiplication (in limbs)
   - with native __uint128_t: NTT modular arithmetic is fast
   - without (MSVC): portable mod128 is slow, need much larger threshold */
#ifndef CFX_NTT_LIMB_THRESHOLD
#  if CFX_HAS_UINT128
#    define CFX_NTT_LIMB_THRESHOLD 7400   /* measured on Mac: crossover 7200-7400 limbs */
#  else
#    define CFX_NTT_LIMB_THRESHOLD 50000  /* measured on MSVC: crossover ~45k limbs */
#  endif
#endif

void cfx_big_mul_auto(cfx_big_t* b, const cfx_big_t* m) {
    const size_t na = b->n;
    const size_t nb = m->n;

    /* trivial cases */
    if (!na || !nb) { cfx_big_from_limb(b, 0); return; }
    if (nb == 1 && m->limb[0] == 1) return;           /* b *= 1 */
    if (na == 1 && b->limb[0] == 1) {                  /* 1 * m */
        cfx_big_reserve(b, nb);
        memcpy(b->limb, m->limb, nb * sizeof(cfx_limb_t));
        b->n = nb;
        return;
    }

    const size_t mn = (na < nb) ? na : nb;

    if (mn < 256) {
        /* small/skinny -> classic single-thread schoolbook */
        cfx_big_mul_eq(b, m);
        return;
    }

    if (mn >= CFX_NTT_LIMB_THRESHOLD) {
        /* very large -> NTT O(n log n) multiplication */
        cfx_big_mul_eq_fft(b, m);
        return;
    }

    /* medium size: threaded row/col/carry */
    /* make sure the *row* operand is the larger one to expose threads */
    if (nb >= na) {
        /* m already larger (or equal): rows = m */
        cfx_big_mul_rows_pthreads(b, m, -1);  /* -1 = auto threads*/
    } else {
        /* m is smaller: compute (m *= b) into a temp, then move into b */
        cfx_big_t tmp;
        cfx_big_init(&tmp);
        cfx_big_copy(&tmp, m);                /* tmp = m */
        cfx_big_mul_rows_pthreads(&tmp, b, -1);
        cfx_big_swap(b, &tmp);
        cfx_big_free(&tmp);
    }
}




/* Precomputed 10^k for k=0..18 */
static const cfx_limb_t POW10U64[CFX_LIMB_DIGITS_DEC + 1] = {
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
    #if (CFX_LIMB_BITS == 64)
    10000000000ULL,
    100000000000ULL,
    1000000000000ULL,
    10000000000000ULL,
    100000000000000ULL,
    1000000000000000ULL,
    10000000000000000ULL,
    100000000000000000ULL,
    1000000000000000000ULL,
    /*10000000000000000000ULL */
    #endif
};

static void flush_chunk(cfx_big_t* out, unsigned base, cfx_limb_t* chunk_val, unsigned* chunk_len) {
    if (*chunk_len == 0) return;

    if (base == 10) {
        cfx_big_mul_sm_eq(out, POW10U64[*chunk_len]);
        cfx_big_add_sm_eq(out, *chunk_val);
    } else if (base == 16) {
        unsigned bits = 4u * (*chunk_len);
        cfx_big_shl_bits_eq(out, bits);
        cfx_big_add_sm_eq(out, *chunk_val);
    } else { /* base 2 */
        unsigned bits = *chunk_len;
        cfx_big_shl_bits_eq(out, bits);
        cfx_big_add_sm_eq(out, *chunk_val);
    }

    *chunk_val = 0;
    *chunk_len = 0;
}

/* Return 0 on success, nonzero on parse error. */
int cfx_big_from_file(cfx_big_t* out, FILE* fp, int base) {
    cfx_big_from_limb(out, 0);

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
    cfx_limb_t chunk_val = 0;
    unsigned chunk_len = 0; /* digits in current chunk */

    unsigned dec_chunk_max = CFX_LIMB_DIGITS_DEC;
    unsigned hex_chunk_max = CFX_LIMB_DIGITS_HEX;

    unsigned char buf[64 * 1024];
    size_t nread;
    while ((nread = fread(buf, 1, sizeof(buf), fp)) > 0) {
        for (size_t i = 0; i < nread; ++i) {
            unsigned char c = buf[i];

            /* allow underscores, quotes, spaces, newlines, tabs as visual separators */
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
                    chunk_val = chunk_val * 10u + (cfx_limb_t)(c - '0');
                    chunk_len++;
                    continue;
                }
            } else if (detected_base == BASE_HEX) {
                int v = hex_table[c];

                if (v != -1) {
                    saw_digit = 1;
                    if (chunk_len == hex_chunk_max) flush_chunk(out, detected_base, &chunk_val, &chunk_len);
                    chunk_val = (chunk_val << 4) | (cfx_limb_t)v;
                    chunk_len++;
                    continue;
                }
            } else { /* BASE_BIN */
                if (c == '0' || c == '1') {
                    saw_digit = 1;
                    if (chunk_len == 60u) flush_chunk(out, detected_base, &chunk_val, &chunk_len); /* keep under 64 bits */
                    chunk_val = (chunk_val << 1) | (cfx_limb_t)(c - '0');
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
                /* errno = EINVAL; */
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
        /* errno = ERANGE; -> or implement signed handling */
        CFX_PRINT_ERR("negative number found");
        return -1;
    }

    return 0;
}

/* ==================== Montgomery ==================== */

/* if x >= n then x -= n */
static inline void cfx_big_cond_sub_n_(cfx_big_t* x, const cfx_big_t* n) {
    if (cfx_big_cmp(x, n) >= 0) cfx_big_sub_eq(x, n);
}

/* n0^{-1} mod 2^64 via Newton; n0 must be odd */
static inline cfx_limb_t inv64_mod2_64_(cfx_limb_t n0) {
    cfx_limb_t x = 1;
    for (int i = 0; i < 6; ++i) x *= (2 - x * n0);
    return x;
}
static inline cfx_limb_t mont_n0inv_(cfx_limb_t n0) {
    return (cfx_limb_t)(0 - inv64_mod2_64_(n0)); /* -n^{-1} mod 2^64 */
}

/* rr = R^2 mod n, with R=2^(64*k). Build rr = 2^(128*k) mod n by repeated doubling. */
static int compute_rr_(cfx_big_t* rr, const cfx_big_t* n, size_t k) {
    cfx_big_from_limb(rr, 1);
    for (size_t bit = 0; bit < 128ull * k; ++bit) {
        /* rr = (rr + rr) mod n (rr < n ⇒ rr+rr < 2n ⇒ one cond. subtract is enough) */
        cfx_big_reserve(rr, rr->n + 1);
        cfx_acc_t carry = 0;
        size_t i = 0;
        for (; i < rr->n; ++i) {
            cfx_acc_t s = (cfx_acc_t)rr->limb[i] * 2 + carry;
            rr->limb[i] = (cfx_limb_t)s;
            carry = s >> CFX_LIMB_BITS;
        }
        if (carry) rr->limb[rr->n++] = (cfx_limb_t)carry;
        cfx_big_cond_sub_n_(rr, n);
    }
    return 1;
}

int cfx_big_mont_ctx_init(cfx_big_mont_ctx_t* ctx, const cfx_big_t* n_in) {
    if (!ctx || !n_in || n_in->n == 0) return 0;
    if ((n_in->limb[0] & 1ull) == 0) return 0; /* modulus must be odd */

    cfx_big_init(&ctx->n);
    cfx_big_init(&ctx->rr);
    cfx_big_init(&ctx->R1);

    cfx_big_assign(&ctx->n, n_in);
    cfx_big_trim(&ctx->n);

    ctx->k = ctx->n.n;
    if (ctx->k == 0) { cfx_big_mont_ctx_free(ctx); return 0; }

    ctx->n0inv = mont_n0inv_(ctx->n.limb[0]);


    if (!compute_rr_(&ctx->rr, &ctx->n, ctx->k)) {
        cfx_big_mont_ctx_free(ctx);
        return 0;
    }
    /* R1 calculation */
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_from_limb(&one, 1);
    int ok = cfx_big_mont_to(&ctx->R1, &one, ctx); /* R1 = 1*RR/R mod n */
    cfx_big_free(&one);

    if (!ok) return 0;
    return 1;
}

void cfx_big_mont_ctx_free(cfx_big_mont_ctx_t* ctx) {
    if (!ctx) return;
    cfx_big_free(&ctx->n);
    cfx_big_free(&ctx->rr);
    memset(ctx, 0, sizeof(*ctx));
}

int cfx_big_mont_mul(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* b, const cfx_big_mont_ctx_t* ctx) {
    if (!out || !a || !b || !ctx) return 0;

    const cfx_big_t* n   = &ctx->n;
    const size_t     k   = ctx->k;       /* number of limbs in n */
    const cfx_limb_t   n0i = ctx->n0inv;   /* n0i = -n[0]^{-1} mod 2^64 */

    cfx_big_t T;
    cfx_big_init(&T);
    cfx_big_reserve(&T, k + 2);          /* we use indices [0..k+1] */
    memset(T.limb, 0, (k + 2) * sizeof(cfx_limb_t));
    T.n = k + 2;

    const cfx_limb_t* an = a->limb;
    const cfx_limb_t* bn = b->limb;
    const cfx_limb_t* nn = n->limb;

    for (size_t i = 0; i < k; ++i) {
        const cfx_limb_t bi = (i < b->n) ? bn[i] : 0;

        /* T += a * bi */
        cfx_acc_t carry = 0;
        for (size_t j = 0; j < k; ++j) {
            const cfx_limb_t aj = (j < a->n) ? an[j] : 0;
            cfx_acc_t sum = (cfx_acc_t)T.limb[j] + (cfx_acc_t)aj * bi + carry;
            T.limb[j] = (cfx_limb_t)sum;
            carry     = sum >> CFX_LIMB_BITS;
        }
        cfx_acc_t top = (cfx_acc_t)T.limb[k] + carry;
        T.limb[k]     = (cfx_limb_t)top;
        T.limb[k + 1] += (cfx_limb_t)(top >> CFX_LIMB_BITS);

        /* m = (T[0] * n0i) mod 2^64  (n0i == -n[0]^{-1} mod 2^64) */
        const cfx_limb_t m = T.limb[0] * n0i;

        /* T += m * n */
        cfx_acc_t carry2 = 0;
        for (size_t j = 0; j < k; ++j) {
            cfx_acc_t sum = (cfx_acc_t)T.limb[j] + (cfx_acc_t)m * nn[j] + carry2;
            T.limb[j] = (cfx_limb_t)sum;
            carry2    = sum >> CFX_LIMB_BITS;
        }
        cfx_acc_t top2 = (cfx_acc_t)T.limb[k] + carry2;
        T.limb[k]     = (cfx_limb_t)top2;
        T.limb[k + 1] += (cfx_limb_t)(top2 >> CFX_LIMB_BITS);

        /* T = (T + m*n) / R: drop limb 0 -> shift 1..k+1 to 0..k */
        memmove(&T.limb[0], &T.limb[1], (k + 1) * sizeof(cfx_limb_t));  /* <<< FIXED */
        T.limb[k + 1] = 0;  /* scratch for next iter */
    }

    /* Final normalization */
    T.n = k + 1;                 /* keep possible top carry for trim/compare */
    cfx_big_trim(&T);
    if (cfx_big_cmp(&T, n) >= 0) {
        cfx_big_sub_eq(&T, n);
    }

    cfx_big_move(out, &T);
    /* cfx_big_free(&T); */
    return 1;
}

/* out = a^e mod n (normal domain). n must be odd. */
int cfx_big_modexp_binary(cfx_big_t* out, const cfx_big_t* a, const cfx_big_t* e, const cfx_big_mont_ctx_t* ctx) {
    /* e == 0 -> 1 mod n */
    if (e->n == 0) {
        return cfx_big_mont_from(out, &ctx->R1, ctx); /* MontFrom(R1) == 1 */
    }

    cfx_big_t baseR, accR;
    cfx_big_init(&baseR);
    cfx_big_init(&accR);

    if (!cfx_big_mont_to(&baseR, a, ctx)) goto FAIL;   /* a in Montgomery */
    cfx_big_assign(&accR, &ctx->R1);                   /* acc = 1 (Montgomery) */

    /* msb index of e */
    size_t msb;
    {
        cfx_limb_t top = e->limb[e->n - 1];
        unsigned lz  = cfx_clz(top);
        msb = (e->n * CFX_LIMB_BITS) - 1u - lz;
    }

    /* LTR binary exponentiation */
    for (size_t i = msb + 1; i-- > 0; ) {
        /* square every round */
        if (!cfx_big_mont_mul(&accR, &accR, &accR, ctx)) goto FAIL;

        /* conditional multiply when bit is 1 */
        size_t limb = i / CFX_LIMB_BITS, bit = i % CFX_LIMB_BITS;
        if (limb < e->n && ((e->limb[limb] >> bit) & 1u)) {
            if (!cfx_big_mont_mul(&accR, &accR, &baseR, ctx)) goto FAIL;
        }
    }

    /* back to normal domain */
    if (!cfx_big_mont_from(out, &accR, ctx)) goto FAIL;
    cfx_big_free(&baseR);
    cfx_big_free(&accR);
    return 1;
FAIL:
    cfx_big_free(&baseR);
    cfx_big_free(&accR);
    return 0;
}


/* Convert to Montgomery: aR mod n = MontMul(a, R^2) */
int cfx_big_mont_to(cfx_big_t* out, const cfx_big_t* a, const cfx_big_mont_ctx_t* ctx) {
    /* fast path: if a fits and < n, skip the mod */
    if (a->n <= ctx->k) {
        if (cfx_big_cmp(a, &ctx->n) >= 0) {
            cfx_big_t t; cfx_big_init(&t);
            cfx_big_assign(&t, a);
            /* At most a couple subs when a.n == k */
            do { cfx_big_sub_eq(&t, &ctx->n); } while (cfx_big_cmp(&t, &ctx->n) >= 0);
            int ok = cfx_big_mont_mul(out, &t, &ctx->rr, ctx);
            cfx_big_free(&t);
            return ok;
        }
        return cfx_big_mont_mul(out, a, &ctx->rr, ctx);
    }

    /*  a is larger than k limbs, take remainder once */
    cfx_big_t rem;
    cfx_big_init(&rem);
    cfx_big_mod(&rem, a, &ctx->n);          /* Knuth D remainder (u mod n) */
    int ok = cfx_big_mont_mul(out, &rem, &ctx->rr, ctx);
    cfx_big_free(&rem);
    return ok;
}

/* Convert from Montgomery: aR * R^{-1} = a mod n = MontMul(aR, 1) */
int cfx_big_mont_from(cfx_big_t* out, const cfx_big_t* aR, const cfx_big_mont_ctx_t* ctx) {
    if (!out || !aR || !ctx) return 0;
    cfx_big_t one;
    cfx_big_init(&one);
    cfx_big_from_limb(&one, 1);
    int ok = cfx_big_mont_mul(out, aR, &one, ctx);
    cfx_big_free(&one);
    return ok;
}

/* ----------------- One-liners that hide the ctx ----------------- */
/* these are pretty useless except for testing */
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
    cfx_big_mont_ctx_t C;
    if (!cfx_big_mont_ctx_init(&C, n)) return 0;
    cfx_big_t aR, r;
    cfx_big_init(&aR);
    cfx_big_init(&r);
    int ok = cfx_big_mont_to(&aR, a, &C)
           && cfx_big_mont_sqr(&r, &aR, &C)
           && cfx_big_mont_from(out, &r, &C);
    cfx_big_free(&aR);
    cfx_big_free(&r);
    cfx_big_mont_ctx_free(&C);
    return ok;
}

int cfx_big_modexp(cfx_big_t* out, const cfx_big_t* base, const cfx_big_t* exp, const cfx_big_t* n) {
    cfx_big_mont_ctx_t C;
    if (!cfx_big_mont_ctx_init(&C, n)) return 0;

    cfx_big_t A, X;
    cfx_big_init(&A);
    cfx_big_init(&X);

    int ok = cfx_big_mont_to(&A, base, &C);
    if (ok) {
        cfx_big_t one;
        cfx_big_init(&one);
        cfx_big_from_limb(&one, 1);
        ok = cfx_big_mont_to(&X, &one, &C);
        cfx_big_free(&one);
    }
    if (ok) {
        for (size_t i = 0; i < exp->n; ++i) {
            cfx_limb_t w = exp->limb[i];
            for (int b = 0; b < 64; ++b) {
                if (w & 1u) {
                    ok = cfx_big_mont_mul(&X, &X, &A, &C);
                    if (!ok) break;
                }
                w >>= 1;
                ok = cfx_big_mont_mul(&A, &A, &A, &C);
                if (!ok) break; /* square */
            }
            if (!ok) break;
        }
        if (ok) ok = cfx_big_mont_from(out, &X, &C);
    }

    cfx_big_free(&A);
    cfx_big_free(&X);
    cfx_big_mont_ctx_free(&C);
    return ok;
}



double cfx_big_log(const cfx_big_t* b, double base) {
    size_t l = b->n;
    cfx_limb_t hi = b->limb[l-1];
    double ln_base = log(base);
    double ln_hi = log(hi);
    double ln_B = 64 * log(2.0);
    return ((l - 1) * ln_B + ln_hi) / ln_base;
}

int cfx_big_to_sci(const cfx_big_t* x, unsigned base,
                                int sig_digits, char* out, size_t outsz) {

    if (!x || !out || outsz == 0 || base < 2) return 0;

    /* zero? */
    if (x->n == 0) {
        snprintf(out, outsz, "0");
        return 1;
    }

    /* find top nonzero limb (defensive) */
    size_t k = x->n;
    while (k > 0 && x->limb[k-1] == 0) --k;
    if (k == 0) { snprintf(out, outsz, "0"); return 1; }

    cfx_limb_t hi = x->limb[k-1];

    const long double lnB   = (long double)CFX_LIMB_BITS * logl(2.0L);
    const long double lnb   = logl((long double)base);

    /* ln(n) ≈ (k-1)*ln(B) + ln(hi) */
    long double ln_n = ((long double)(k - 1)) * lnB + logl((long double)hi);

    /* e = floor( ln(n) / ln(b) ) */
    long double logb_n = ln_n / lnb;
    long long e = (long long) floorl(logb_n);

    /* mantissa m = b^( fractional part ) */
    long double frac = logb_n - (long double)e;
    long double m = expl(frac * lnb);   /* == powl(base, frac);  guarantees 1 <= m < b (up to rounding) */

    /* rounding guard: if m rounds to exactly b, bump e and renormalize */
    if (!(m < (long double)base)) { m /= (long double)base; ++e; }

    /* format: decimal mantissa * base^e */
    if (sig_digits < 1) sig_digits = 1;
    int after_decimal = sig_digits - 1;

    if (base == 10) {
        /* 3.236e8930 */
        snprintf(out, outsz, "%.*Lf" "e%lld", after_decimal, m, e);
    } else {
        /* 3.236 * 7^12345 */
        snprintf(out, outsz, "%.*Lf x %u^%lld", after_decimal, m, base, e);
    }
    return 1;
}
