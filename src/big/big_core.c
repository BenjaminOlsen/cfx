/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include <inttypes.h>


/* cfx_big_init and cfx_big_free are provided by the memory backend */
/* See: src/big/mem/dynamic/ or src/big/mem/static/ */

void cfx_big_clear(cfx_big_t *b) {
    b->n = 0;
}

int cfx_big_copy(cfx_big_t *dst, const cfx_big_t *src) {
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

void cfx_big_assign(cfx_big_t *dst, const cfx_big_t *src) {
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

void cfx_big_assign_sm(cfx_big_t *dst, const cfx_limb_t src) {
    if (src == 0) {
        dst->n = 0;
        if (dst->cap) memset(dst->limb, 0, dst->cap * sizeof(cfx_limb_t));
    }

    cfx_big_reserve(dst, 1); /* todo: error if != 0 */

    memset(dst->limb, 0, dst->cap * sizeof(cfx_limb_t));
    dst->n = 1;
    dst->limb[0] = src;
}

void cfx_big_assign_zero(cfx_big_t *b) {
    b->n = 0;
    if (b->cap == 0) cfx_big_reserve(b, 1);
    b->limb[0] = 0;
}

void cfx_big_assign_one(cfx_big_t *b) {
    if (b->cap == 0) cfx_big_reserve(b, 1);
    b->n = 1;
    b->limb[0] = 1;
}


void cfx_big_move(cfx_big_t *dst, cfx_big_t *src) {
    if (dst == src) return;
    cfx_big_free(dst);
    *dst = *src;
    cfx_big_init(src);
}

int cfx_big_is_zero(const cfx_big_t *b) {
    return (b->n == 0) || (b->n == 1 && b->limb[0] == 0);
}

int cfx_big_is_one(const cfx_big_t *b) {
    return (b->n == 1 && b->limb[0] == 1);
}

int cfx_big_is_even(const cfx_big_t *b) {
    /* zero is even (0 = 2*0), otherwise check least significant bit */
    return (b->n == 0) || !(b->limb[0] & 0x1);
}

int cfx_big_eq_u64(const cfx_big_t *b, cfx_limb_t n) {
    return (n == 0 && cfx_big_is_zero(b)) || (b->n == 1 && b->limb[0] == n);
}

int cfx_big_eq(const cfx_big_t *a, const cfx_big_t *b) {
    size_t max_n = (a->n > b->n) ? a->n : b->n;
    cfx_limb_t diff = 0;
    for (size_t i = 0; i < max_n; ++i) {
        cfx_limb_t ai = (i < a->n) ? a->limb[i] : 0;
        cfx_limb_t bi = (i < b->n) ? b->limb[i] : 0;
        diff |= (ai ^ bi);
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
int cfx_big_cmp(const cfx_big_t *a, const cfx_big_t *b) {
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
int cfx_big_cmp_sm(const cfx_big_t *a, cfx_limb_t b) {
    if (a->n == 0) /* a == 0 */ return b == 0 ? 0 : -1;
    if (a->n > 1) return 1;
    if (a->limb[0] != b) return (a->limb[0] < b) ? -1 : 1;
    return 0;
}

/* compare with uint64_t - handles 32-bit limb builds too */
int cfx_big_cmp_u64(const cfx_big_t *a, uint64_t b) {
#if CFX_LIMB_BITS == 64
    return cfx_big_cmp_sm(a, b);
#else
    /* 32-bit limbs: b might need 2 limbs */
    uint32_t lo = (uint32_t)b;
    uint32_t hi = (uint32_t)(b >> 32);
    if (hi == 0) {
        return cfx_big_cmp_sm(a, lo);
    }
    /* b needs 2 limbs */
    if (a->n == 0) return -1;
    if (a->n > 2) return 1;
    if (a->n == 1) return -1;  /* a has 1 limb, b needs 2 */
    /* a->n == 2 */
    if (a->limb[1] != hi) return (a->limb[1] < hi) ? -1 : 1;
    if (a->limb[0] != lo) return (a->limb[0] < lo) ? -1 : 1;
    return 0;
#endif
}

void cfx_big_swap(cfx_big_t *a, cfx_big_t *b) {
    if (a == b) return;
    cfx_big_t tmp = *a;
    *a = *b;
    *b = tmp;
}

/*
 * Constant-time (maybe) conditional swap: if condition != 0, swap a and b.
 */
void cfx_big_cswap(cfx_big_t *a, cfx_big_t *b, int condition) {
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

size_t cfx_big_bitlen(const cfx_big_t *b) {
    /* Returns 0 for zero; otherwise assumes no leading zero limbs */
    if (b->n == 0) return 0;
    const size_t limb_bits = 8u * sizeof(cfx_limb_t);
    cfx_limb_t v = b->limb[b->n - 1];
    size_t lz = cfx_clz(v);
    return b->n * limb_bits - lz;
}

/* cfx_big_reserve is provided by the memory backend */


int cfx_big_from_u64(cfx_big_t *b, uint64_t v) {
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
int cfx_big_from_limb(cfx_big_t *b, cfx_limb_t v) {
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
int cfx_big_to_bytes_be(uint8_t *out, size_t *outlen, const cfx_big_t *b) {

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
int cfx_big_from_bytes_be(cfx_big_t *out, const uint8_t *be, size_t len) {
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
        size_t take = (src_end > off) ? ((src_end - off) < lb ? (src_end - off) : lb) : 0;

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

int cfx_big_endswith_u64(const cfx_big_t *x, uint64_t value) {
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

void cfx_big_from_limbs(cfx_big_t *b, const cfx_limb_t *limbs, size_t n) {
    cfx_big_free(b);
    cfx_big_reserve(b, n);
    memcpy(b->limb, limbs, n * sizeof(cfx_limb_t));
    b->n = n;
    cfx_big_trim(b);
}
