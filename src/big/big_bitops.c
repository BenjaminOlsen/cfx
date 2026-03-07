/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

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
static void big_bitop(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b, bitop_t op) {
    /* XOR with self is zero */
    if (op == BITOP_XOR && a == b) {
        cfx_big_assign_zero(out);
        return;
    }

    const size_t an = a->n, bn = b->n;
    const size_t min_n = (an < bn) ? an : bn;
    const size_t max_n = (an > bn) ? an : bn;
    const cfx_big_t *lg = (an >= bn) ? a : b;

    /* AND shrinks to min, OR/XOR grow to max */
    const size_t out_n = (op == BITOP_AND) ? min_n : max_n;

    /* Handle aliasing: copy inputs if out overlaps */
    cfx_big_t tmp_a, tmp_b;
    int alias_a = (out == a), alias_b = (out == b);
    if (alias_a) {
        cfx_big_init(&tmp_a); cfx_big_copy(&tmp_a, a); a = &tmp_a;
    }
    if (alias_b) {
        cfx_big_init(&tmp_b); cfx_big_copy(&tmp_b, b); b = &tmp_b;
    }

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

void cfx_big_and(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
    big_bitop(out, a, b, BITOP_AND);
}

void cfx_big_or(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
    big_bitop(out, a, b, BITOP_OR);
}

void cfx_big_xor(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
    big_bitop(out, a, b, BITOP_XOR);
}



/* In-place variants */
void cfx_big_and_eq(cfx_big_t *a, const cfx_big_t *b) {
    big_bitop(a, a, b, BITOP_AND);
}

void cfx_big_or_eq(cfx_big_t *a, const cfx_big_t *b) {
    big_bitop(a, a, b, BITOP_OR);
}

void cfx_big_xor_eq(cfx_big_t *a, const cfx_big_t *b) {
    big_bitop(a, a, b, BITOP_XOR);
}

void cfx_big_rotl_w(cfx_big_t *out, const cfx_big_t *b, unsigned r, unsigned w) {
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

void cfx_big_rotl(cfx_big_t *out, const cfx_big_t *b, unsigned r) {
    cfx_big_rotl_w(out, b, r, cfx_big_bitlen(b));
}

void cfx_big_rotr_w(cfx_big_t *out, const cfx_big_t *b, unsigned r, unsigned w) {
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

void cfx_big_rotr(cfx_big_t *out, const cfx_big_t *b, unsigned r) {
    cfx_big_rotr_w(out, b, r, cfx_big_bitlen(b));
}

void cfx_big_mask_bits(cfx_big_t *a, unsigned nbits) {
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

/* Single-bit operations */

int cfx_big_bit_is_set(const cfx_big_t *x, size_t bit) {
    if (!x || x->n == 0) return 0;

    const size_t limb_idx = bit / CFX_LIMB_BITS;
    const size_t bit_idx  = bit % CFX_LIMB_BITS;

    if (limb_idx >= x->n) return 0;

    return (x->limb[limb_idx] >> bit_idx) & 1;
}

void cfx_big_bit_set(cfx_big_t *x, size_t bit) {
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

void cfx_big_bit_clear(cfx_big_t *x, size_t bit) {
    const size_t limb_idx = bit / CFX_LIMB_BITS;
    const size_t bit_idx  = bit % CFX_LIMB_BITS;

    if (limb_idx >= x->n) return; /* already zero */

    x->limb[limb_idx] &= ~((cfx_limb_t)1 << bit_idx);
    cfx_big_trim(x);
}

void cfx_big_bit_flip(cfx_big_t *x, size_t bit) {
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

size_t cfx_big_popcount(const cfx_big_t *x) {
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

/* ------------------------------------------------------------- */
/* NOTE: Use cfx_clz() from algo.h for limb-aware leading-zero count */

/* b <<= s (in-place, no temp allocation) */
void cfx_big_shl_bits_eq(cfx_big_t *b, unsigned s) {
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
void cfx_big_shr_bits_eq(cfx_big_t *b, unsigned s) {
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
void cfx_big_shl_bits(cfx_big_t *out, const cfx_big_t *a, unsigned s) {
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

    /* (optional but good hygiene) clear destination region you'll use */
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
void cfx_big_shr_bits(cfx_big_t *out, const cfx_big_t *a, unsigned s) {
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
