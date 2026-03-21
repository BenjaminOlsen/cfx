/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "big_internal.h"

#include "cfx/algo.h"
#include "cfx/ntt.h"
#include "cfx/compat.h"

/* Platform-specific includes */
#if defined(_WIN32) || defined(_WIN64)
    #include <intrin.h>
#else
    #include <unistd.h>
    #ifdef _POSIX_THREADS
        #define CFX_HAS_PTHREAD 1
        #include <pthread.h>
    #endif
#endif


void cfx_big_sq_eq(cfx_big_t *b) {

    const size_t n = b->n;
    cfx_big_t ret;
    cfx_big_init(&ret);

    if (n == 0) {
        return;
    }

    cfx_big_reserve(&ret, 2*n);
    memset(ret.limb, 0, 2*n * sizeof(cfx_limb_t));
    ret.n = 2*n;

    /* 1) Cross terms: for i < j, add b[i]*b[j] into ret[i+j] (single) */
    for (size_t i = 0; i < n; ++i) {
        cfx_acc_t carry = 0;
        for (size_t j = i + 1; j < n; ++j) {
            cfx_acc_t p = (cfx_acc_t)b->limb[i] * b->limb[j];
            cfx_acc_t t = (cfx_acc_t)ret.limb[i + j]
                          + (cfx_limb_t)p
                          + (cfx_limb_t)carry;
            ret.limb[i + j] = (cfx_limb_t)t;
            carry = (carry >> CFX_LIMB_BITS) + (t >> CFX_LIMB_BITS) + (p >> CFX_LIMB_BITS);
        }

        size_t k = i + n;
        while (carry) {
            cfx_acc_t t = (cfx_acc_t)ret.limb[k] + (cfx_limb_t)carry;
            ret.limb[k] = (cfx_limb_t)t;
            carry = (carry >> CFX_LIMB_BITS) + (t >> CFX_LIMB_BITS);
            ++k;
        }
    }

    /* 2) Double cross terms: left-shift ret by 1 bit */
    cfx_limb_t shl_carry = 0;
    for (size_t i = 0; i < 2*n; ++i) {
        cfx_limb_t new_carry = ret.limb[i] >> (CFX_LIMB_BITS - 1);
        ret.limb[i] = (ret.limb[i] << 1) | shl_carry;
        shl_carry = new_carry;
    }

    /* 3) Diagonals: add b[i]^2 at ret[2*i] */
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

void cfx_big_mul_fft(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
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

void cfx_big_mul_eq_fft(cfx_big_t *b, const cfx_big_t *m) {
    cfx_big_t result;
    cfx_big_init(&result);
    cfx_big_mul_fft(&result, b, m);
    cfx_big_swap(&result, b);
    cfx_big_free(&result);
}

void cfx_big_mul_eq_csa(cfx_big_t *b, const cfx_big_t *m) {
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
void cfx_big_mul_csa_scratch(cfx_big_t *b, const cfx_big_t *m, cfx_mul_scratch_t *scratch) {
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
void cfx_big_mul(cfx_big_t *out, const cfx_big_t *a, const cfx_big_t *b) {
    /* Edge case: zero inputs */
    if (cfx_big_is_zero(a) || cfx_big_is_zero(b)) {
        cfx_big_from_limb(out, 0);
        return;
    }

    /* Delegate to backend implementation */
    cfx_big_mul_impl(out, a, b);
}

/* In-place multiplication: b *= m */
void cfx_big_mul_eq(cfx_big_t *b, const cfx_big_t *m) {
    /* Edge case: zero inputs */
    if (cfx_big_is_zero(b) || cfx_big_is_zero(m)) {
        cfx_big_from_limb(b, 0);
        return;
    }

    /* Delegate to backend (handles aliasing internally) */
    cfx_big_mul_impl(b, b, m);
}

/* 128-bit + wrap counter accumulator per column. */
/* Value represented = hi * 2^128 + lo */
typedef struct {
    cfx_acc_t lo;
    cfx_limb_t hi;
    cfx_limb_t pad; /* pad to 32 bytes to reduce 0 sharing during reductions */
} acc128p_t;

/* add 64-bit x into accumulator (as a 128-bit add) */
static inline void acc_add_u64(acc128p_t *a, cfx_limb_t x) {
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
static inline void acc_vec_add(acc128p_t *dst, const acc128p_t *src, size_t n) {
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
    const cfx_limb_t *a; size_t na;
    const cfx_limb_t *b; size_t nb;
    size_t j_begin, j_end;          /* range of rows (limbs of b) to process */
    acc128p_t *local_acc;           /* per-thread accumulator array (length = ncols) */
    size_t ncols;                   /* ncols = na + nb */
} rc_worker_args_t;



static void * worker_rowblock(void *vp) {
    rc_worker_args_t *w = (rc_worker_args_t *)vp;
    const cfx_limb_t *a = w->a;
    const cfx_limb_t *b = w->b;
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
static void expand_and_carry(const acc128p_t *acc, size_t ncols, cfx_limb_t *out){
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
void cfx_big_mul_rows_pthreads(cfx_big_t *b, const cfx_big_t *m, int threads){
    #if CFX_HAS_PTHREAD
    const size_t na = b->n;
    const size_t nb = m->n;

    if (!na || !nb) {
        cfx_big_from_limb(b, 0); return;
    }
    if (nb == 1 && m->limb[0] == 1) {
        return;
    }                                           /* b *= 1 */

    /* Copy multiplicand 'a' because we'll overwrite b. */
    cfx_limb_t *a_copy = (cfx_limb_t *)malloc(na * sizeof(cfx_limb_t));
    if (!a_copy) { /* handle OOM */
        abort();
    }
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
    acc128p_t **locals = (acc128p_t **)malloc(threads * sizeof(acc128p_t *));
    rc_worker_args_t *args = (rc_worker_args_t *)malloc(threads * sizeof(rc_worker_args_t));
    pthread_t *tids = (pthread_t *)malloc(threads * sizeof(pthread_t));
    if (!locals || !args || !tids) {
        abort();
    }

    for (int t = 0; t < threads; ++t) {
        /* 64B-align to reduce 0 sharing during reduction */
        void *p = NULL;
#if defined(_ISOC11_SOURCE)
        p = aligned_alloc(64, ((acc_len * sizeof(acc128p_t) + 63) / 64) * 64);
        if (!p) {
            abort();
        }
#else
        if (posix_memalign(&p, 64, acc_len * sizeof(acc128p_t)) != 0) {
            abort();
        }
#endif
        locals[t] = (acc128p_t *)p;

        size_t chunk = (nb + threads - 1) / threads;
        size_t j0 = (size_t)t * chunk;
        size_t j1 = j0 + chunk; if (j1 > nb) j1 = nb;

        args[t].a = a_copy; args[t].na = na;
        args[t].b = m->limb; args[t].nb = nb;
        args[t].j_begin = j0; args[t].j_end = j1;
        args[t].local_acc = locals[t];
        args[t].ncols = ncols;

        int rc = pthread_create(&tids[t], NULL, worker_rowblock, &args[t]);
        if (rc != 0) {
            abort();
        }
    }

    for (int t = 0; t < threads; ++t) {
        pthread_join(tids[t], NULL);
    }

    /* Reduce locals -> global accumulator */
    acc128p_t *acc = NULL;
#if defined(_ISOC11_SOURCE)
    acc = (acc128p_t *)aligned_alloc(64, ((acc_len * sizeof(acc128p_t) + 63) / 64) * 64);
    if (!acc) {
        abort();
    }
#else
    void *accp = NULL;
    if (posix_memalign(&accp, 64, acc_len * sizeof(acc128p_t)) != 0) {
        abort();
    }
    acc = (acc128p_t *)accp;
#endif
    memset(acc, 0, acc_len * sizeof(acc128p_t));

    for (int t = 0; t < threads; ++t) {
        acc_vec_add(acc, locals[t], acc_len);
    }

    /* Expand to base-2^64 lanes and do one global carry pass */
    cfx_limb_t *out = (cfx_limb_t *)calloc(out_len, sizeof(cfx_limb_t));
    if (!out) {
        abort();
    }
    expand_and_carry(acc, acc_len, out);

    /* Normalize: find actual limb length (trim leading zeros) */
    size_t rn = out_len;
    while (rn > 0 && out[rn - 1] == 0) {
        --rn;
    }
    if (rn == 0) {
        cfx_big_from_limb(b, 0);
    } else {
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

void cfx_big_mul_auto(cfx_big_t *b, const cfx_big_t *m) {
    const size_t na = b->n;
    const size_t nb = m->n;

    /* trivial cases */
    if (!na || !nb) {
        cfx_big_from_limb(b, 0); return;
    }
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
