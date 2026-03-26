/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */
/*
 * Poly1271: Polynomial MAC using the Mersenne prime 2^127 - 1
 *
 * Based on the reference implementation from github.com/poly1271/poly1271
 */

#include "cfx/poly1271.h"
#include "cfx/memory.h"
#include <string.h>

#define load64_le cfx_load64_le
#define store64_le cfx_store64_le

/* CONSTANTS */

#define M127_HI 0x7FFFFFFFFFFFFFFFULL
#define M127_LO 0xFFFFFFFFFFFFFFFFULL

/* 64x64 -> 128 BIT MULTIPLICATION */

#if defined(__SIZEOF_INT128__)
typedef unsigned __int128 uint128_t;
static inline void mul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi) {
    uint128_t r = (uint128_t)a * b;
    *lo = (uint64_t)r;
    *hi = (uint64_t)(r >> 64);
}
#elif defined(_MSC_VER) && defined(_M_X64)
#include <intrin.h>
static inline void mul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi) {
    *lo = _umul128(a, b, hi);
}
#else
static inline void mul64(uint64_t a, uint64_t b, uint64_t *lo, uint64_t *hi) {
    uint32_t a0 = (uint32_t)a, a1 = (uint32_t)(a >> 32);
    uint32_t b0 = (uint32_t)b, b1 = (uint32_t)(b >> 32);
    uint64_t p00 = (uint64_t)a0 * b0;
    uint64_t p01 = (uint64_t)a0 * b1;
    uint64_t p10 = (uint64_t)a1 * b0;
    uint64_t p11 = (uint64_t)a1 * b1;
    uint64_t mid = p01 + (p00 >> 32) + (uint32_t)p10;
    *lo = (p00 & 0xFFFFFFFF) | (mid << 32);
    *hi = p11 + (mid >> 32) + (p10 >> 32);
}
#endif

/* CONSTANT-TIME ADDITION WITH CARRY */

static inline uint64_t adc(uint64_t *a, uint64_t b) {
    uint64_t old = *a;
    *a += b;
    return (*a < old);
}

/* MODULAR REDUCTION */

static void mul_reduce(uint64_t out[3], const uint64_t a[3], const uint64_t r[2]) {
    uint64_t p0l, p0h, p1l, p1h, p2l, p2h, p3l, p3h, p4l, p4h, p5l, p5h;
    mul64(a[0], r[0], &p0l, &p0h);
    mul64(a[0], r[1], &p1l, &p1h);
    mul64(a[1], r[0], &p2l, &p2h);
    mul64(a[1], r[1], &p3l, &p3h);
    mul64(a[2], r[0], &p4l, &p4h);
    mul64(a[2], r[1], &p5l, &p5h);

    uint64_t r0 = p0l, c = 0;
    uint64_t r1 = p0h;
    c = adc(&r1, p1l);
    c += adc(&r1, p2l);
    uint64_t r2 = p1h + c;
    c = (r2 < p1h);
    c += adc(&r2, p2h);
    c += adc(&r2, p3l);
    c += adc(&r2, p4l);
    uint64_t r3 = p3h + c;
    c = (r3 < p3h);
    c += adc(&r3, p4h);
    c += adc(&r3, p5l);
    uint64_t r4 = p5h + c;

    uint64_t h0 = (r1 >> 63) | (r2 << 1);
    uint64_t h1 = (r2 >> 63) | (r3 << 1);
    uint64_t h2 = (r3 >> 63) | (r4 << 1);
    uint64_t h3 = r4 >> 63;

    out[0] = r0;
    out[1] = r1 & M127_HI;
    c = adc(&out[0], h0);
    c = adc(&out[1], h1 + c);
    out[2] = h2 + h3 + c;

    uint64_t fold = out[2] << 1;
    c = adc(&out[0], fold);
    out[1] += c;
    out[2] = 0;
}

static void mul_2x2(uint64_t out[3], const uint64_t a[2], const uint64_t r[2]) {
    uint64_t p0l, p0h, p1l, p1h, p2l, p2h, p3l, p3h;
    mul64(a[0], r[0], &p0l, &p0h);
    mul64(a[0], r[1], &p1l, &p1h);
    mul64(a[1], r[0], &p2l, &p2h);
    mul64(a[1], r[1], &p3l, &p3h);

    uint64_t r0 = p0l, c = 0;
    uint64_t r1 = p0h;
    c = adc(&r1, p1l);
    c += adc(&r1, p2l);
    uint64_t r2 = p1h + c;
    c = (r2 < p1h);
    c += adc(&r2, p2h);
    c += adc(&r2, p3l);
    uint64_t r3 = p3h + c;

    uint64_t h0 = (r1 >> 63) | (r2 << 1);
    uint64_t h1 = (r2 >> 63) | (r3 << 1);
    uint64_t h2 = r3 >> 63;
    out[0] = r0;
    out[1] = r1 & M127_HI;
    c = adc(&out[0], h0);
    c = adc(&out[1], h1 + c);
    out[2] = h2 + c;
}

static void reduce_full(uint64_t acc[3]) {
    uint64_t hi = (acc[1] >> 63) | (acc[2] << 1);
    uint64_t hh = acc[2] >> 63;
    uint64_t lo = acc[0];
    uint64_t mid = acc[1] & M127_HI;
    uint64_t c = adc(&lo, hi);
    mid += hh + c;

    uint64_t ov = mid >> 63;
    mid &= M127_HI;
    c = adc(&lo, ov);
    mid += c;
    ov = mid >> 63;
    mid &= M127_HI;
    lo += ov;

    acc[0] = lo;
    acc[1] = mid;
    acc[2] = 0;
}

static void mul_mod(uint64_t out[2], const uint64_t a[2], const uint64_t b[2]) {
    uint64_t t[3] = {a[0], a[1], 0};
    mul_reduce(t, t, b);
    reduce_full(t);
    out[0] = t[0];
    out[1] = t[1];
}

static void square_mod(uint64_t out[2], const uint64_t r[2]) {
    uint64_t p0l, p0h, p1l, p1h, p2l, p2h;
    mul64(r[0], r[0], &p0l, &p0h);
    mul64(r[0], r[1], &p1l, &p1h);
    mul64(r[1], r[1], &p2l, &p2h);

    uint64_t p1l2 = p1l << 1;
    uint64_t p1h2 = (p1h << 1) | (p1l >> 63);

    uint64_t r0 = p0l, c = 0;
    uint64_t r1 = p0h;
    c = adc(&r1, p1l2);
    uint64_t r2 = p1h2 + c;
    c = (r2 < p1h2);
    c += adc(&r2, p2l);
    uint64_t r3 = p2h + c;

    uint64_t h0 = (r1 >> 63) | (r2 << 1);
    uint64_t h1 = (r2 >> 63) | (r3 << 1);
    out[0] = r0;
    out[1] = r1 & M127_HI;
    c = adc(&out[0], h0);
    c = adc(&out[1], h1 + c);

    uint64_t ov = (out[1] >> 63) + c;
    out[1] &= M127_HI;
    c = adc(&out[0], ov);
    out[1] += c;
    ov = out[1] >> 63;
    out[1] &= M127_HI;
    out[0] += ov;
}

/* BLOCK PROCESSING */

static inline void load_block(uint64_t b[2], const uint8_t *m) {
    b[0] = load64_le(m);
    b[1] = (load64_le(m + 7) >> 8) | (0x01ULL << 56);
}

static inline void add_to_acc(uint64_t acc[3], const uint64_t b[2]) {
    uint64_t c = adc(&acc[0], b[0]);
    c = adc(&acc[1], b[1] + c);
    acc[2] += c;
}

static void add_block_partial(uint64_t acc[3], const uint8_t *m, size_t len) {
    uint8_t pad[16] = {0};
    memcpy(pad, m, len);
    pad[len] = 0x01;
    uint64_t b[2] = {load64_le(pad), load64_le(pad + 8)};
    add_to_acc(acc, b);
}

/* KEY CLAMPING */

static void clamp_r(uint8_t r[16]) {
    r[3]  &= 0x0F;
    r[7]  &= 0x0F;
    r[11] &= 0x0F;
    r[15] &= 0x07;
    r[4]  &= 0xFC;
    r[8]  &= 0xFC;
    r[12] &= 0xFC;
}

/* MULTI-BLOCK PROCESSING */

static void process_2blocks(uint64_t acc[3], const uint8_t *m1, const uint8_t *m2,
    const uint64_t r[2], const uint64_t r2[2]) {
    uint64_t b1[2], b2[2];
    load_block(b1, m1);
    load_block(b2, m2);

    add_to_acc(acc, b1);
    mul_reduce(acc, acc, r2);

    uint64_t p2[3];
    mul_2x2(p2, b2, r);
    uint64_t c = adc(&acc[0], p2[0]);
    c = adc(&acc[1], p2[1] + c);
    acc[2] += p2[2] + c;
}

static void process_4blocks(uint64_t acc[3], const uint8_t *m1, const uint8_t *m2,
    const uint8_t *m3, const uint8_t *m4,
    const uint64_t r[2], const uint64_t r2[2],
    const uint64_t r3[2], const uint64_t r4[2]) {
    uint64_t b1[2], b2[2], b3[2], b4[2];
    load_block(b1, m1);
    load_block(b2, m2);
    load_block(b3, m3);
    load_block(b4, m4);

    add_to_acc(acc, b1);
    mul_reduce(acc, acc, r4);

    uint64_t p2[3], p3[3], p4[3];
    mul_2x2(p2, b2, r3);
    mul_2x2(p3, b3, r2);
    mul_2x2(p4, b4, r);

    uint64_t c = adc(&acc[0], p2[0]);
    c = adc(&acc[1], p2[1] + c);
    acc[2] += p2[2] + c;
    c = adc(&acc[0], p3[0]);
    c = adc(&acc[1], p3[1] + c);
    acc[2] += p3[2] + c;
    c = adc(&acc[0], p4[0]);
    c = adc(&acc[1], p4[1] + c);
    acc[2] += p4[2] + c;
}

/* FINALIZATION */

static void finalize(uint64_t acc[3], const uint64_t s[2], uint8_t tag[16]) {
    reduce_full(acc);

    uint64_t eq_lo = ~(acc[0] ^ M127_LO);
    eq_lo &= eq_lo >> 32;
    eq_lo &= eq_lo >> 16;
    eq_lo &= eq_lo >> 8;
    eq_lo &= eq_lo >> 4;
    eq_lo &= eq_lo >> 2;
    eq_lo &= eq_lo >> 1;

    uint64_t eq_hi = ~(acc[1] ^ M127_HI);
    eq_hi &= eq_hi >> 32;
    eq_hi &= eq_hi >> 16;
    eq_hi &= eq_hi >> 8;
    eq_hi &= eq_hi >> 4;
    eq_hi &= eq_hi >> 2;
    eq_hi &= eq_hi >> 1;

    uint64_t wrap = (uint64_t)0 - (eq_lo & eq_hi & 1);
    uint64_t c = adc(&acc[0], wrap & 1);
    acc[1] = (acc[1] + c) & M127_HI;

    uint64_t t0 = acc[0] + s[0];
    c = (t0 < acc[0]);
    uint64_t t1 = acc[1] + s[1] + c;

    store64_le(tag, t0);
    store64_le(tag + 8, t1);
}

#define LAZY_INTERVAL 8

/* PUBLIC API */

void cfx_poly1271_init(cfx_poly1271_ctx_t *ctx, const uint8_t key[32]) {
    uint8_t rc[16];
    memcpy(rc, key, 16);
    clamp_r(rc);

    ctx->r[0] = load64_le(rc);
    ctx->r[1] = load64_le(rc + 8);

    square_mod(ctx->r2, ctx->r);
    mul_mod(ctx->r3, ctx->r2, ctx->r);
    mul_mod(ctx->r4, ctx->r2, ctx->r2);

    ctx->s[0] = load64_le(key + 16);
    ctx->s[1] = load64_le(key + 24);

    ctx->acc[0] = ctx->acc[1] = ctx->acc[2] = 0;
    ctx->buflen = ctx->blkcnt = 0;

    CFX_MEMZERO_S(rc, 16);
}

void cfx_poly1271_update(cfx_poly1271_ctx_t *ctx, const uint8_t *msg, size_t len) {
    if (ctx->buflen) {
        size_t need = CFX_POLY1271_BLOCK_SIZE - ctx->buflen;
        if (len < need) {
            memcpy(ctx->buf + ctx->buflen, msg, len);
            ctx->buflen += (uint8_t)len;
            return;
        }
        memcpy(ctx->buf + ctx->buflen, msg, need);

        uint64_t b[2];
        load_block(b, ctx->buf);
        add_to_acc(ctx->acc, b);
        mul_reduce(ctx->acc, ctx->acc, ctx->r);

        if (++ctx->blkcnt >= LAZY_INTERVAL) {
            reduce_full(ctx->acc);
            ctx->blkcnt = 0;
        }

        msg += need;
        len -= need;
        ctx->buflen = 0;
    }

    while (len >= 4 * CFX_POLY1271_BLOCK_SIZE) {
        process_4blocks(ctx->acc, msg,
            msg + CFX_POLY1271_BLOCK_SIZE,
            msg + 2 * CFX_POLY1271_BLOCK_SIZE,
            msg + 3 * CFX_POLY1271_BLOCK_SIZE,
            ctx->r, ctx->r2, ctx->r3, ctx->r4);
        ctx->blkcnt += 4;
        if (ctx->blkcnt >= LAZY_INTERVAL) {
            reduce_full(ctx->acc);
            ctx->blkcnt = 0;
        }
        msg += 4 * CFX_POLY1271_BLOCK_SIZE;
        len -= 4 * CFX_POLY1271_BLOCK_SIZE;
    }

    while (len >= 2 * CFX_POLY1271_BLOCK_SIZE) {
        process_2blocks(ctx->acc, msg, msg + CFX_POLY1271_BLOCK_SIZE, ctx->r, ctx->r2);
        ctx->blkcnt += 2;
        if (ctx->blkcnt >= LAZY_INTERVAL) {
            reduce_full(ctx->acc);
            ctx->blkcnt = 0;
        }
        msg += 2 * CFX_POLY1271_BLOCK_SIZE;
        len -= 2 * CFX_POLY1271_BLOCK_SIZE;
    }

    while (len >= CFX_POLY1271_BLOCK_SIZE) {
        uint64_t b[2];
        load_block(b, msg);
        add_to_acc(ctx->acc, b);
        mul_reduce(ctx->acc, ctx->acc, ctx->r);
        if (++ctx->blkcnt >= LAZY_INTERVAL) {
            reduce_full(ctx->acc);
            ctx->blkcnt = 0;
        }
        msg += CFX_POLY1271_BLOCK_SIZE;
        len -= CFX_POLY1271_BLOCK_SIZE;
    }

    if (len) {
        memcpy(ctx->buf, msg, len);
        ctx->buflen = (uint8_t)len;
    }
}

void cfx_poly1271_finish(cfx_poly1271_ctx_t *ctx, uint8_t tag[16]) {
    if (ctx->buflen) {
        add_block_partial(ctx->acc, ctx->buf, ctx->buflen);
        mul_reduce(ctx->acc, ctx->acc, ctx->r);
    }

    finalize(ctx->acc, ctx->s, tag);
    CFX_MEMZERO_S(ctx, sizeof(*ctx));
}

void cfx_poly1271(uint8_t tag[16], const uint8_t *msg, size_t len, const uint8_t key[32]) {
    cfx_poly1271_ctx_t ctx;
    cfx_poly1271_init(&ctx, key);
    cfx_poly1271_update(&ctx, msg, len);
    cfx_poly1271_finish(&ctx, tag);
}

int cfx_poly1271_verify(const uint8_t tag[16], const uint8_t *msg, size_t len, const uint8_t key[32]) {
    uint8_t computed[16];
    cfx_poly1271(computed, msg, len, key);

    uint8_t diff = 0;
    for (int i = 0; i < 16; i++) {
        diff |= computed[i] ^ tag[i];
    }

    CFX_MEMZERO_S(computed, 16);
    return diff ? -1 : 0;
}
