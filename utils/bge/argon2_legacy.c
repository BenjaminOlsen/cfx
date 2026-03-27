/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * argon2_legacy.c — Argon2id with the bugs present in cfx <= v2.6.0
 *
 * This is intentionally buggy code kept for backward compatibility with
 * existing BGE v2/v4 stores.  Three bugs are preserved:
 *   1. permute() applies two rounds instead of one
 *   2. Column extraction uses stride 8 instead of stride-16 pairs
 *   3. index_alpha has a spurious same_lane bonus of +seg_len
 *   4. Data-independent address block is not pre-generated
 *
 * New stores (v5+) must use cfx_argon2id() from the library.
 * DO NOT use this for anything other than legacy store decryption.
 */

#include "cfx/blake2.h"
#include "cfx/memory.h"
#include "cfx/bits.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define BLOCK_SZ    1024
#define BLOCK_QW    128
#define SYNC_POINTS 4
#define VERSION     0x13

typedef struct { uint64_t v[BLOCK_QW]; } leg_block_t;

typedef struct {
    leg_block_t *mem;
    uint32_t m, t, lanes;
    uint32_t lane_len, seg_len;
} leg_ctx_t;

#define rotr64 cfx_rotr64

#define LGB(a, b, c, d) do {                                            \
    a = a + b + 2 * (uint64_t)(uint32_t)a * (uint64_t)(uint32_t)b;     \
    d = rotr64(d ^ a, 32);                                              \
    c = c + d + 2 * (uint64_t)(uint32_t)c * (uint64_t)(uint32_t)d;     \
    b = rotr64(b ^ c, 24);                                              \
    a = a + b + 2 * (uint64_t)(uint32_t)a * (uint64_t)(uint32_t)b;     \
    d = rotr64(d ^ a, 16);                                              \
    c = c + d + 2 * (uint64_t)(uint32_t)c * (uint64_t)(uint32_t)d;     \
    b = rotr64(b ^ c, 63);                                              \
} while (0)

/* BUG 1: two rounds instead of one */
static void leg_permute(uint64_t *v) {
    LGB(v[0], v[4], v[8],  v[12]);
    LGB(v[1], v[5], v[9],  v[13]);
    LGB(v[2], v[6], v[10], v[14]);
    LGB(v[3], v[7], v[11], v[15]);
    LGB(v[0], v[5], v[10], v[15]);
    LGB(v[1], v[6], v[11], v[12]);
    LGB(v[2], v[7], v[8],  v[13]);
    LGB(v[3], v[4], v[9],  v[14]);

    LGB(v[0], v[4], v[8],  v[12]);
    LGB(v[1], v[5], v[9],  v[13]);
    LGB(v[2], v[6], v[10], v[14]);
    LGB(v[3], v[7], v[11], v[15]);
    LGB(v[0], v[5], v[10], v[15]);
    LGB(v[1], v[6], v[11], v[12]);
    LGB(v[2], v[7], v[8],  v[13]);
    LGB(v[3], v[4], v[9],  v[14]);
}

static void leg_fill_block(const leg_block_t *prev, const leg_block_t *ref,
                            leg_block_t *next, int do_xor) {
    leg_block_t tmp;
    uint64_t r[BLOCK_QW];

    for (int i = 0; i < BLOCK_QW; i++)
        r[i] = prev->v[i] ^ ref->v[i];
    memcpy(&tmp, r, sizeof(tmp));

    for (int i = 0; i < 8; i++)
        leg_permute(&tmp.v[i * 16]);

    /* BUG 2: wrong column stride (j*8+i instead of stride-16 pairs) */
    for (int i = 0; i < 8; i++) {
        uint64_t col[16];
        for (int j = 0; j < 16; j++)
            col[j] = tmp.v[j * 8 + i];
        leg_permute(col);
        for (int j = 0; j < 16; j++)
            tmp.v[j * 8 + i] = col[j];
    }

    if (do_xor) {
        for (int i = 0; i < BLOCK_QW; i++)
            next->v[i] ^= r[i] ^ tmp.v[i];
    } else {
        for (int i = 0; i < BLOCK_QW; i++)
            next->v[i] = r[i] ^ tmp.v[i];
    }

    cfx_memzero_s(r, sizeof(r));
    cfx_memzero_s(&tmp, sizeof(tmp));
}

/* BUG 3: same_lane bonus of +seg_len */
static uint32_t leg_index_alpha(const leg_ctx_t *c, uint32_t pass,
                                 uint32_t slice, uint32_t idx,
                                 uint64_t rand, int same_lane) {
    uint32_t area, start;
    if (pass == 0) {
        area = (slice == 0)
            ? idx - 1
            : slice * c->seg_len + idx - 1 + (same_lane ? c->seg_len : 0);
        start = 0;
    } else {
        area = c->lane_len - c->seg_len + idx - 1
             + (same_lane ? c->seg_len : 0);
        start = (slice + 1) * c->seg_len;
        if (start >= c->lane_len) start = 0;
    }
    uint64_t rel = rand & 0xFFFFFFFF;
    rel = (rel * rel) >> 32;
    rel = area - 1 - ((area * rel) >> 32);
    return (start + (uint32_t)rel) % c->lane_len;
}

static void leg_fill_segment(leg_ctx_t *c, uint32_t pass, uint32_t lane,
                              uint32_t slice) {
    uint32_t si  = (pass == 0 && slice == 0) ? 2 : 0;
    uint32_t cur = lane * c->lane_len + slice * c->seg_len + si;
    uint32_t prv = cur - 1;
    if (cur % c->lane_len == 0) prv += c->lane_len;

    int di = (pass == 0 && slice < SYNC_POINTS / 2);

    leg_block_t addr, zero, input;
    memset(&zero, 0, sizeof(zero));
    memset(&input, 0, sizeof(input));

    if (di) {
        input.v[0] = pass;
        input.v[1] = lane;
        input.v[2] = slice;
        input.v[3] = c->lane_len * c->lanes;
        input.v[4] = c->t;
        input.v[5] = 2; /* CFX_ARGON2ID */
    }

    /* BUG 4: no pre-generation of address block */
    uint32_t ai = 0;

    for (uint32_t i = si; i < c->seg_len; i++, cur++, prv++) {
        if (cur % c->lane_len == 1) prv = cur - 1;

        uint64_t rand;
        if (di) {
            if (ai == 0) {
                input.v[6]++;
                leg_fill_block(&zero, &input, &addr, 0);
                leg_fill_block(&zero, &addr, &addr, 0);
            }
            rand = addr.v[ai];
            ai = (ai + 1) % BLOCK_QW;
        } else {
            rand = c->mem[prv].v[0];
        }

        uint32_t rl = (uint32_t)(rand >> 32) % c->lanes;
        if (pass == 0 && slice == 0) rl = lane;

        uint32_t ri = leg_index_alpha(c, pass, slice, i, rand, rl == lane);
        uint32_t ro = rl * c->lane_len + ri;

        leg_fill_block(&c->mem[prv], &c->mem[ro], &c->mem[cur], pass > 0);
    }
}

static void leg_hash_long(uint8_t *out, size_t n,
                           const uint8_t *in, size_t inlen) {
    uint8_t nb[4];
    cfx_store32_le(nb, (uint32_t)n);

    if (n <= 64) {
        cfx_blake2b_ctx_t ctx;
        cfx_blake2b_init(&ctx, n);
        cfx_blake2b_update(&ctx, nb, 4);
        cfx_blake2b_update(&ctx, in, inlen);
        cfx_blake2b_final(&ctx, out);
        return;
    }

    uint8_t v[64];
    cfx_blake2b_ctx_t ctx;
    cfx_blake2b_init(&ctx, 64);
    cfx_blake2b_update(&ctx, nb, 4);
    cfx_blake2b_update(&ctx, in, inlen);
    cfx_blake2b_final(&ctx, v);
    memcpy(out, v, 32);
    out += 32; n -= 32;

    while (n > 64) {
        cfx_blake2b_init(&ctx, 64);
        cfx_blake2b_update(&ctx, v, 64);
        cfx_blake2b_final(&ctx, v);
        memcpy(out, v, 32);
        out += 32; n -= 32;
    }

    cfx_blake2b_init(&ctx, n);
    cfx_blake2b_update(&ctx, v, 64);
    cfx_blake2b_final(&ctx, out);

    cfx_memzero_s(v, sizeof(v));
}

int bge_argon2_legacy(uint8_t *out, size_t outlen,
                       const uint8_t *pwd, size_t pwdlen,
                       const uint8_t *salt, size_t saltlen,
                       uint32_t m, uint32_t t, uint32_t p) {
    if (!out || outlen < 4) return -1;
    if (!salt || saltlen < 8) return -1;
    if (t < 1 || p < 1) return -1;

    uint32_t nblocks = m;
    if (nblocks < 8 * p) nblocks = 8 * p;
    uint32_t seg_len = nblocks / (4 * p);
    nblocks = seg_len * 4 * p;

    leg_ctx_t c;
    c.mem = (leg_block_t *)calloc(nblocks, sizeof(leg_block_t));
    if (!c.mem) return -1;

    c.m = m; c.t = t; c.lanes = p;
    c.lane_len = nblocks / p;
    c.seg_len = seg_len;

    /* H0 — identical to correct implementation */
    uint8_t h0[64];
    {
        cfx_blake2b_ctx_t ctx;
        uint8_t buf[4];
        cfx_blake2b_init(&ctx, 64);
        cfx_store32_le(buf, p);             cfx_blake2b_update(&ctx, buf, 4);
        cfx_store32_le(buf, (uint32_t)outlen); cfx_blake2b_update(&ctx, buf, 4);
        cfx_store32_le(buf, m);             cfx_blake2b_update(&ctx, buf, 4);
        cfx_store32_le(buf, t);             cfx_blake2b_update(&ctx, buf, 4);
        cfx_store32_le(buf, VERSION);       cfx_blake2b_update(&ctx, buf, 4);
        cfx_store32_le(buf, 2);             cfx_blake2b_update(&ctx, buf, 4);
        cfx_store32_le(buf, (uint32_t)pwdlen); cfx_blake2b_update(&ctx, buf, 4);
        if (pwdlen > 0) cfx_blake2b_update(&ctx, pwd, pwdlen);
        cfx_store32_le(buf, (uint32_t)saltlen); cfx_blake2b_update(&ctx, buf, 4);
        cfx_blake2b_update(&ctx, salt, saltlen);
        cfx_store32_le(buf, 0); cfx_blake2b_update(&ctx, buf, 4);
        cfx_store32_le(buf, 0); cfx_blake2b_update(&ctx, buf, 4);
        cfx_blake2b_final(&ctx, h0);
    }

    uint8_t bhi[72];
    memcpy(bhi, h0, 64);
    for (uint32_t lane = 0; lane < p; lane++) {
        cfx_store32_le(bhi + 64, 0);
        cfx_store32_le(bhi + 68, lane);
        leg_hash_long((uint8_t *)c.mem[lane * c.lane_len].v, BLOCK_SZ, bhi, 72);
        cfx_store32_le(bhi + 64, 1);
        leg_hash_long((uint8_t *)c.mem[lane * c.lane_len + 1].v, BLOCK_SZ, bhi, 72);
    }

    for (uint32_t pass = 0; pass < t; pass++) {
        for (uint32_t slice = 0; slice < SYNC_POINTS; slice++) {
            for (uint32_t lane = 0; lane < p; lane++) {
                leg_fill_segment(&c, pass, lane, slice);
            }
        }
    }

    leg_block_t fin;
    memcpy(&fin, &c.mem[c.lane_len - 1], sizeof(leg_block_t));
    for (uint32_t lane = 1; lane < p; lane++) {
        uint32_t last = lane * c.lane_len + c.lane_len - 1;
        for (int i = 0; i < BLOCK_QW; i++)
            fin.v[i] ^= c.mem[last].v[i];
    }

    leg_hash_long(out, outlen, (uint8_t *)fin.v, BLOCK_SZ);

    cfx_memzero_s(c.mem, nblocks * sizeof(leg_block_t));
    cfx_memzero_s(h0, sizeof(h0));
    cfx_memzero_s(&fin, sizeof(fin));
    free(c.mem);
    return 0;
}
