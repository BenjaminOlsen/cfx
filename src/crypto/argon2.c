/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/* argon2.c - Argon2 password hashing (RFC 9106) */

#include "cfx/argon2.h"
#include "cfx/blake2.h"
#include "cfx/base64.h"
#include "cfx/memory.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define BLOCK_SZ    1024
#define BLOCK_QW    128         /* BLOCK_SZ / sizeof(uint64_t) */
#define SYNC_POINTS 4
#define VERSION     0x13        /* RFC 9106 v19 */
#define MIN_OUTLEN  4
#define MIN_SALT    8
#define MIN_MEM     8
#define MIN_TIME    1
#define MIN_LANES   1

typedef struct {
    uint64_t v[BLOCK_QW];
} block_t;

typedef struct {
    block_t  *mem;
    uint32_t  m, t, lanes;
    uint32_t  lane_len, seg_len;
    int       type;
} ctx_t;

static inline uint64_t rotr64(uint64_t x, unsigned n) {
    return (x >> n) | (x << (64 - n));
}

/* H' variable-length hash (RFC 9106 §3.2) */
static void hash_long(uint8_t *out, size_t n, const uint8_t *in, size_t inlen) {
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
    out += 32;
    n   -= 32;

    while (n > 64) {
        cfx_blake2b_init(&ctx, 64);
        cfx_blake2b_update(&ctx, v, 64);
        cfx_blake2b_final(&ctx, v);
        memcpy(out, v, 32);
        out += 32;
        n   -= 32;
    }

    cfx_blake2b_init(&ctx, n);
    cfx_blake2b_update(&ctx, v, 64);
    cfx_blake2b_final(&ctx, out);
}

/* GB: argon2 mixing function (modified blake2b G with mul) */
#define GB(a, b, c, d) do {                                             \
    a = a + b + 2 * (uint64_t)(uint32_t)a * (uint64_t)(uint32_t)b;     \
    d = rotr64(d ^ a, 32);                                              \
    c = c + d + 2 * (uint64_t)(uint32_t)c * (uint64_t)(uint32_t)d;     \
    b = rotr64(b ^ c, 24);                                              \
    a = a + b + 2 * (uint64_t)(uint32_t)a * (uint64_t)(uint32_t)b;     \
    d = rotr64(d ^ a, 16);                                              \
    c = c + d + 2 * (uint64_t)(uint32_t)c * (uint64_t)(uint32_t)d;     \
    b = rotr64(b ^ c, 63);                                              \
} while (0)

static void permute(uint64_t *v) {
    /* two rounds of blake2b mixing */
    for (int r = 0; r < 2; r++) {
        GB(v[0], v[4], v[8],  v[12]);
        GB(v[1], v[5], v[9],  v[13]);
        GB(v[2], v[6], v[10], v[14]);
        GB(v[3], v[7], v[11], v[15]);
        GB(v[0], v[5], v[10], v[15]);
        GB(v[1], v[6], v[11], v[12]);
        GB(v[2], v[7], v[8],  v[13]);
        GB(v[3], v[4], v[9],  v[14]);
    }
}

static void fill_block(const block_t *prev, const block_t *ref,
                        block_t *next, int xor) {
    block_t tmp;
    uint64_t r[BLOCK_QW];

    for (int i = 0; i < BLOCK_QW; i++) {
        r[i] = prev->v[i] ^ ref->v[i];
    }
    memcpy(&tmp, r, sizeof(tmp));

    /* row-wise permutation */
    for (int i = 0; i < 8; i++)
        permute(&tmp.v[i * 16]);

    /* column-wise permutation */
    for (int i = 0; i < 8; i++) {
        uint64_t col[16];
        for (int j = 0; j < 16; j++) col[j] = tmp.v[j * 8 + i];
        permute(col);
        for (int j = 0; j < 16; j++) tmp.v[j * 8 + i] = col[j];
    }

    if (xor) {
        for (int i = 0; i < BLOCK_QW; i++)
            next->v[i] ^= r[i] ^ tmp.v[i];
    } else {
        for (int i = 0; i < BLOCK_QW; i++)
            next->v[i] = r[i] ^ tmp.v[i];
    }
}

static uint32_t index_alpha(const ctx_t *c, uint32_t pass, uint32_t slice,
                            uint32_t idx, uint64_t rand, int same_lane) {
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

static void fill_segment(ctx_t *c, uint32_t pass, uint32_t lane, uint32_t slice) {
    uint32_t si  = (pass == 0 && slice == 0) ? 2 : 0;
    uint32_t cur = lane * c->lane_len + slice * c->seg_len + si;
    uint32_t prv = cur - 1;

    if (cur % c->lane_len == 0)
        prv += c->lane_len;

    int di = (c->type == CFX_ARGON2I) ||
             (c->type == CFX_ARGON2ID && pass == 0 && slice < SYNC_POINTS / 2);

    block_t addr, zero, input;
    memset(&zero, 0, sizeof(zero));
    memset(&input, 0, sizeof(input));

    if (di) {
        input.v[0] = pass;
        input.v[1] = lane;
        input.v[2] = slice;
        input.v[3] = c->lane_len * c->lanes;
        input.v[4] = c->t;
        input.v[5] = c->type;
    }

    uint32_t ai = 0;

    for (uint32_t i = si; i < c->seg_len; i++, cur++, prv++) {
        if (cur % c->lane_len == 1)
            prv = cur - 1;

        uint64_t rand;
        if (di) {
            if (ai == 0) {
                input.v[6]++;
                fill_block(&zero, &input, &addr, 0);
                fill_block(&zero, &addr, &addr, 0);
            }
            rand = addr.v[ai];
            ai = (ai + 1) % BLOCK_QW;
        } else {
            rand = c->mem[prv].v[0];
        }

        uint32_t rl = (uint32_t)(rand >> 32) % c->lanes;
        if (pass == 0 && slice == 0) rl = lane;

        uint32_t ri = index_alpha(c, pass, slice, i, rand, rl == lane);
        uint32_t ro = rl * c->lane_len + ri;

        fill_block(&c->mem[prv], &c->mem[ro], &c->mem[cur], pass > 0);
    }
}

/* H0: initial 64-byte hash from all parameters (RFC 9106 §3.3) */
static void initial_hash(uint8_t *h0,
                         const uint8_t *pwd, size_t pwdlen,
                         const uint8_t *salt, size_t saltlen,
                         uint32_t m, uint32_t t, uint32_t p,
                         size_t outlen, int type) {
    cfx_blake2b_ctx_t ctx;
    uint8_t buf[4];

    cfx_blake2b_init(&ctx, 64);

    cfx_store32_le(buf, p);            
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_store32_le(buf, (uint32_t)outlen);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_store32_le(buf, m);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_store32_le(buf, t);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_store32_le(buf, VERSION);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_store32_le(buf, (uint32_t)type);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_store32_le(buf, (uint32_t)pwdlen);
    cfx_blake2b_update(&ctx, buf, 4);
    if (pwdlen > 0) cfx_blake2b_update(&ctx, pwd, pwdlen);

    cfx_store32_le(buf, (uint32_t)saltlen);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_blake2b_update(&ctx, salt, saltlen);

    /* secret + AD not supported (length = 0) */
    cfx_store32_le(buf, 0);
    cfx_blake2b_update(&ctx, buf, 4);
    cfx_blake2b_update(&ctx, buf, 4);

    cfx_blake2b_final(&ctx, h0);
}

int cfx_argon2(uint8_t *out, size_t outlen,
               const uint8_t *pwd, size_t pwdlen,
               const uint8_t *salt, size_t saltlen,
               uint32_t m, uint32_t t, uint32_t p, int type) {

    if (!out || outlen < MIN_OUTLEN)         return -1;
    if (!salt || saltlen < MIN_SALT)         return -1;
    if (t < MIN_TIME || p < MIN_LANES)       return -1;
    if (type < CFX_ARGON2D || type > CFX_ARGON2ID) return -1;

    uint32_t nblocks = m;
    if (nblocks < 8 * p) nblocks = 8 * p;

    uint32_t seg_len = nblocks / (4 * p);
    nblocks = seg_len * 4 * p;

    ctx_t c;
    c.mem = (block_t *)calloc(nblocks, sizeof(block_t));
    if (!c.mem) return -1;

    c.m        = m;
    c.t        = t;
    c.lanes    = p;
    c.lane_len = nblocks / p;
    c.seg_len  = seg_len;
    c.type     = type;

    /* H0 */
    uint8_t h0[64];
    initial_hash(h0, pwd, pwdlen, salt, saltlen, m, t, p, outlen, type);

    /* fill first two blocks of each lane */
    uint8_t bhi[72];   /* block hash input = H0 || LE32(index) || LE32(lane) */
    memcpy(bhi, h0, 64);

    for (uint32_t lane = 0; lane < p; lane++) {
        cfx_store32_le(bhi + 64, 0);
        cfx_store32_le(bhi + 68, lane);
        hash_long((uint8_t *)c.mem[lane * c.lane_len].v, BLOCK_SZ, bhi, 72);

        cfx_store32_le(bhi + 64, 1);
        hash_long((uint8_t *)c.mem[lane * c.lane_len + 1].v, BLOCK_SZ, bhi, 72);
    }

    /* fill remaining blocks */
    for (uint32_t pass = 0; pass < t; pass++)
        for (uint32_t slice = 0; slice < SYNC_POINTS; slice++)
            for (uint32_t lane = 0; lane < p; lane++)
                fill_segment(&c, pass, lane, slice);

    /* XOR last block of each lane → final block → H' → output */
    block_t fin;
    memcpy(&fin, &c.mem[c.lane_len - 1], sizeof(block_t));
    for (uint32_t lane = 1; lane < p; lane++) {
        uint32_t last = lane * c.lane_len + c.lane_len - 1;
        for (int i = 0; i < BLOCK_QW; i++)
            fin.v[i] ^= c.mem[last].v[i];
    }

    hash_long(out, outlen, (uint8_t *)fin.v, BLOCK_SZ);

    cfx_memzero_s(c.mem, nblocks * sizeof(block_t));
    cfx_memzero_s(h0, sizeof(h0));
    cfx_memzero_s(&fin, sizeof(fin));
    free(c.mem);

    return 0;
}

int cfx_argon2id(uint8_t *out, size_t outlen, const uint8_t *pwd, size_t pwdlen,
                 const uint8_t *salt, size_t saltlen, uint32_t m, uint32_t t, uint32_t p) {
    return cfx_argon2(out, outlen, pwd, pwdlen, salt, saltlen, m, t, p, CFX_ARGON2ID);
}

int cfx_argon2d(uint8_t *out, size_t outlen, const uint8_t *pwd, size_t pwdlen,
                const uint8_t *salt, size_t saltlen, uint32_t m, uint32_t t, uint32_t p) {
    return cfx_argon2(out, outlen, pwd, pwdlen, salt, saltlen, m, t, p, CFX_ARGON2D);
}

int cfx_argon2i(uint8_t *out, size_t outlen, const uint8_t *pwd, size_t pwdlen,
                const uint8_t *salt, size_t saltlen, uint32_t m, uint32_t t, uint32_t p) {
    return cfx_argon2(out, outlen, pwd, pwdlen, salt, saltlen, m, t, p, CFX_ARGON2I);
}

/* ── PHC string encoding ─────────────────────────────────────────── */

static const char *type_str(int type) {
    switch (type) {
    case CFX_ARGON2D:  return "argon2d";
    case CFX_ARGON2I:  return "argon2i";
    case CFX_ARGON2ID: return "argon2id";
    default:           return "argon2";
    }
}

int cfx_argon2_encode(char *out, size_t outlen,
                      const uint8_t *hash, size_t hashlen,
                      const uint8_t *salt, size_t saltlen,
                      uint32_t m, uint32_t t, uint32_t p, int type) {

    size_t sb64 = cfx_base64_enc_len(saltlen);
    size_t hb64 = cfx_base64_enc_len(hashlen);
    size_t need = 1 + strlen(type_str(type)) + 1 + 5 + 3
                + 2 + 10 + 1 + 2 + 10 + 1 + 2 + 10 + 1
                + sb64 + 1 + hb64 + 1;

    if (outlen < need) return -1;

    char sb[256], hb[256];
    size_t sn = sizeof(sb), hn = sizeof(hb);
    cfx_base64_encode(sb, &sn, salt, saltlen);
    cfx_base64_encode(hb, &hn, hash, hashlen);

    /* strip padding (PHC uses unpadded base64) */
    while (sn > 0 && sb[sn - 1] == '=') sn--;
    while (hn > 0 && hb[hn - 1] == '=') hn--;
    sb[sn] = '\0';
    hb[hn] = '\0';

    return snprintf(out, outlen, "$%s$v=%d$m=%u,t=%u,p=%u$%s$%s",
                    type_str(type), VERSION, m, t, p, sb, hb);
}

static int b64_decode_unpadded(uint8_t *out, size_t *outlen,
                               const char *in, size_t inlen) {
    char padded[512];
    if (inlen >= sizeof(padded) - 4) return -1;

    memcpy(padded, in, inlen);
    size_t pad = (4 - (inlen % 4)) % 4;
    for (size_t i = 0; i < pad; i++) padded[inlen + i] = '=';
    padded[inlen + pad] = '\0';

    return cfx_base64_decode(out, outlen, padded, inlen + pad);
}

int cfx_argon2_verify(const char *enc, const uint8_t *pwd, size_t pwdlen) {
    if (!enc || enc[0] != '$') return -1;

    int type;
    if      (strncmp(enc + 1, "argon2id$", 9) == 0) { type = CFX_ARGON2ID; enc += 10; }
    else if (strncmp(enc + 1, "argon2i$",  8) == 0) { type = CFX_ARGON2I;  enc += 9;  }
    else if (strncmp(enc + 1, "argon2d$",  8) == 0) { type = CFX_ARGON2D;  enc += 9;  }
    else return -1;

    int ver;
    if (sscanf(enc, "v=%d$", &ver) != 1) return -1;
    enc = strchr(enc, '$');
    if (!enc) return -1;
    enc++;

    uint32_t m, t, p;
    if (sscanf(enc, "m=%u,t=%u,p=%u$", &m, &t, &p) != 3) return -1;
    enc = strchr(enc, '$');
    if (!enc) return -1;
    enc++;

    /* decode salt */
    const char *sep = strchr(enc, '$');
    if (!sep) return -1;
    size_t sb64len = sep - enc;

    uint8_t salt[256];
    size_t saltlen = sizeof(salt);
    if (b64_decode_unpadded(salt, &saltlen, enc, sb64len) != 0)
        return -1;

    enc = sep + 1;

    /* decode stored hash */
    uint8_t stored[256];
    size_t hashlen = sizeof(stored);
    if (b64_decode_unpadded(stored, &hashlen, enc, strlen(enc)) != 0)
        return -1;

    /* recompute */
    uint8_t computed[256];
    if (cfx_argon2(computed, hashlen, pwd, pwdlen, salt, saltlen, m, t, p, type) != 0)
        return -1;

    /* constant-time compare */
    uint8_t diff = 0;
    for (size_t i = 0; i < hashlen; i++)
        diff |= stored[i] ^ computed[i];

    cfx_memzero_s(computed, sizeof(computed));

    return diff != 0;   /* 0 = match, 1 = mismatch */
}
