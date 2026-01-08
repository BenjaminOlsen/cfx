/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * sha3.c - SHA-3 (Keccak) implementation (FIPS 202)
 */

#include "cfx/sha3.h"
#include "cfx/macros.h"
#include <string.h>

/* internal state */
typedef struct {
    uint64_t state[25];  /* 1600-bit state as 5x5 lanes */
    uint8_t  buf[168];   /* max rate = 1344 bits = 168 bytes (SHAKE128) */
    size_t   buflen;
    size_t   rate;       /* rate in bytes */
    uint8_t  domain;     /* domain separator: 0x06 for SHA3, 0x1F for SHAKE */
    uint8_t  finalized;  /* 1 if in squeeze phase */
} sha3_state_t;

CFX_STATIC_ASSERT(sizeof(sha3_state_t) <= CFX_SHA3_CTX_SIZE, sha3_ctx_too_small);

/* Keccak round constants */
static const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL,
    0x800000000000808aULL, 0x8000000080008000ULL,
    0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008aULL, 0x0000000000000088ULL,
    0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL,
    0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL,
    0x8000000080008081ULL, 0x8000000000008080ULL,
    0x0000000080000001ULL, 0x8000000080008008ULL
};

/* rotation offsets for rho step */
static const int RHO[25] = {
     0,  1, 62, 28, 27,
    36, 44,  6, 55, 20,
     3, 10, 43, 25, 39,
    41, 45, 15, 21,  8,
    18,  2, 61, 56, 14
};

/* pi step permutation: pi[x + 5*y] = source index */
static const int PI[25] = {
     0, 6, 12, 18, 24,
     3, 9, 10, 16, 22,
     1, 7, 13, 19, 20,
     4, 5, 11, 17, 23,
     2, 8, 14, 15, 21
};

static inline uint64_t rotl64(uint64_t x, int n) {
    return (x << n) | (x >> (64 - n));
}

static inline uint64_t load64_le(const void* src) {
    const uint8_t* p = (const uint8_t*)src;
    return ((uint64_t)p[0])       | ((uint64_t)p[1] << 8)  |
           ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}


/*
 * Keccak-f[1600] permutation - 24 rounds
 */
static void keccak_f1600(uint64_t state[25]) {
    uint64_t C[5], D[5], B[25];

    for (int round = 0; round < 24; round++) {
        /* theta: column parity */
        for (int x = 0; x < 5; x++) {
            C[x] = state[x] ^ state[x + 5] ^ state[x + 10] ^ state[x + 15] ^ state[x + 20];
        }
        for (int x = 0; x < 5; x++) {
            D[x] = C[(x + 4) % 5] ^ rotl64(C[(x + 1) % 5], 1);
        }
        for (int i = 0; i < 25; i++) {
            state[i] ^= D[i % 5];
        }

        /* rho and pi combined */
        for (int i = 0; i < 25; i++) {
            B[i] = rotl64(state[PI[i]], RHO[PI[i]]);
        }

        /* chi: non-linear mixing */
        for (int y = 0; y < 5; y++) {
            for (int x = 0; x < 5; x++) {
                int i = x + 5 * y;
                state[i] = B[i] ^ ((~B[(x + 1) % 5 + 5 * y]) & B[(x + 2) % 5 + 5 * y]);
            }
        }

        /* iota: round constant */
        state[0] ^= RC[round];
    }
}

/*
 * Absorb data into sponge
 */
static void sha3_absorb(sha3_state_t* S, const uint8_t* data, size_t len) {
    size_t rate = S->rate;

    while (len > 0) {
        size_t chunk = rate - S->buflen;
        if (chunk > len) chunk = len;

        memcpy(S->buf + S->buflen, data, chunk);
        S->buflen += chunk;
        data += chunk;
        len -= chunk;

        if (S->buflen == rate) {
            /* XOR block into state */
            for (size_t i = 0; i < rate / 8; i++) {
                S->state[i] ^= load64_le(S->buf + i * 8);
            }
            keccak_f1600(S->state);
            S->buflen = 0;
        }
    }
}

/*
 * Finalize absorbing phase with padding
 */
static void sha3_finalize(sha3_state_t* S) {
    size_t rate = S->rate;

    /* pad: domain || 10*1 */
    S->buf[S->buflen] = S->domain;
    memset(S->buf + S->buflen + 1, 0, rate - S->buflen - 1);
    S->buf[rate - 1] |= 0x80;

    /* absorb final block */
    for (size_t i = 0; i < rate / 8; i++) {
        S->state[i] ^= load64_le(S->buf + i * 8);
    }
    keccak_f1600(S->state);

    S->buflen = 0;
    S->finalized = 1;
}

/*
 * Squeeze output from sponge
 */
static void sha3_squeeze(sha3_state_t* S, uint8_t* out, size_t outlen) {
    size_t rate = S->rate;
    size_t offset = S->buflen;

    while (outlen > 0) {
        if (offset == rate) {
            keccak_f1600(S->state);
            offset = 0;
        }

        size_t chunk = rate - offset;
        if (chunk > outlen) chunk = outlen;

        /* extract from state */
        for (size_t i = 0; i < chunk; i++) {
            size_t lane = (offset + i) / 8;
            size_t byte_in_lane = (offset + i) % 8;
            out[i] = (uint8_t)(S->state[lane] >> (8 * byte_in_lane));
        }

        out += chunk;
        outlen -= chunk;
        offset += chunk;
    }

    S->buflen = offset;
}

/* ========================================================================== */
/* SHA3-224 */
/* ========================================================================== */

void cfx_sha3_224_init(cfx_sha3_ctx_t* ctx) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    S->rate = 144;  /* (1600 - 448) / 8 = 144 bytes */
    S->domain = 0x06;
}

void cfx_sha3_224_update(cfx_sha3_ctx_t* ctx, const void* data, size_t len) {
    sha3_absorb((sha3_state_t*)ctx, (const uint8_t*)data, len);
}

void cfx_sha3_224_final(cfx_sha3_ctx_t* ctx, uint8_t out[28]) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    sha3_finalize(S);
    sha3_squeeze(S, out, 28);
    memset(S, 0, sizeof(*S));
}

void cfx_sha3_224(uint8_t out[28], const void* data, size_t len) {
    cfx_sha3_ctx_t ctx;
    cfx_sha3_224_init(&ctx);
    cfx_sha3_224_update(&ctx, data, len);
    cfx_sha3_224_final(&ctx, out);
}

/* ========================================================================== */
/* SHA3-256 */
/* ========================================================================== */

void cfx_sha3_256_init(cfx_sha3_ctx_t* ctx) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    S->rate = 136;  /* (1600 - 512) / 8 = 136 bytes */
    S->domain = 0x06;
}

void cfx_sha3_256_update(cfx_sha3_ctx_t* ctx, const void* data, size_t len) {
    sha3_absorb((sha3_state_t*)ctx, (const uint8_t*)data, len);
}

void cfx_sha3_256_final(cfx_sha3_ctx_t* ctx, uint8_t out[32]) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    sha3_finalize(S);
    sha3_squeeze(S, out, 32);
    memset(S, 0, sizeof(*S));
}

void cfx_sha3_256(uint8_t out[32], const void* data, size_t len) {
    cfx_sha3_ctx_t ctx;
    cfx_sha3_256_init(&ctx);
    cfx_sha3_256_update(&ctx, data, len);
    cfx_sha3_256_final(&ctx, out);
}

/* ========================================================================== */
/* SHA3-384 */
/* ========================================================================== */

void cfx_sha3_384_init(cfx_sha3_ctx_t* ctx) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    S->rate = 104;  /* (1600 - 768) / 8 = 104 bytes */
    S->domain = 0x06;
}

void cfx_sha3_384_update(cfx_sha3_ctx_t* ctx, const void* data, size_t len) {
    sha3_absorb((sha3_state_t*)ctx, (const uint8_t*)data, len);
}

void cfx_sha3_384_final(cfx_sha3_ctx_t* ctx, uint8_t out[48]) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    sha3_finalize(S);
    sha3_squeeze(S, out, 48);
    memset(S, 0, sizeof(*S));
}

void cfx_sha3_384(uint8_t out[48], const void* data, size_t len) {
    cfx_sha3_ctx_t ctx;
    cfx_sha3_384_init(&ctx);
    cfx_sha3_384_update(&ctx, data, len);
    cfx_sha3_384_final(&ctx, out);
}

/* ========================================================================== */
/* SHA3-512 */
/* ========================================================================== */

void cfx_sha3_512_init(cfx_sha3_ctx_t* ctx) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    S->rate = 72;  /* (1600 - 1024) / 8 = 72 bytes */
    S->domain = 0x06;
}

void cfx_sha3_512_update(cfx_sha3_ctx_t* ctx, const void* data, size_t len) {
    sha3_absorb((sha3_state_t*)ctx, (const uint8_t*)data, len);
}

void cfx_sha3_512_final(cfx_sha3_ctx_t* ctx, uint8_t out[64]) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    sha3_finalize(S);
    sha3_squeeze(S, out, 64);
    memset(S, 0, sizeof(*S));
}

void cfx_sha3_512(uint8_t out[64], const void* data, size_t len) {
    cfx_sha3_ctx_t ctx;
    cfx_sha3_512_init(&ctx);
    cfx_sha3_512_update(&ctx, data, len);
    cfx_sha3_512_final(&ctx, out);
}

/* ========================================================================== */
/* SHAKE128 (XOF) */
/* ========================================================================== */

void cfx_shake128_init(cfx_sha3_ctx_t* ctx) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    S->rate = 168;  /* (1600 - 256) / 8 = 168 bytes */
    S->domain = 0x1F;
}

void cfx_shake128_update(cfx_sha3_ctx_t* ctx, const void* data, size_t len) {
    sha3_absorb((sha3_state_t*)ctx, (const uint8_t*)data, len);
}

void cfx_shake128_final(cfx_sha3_ctx_t* ctx) {
    sha3_finalize((sha3_state_t*)ctx);
}

void cfx_shake128_squeeze(cfx_sha3_ctx_t* ctx, void* out, size_t outlen) {
    sha3_squeeze((sha3_state_t*)ctx, (uint8_t*)out, outlen);
}

void cfx_shake128(void* out, size_t outlen, const void* data, size_t len) {
    cfx_sha3_ctx_t ctx;
    cfx_shake128_init(&ctx);
    cfx_shake128_update(&ctx, data, len);
    cfx_shake128_final(&ctx);
    cfx_shake128_squeeze(&ctx, out, outlen);
}

/* ========================================================================== */
/* SHAKE256 (XOF) */
/* ========================================================================== */

void cfx_shake256_init(cfx_sha3_ctx_t* ctx) {
    sha3_state_t* S = (sha3_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    S->rate = 136;  /* (1600 - 512) / 8 = 136 bytes */
    S->domain = 0x1F;
}

void cfx_shake256_update(cfx_sha3_ctx_t* ctx, const void* data, size_t len) {
    sha3_absorb((sha3_state_t*)ctx, (const uint8_t*)data, len);
}

void cfx_shake256_final(cfx_sha3_ctx_t* ctx) {
    sha3_finalize((sha3_state_t*)ctx);
}

void cfx_shake256_squeeze(cfx_sha3_ctx_t* ctx, void* out, size_t outlen) {
    sha3_squeeze((sha3_state_t*)ctx, (uint8_t*)out, outlen);
}

void cfx_shake256(void* out, size_t outlen, const void* data, size_t len) {
    cfx_sha3_ctx_t ctx;
    cfx_shake256_init(&ctx);
    cfx_shake256_update(&ctx, data, len);
    cfx_shake256_final(&ctx);
    cfx_shake256_squeeze(&ctx, out, outlen);
}
