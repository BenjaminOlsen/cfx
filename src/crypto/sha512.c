/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * sha512.c - SHA-512 hash function
 *
 * Implementation of SHA-512 as specified in FIPS 180-4.
 *
 * Self-guarding: compiles when no optimized backend is selected.
 */

#include "cfx/sha512.h"
#include "cfx/macros.h"
#include <string.h>

typedef struct {
    uint64_t state[8];
    uint64_t count[2];
    uint8_t buffer[128];
} sha512_state_t;

CFX_STATIC_ASSERT(
    sizeof(sha512_state_t) <= CFX_SHA512_CTX_SIZE,
    CFX_SHA512_CTX_SIZE_too_small_for_sha512_state
    );

#define CTX(p) ((sha512_state_t *)(p)->opaque)

static const uint64_t K[80] = {
    0x428a2f98d728ae22ULL, 0x7137449123ef65cdULL, 0xb5c0fbcfec4d3b2fULL, 0xe9b5dba58189dbbcULL,
    0x3956c25bf348b538ULL, 0x59f111f1b605d019ULL, 0x923f82a4af194f9bULL, 0xab1c5ed5da6d8118ULL,
    0xd807aa98a3030242ULL, 0x12835b0145706fbeULL, 0x243185be4ee4b28cULL, 0x550c7dc3d5ffb4e2ULL,
    0x72be5d74f27b896fULL, 0x80deb1fe3b1696b1ULL, 0x9bdc06a725c71235ULL, 0xc19bf174cf692694ULL,
    0xe49b69c19ef14ad2ULL, 0xefbe4786384f25e3ULL, 0x0fc19dc68b8cd5b5ULL, 0x240ca1cc77ac9c65ULL,
    0x2de92c6f592b0275ULL, 0x4a7484aa6ea6e483ULL, 0x5cb0a9dcbd41fbd4ULL, 0x76f988da831153b5ULL,
    0x983e5152ee66dfabULL, 0xa831c66d2db43210ULL, 0xb00327c898fb213fULL, 0xbf597fc7beef0ee4ULL,
    0xc6e00bf33da88fc2ULL, 0xd5a79147930aa725ULL, 0x06ca6351e003826fULL, 0x142929670a0e6e70ULL,
    0x27b70a8546d22ffcULL, 0x2e1b21385c26c926ULL, 0x4d2c6dfc5ac42aedULL, 0x53380d139d95b3dfULL,
    0x650a73548baf63deULL, 0x766a0abb3c77b2a8ULL, 0x81c2c92e47edaee6ULL, 0x92722c851482353bULL,
    0xa2bfe8a14cf10364ULL, 0xa81a664bbc423001ULL, 0xc24b8b70d0f89791ULL, 0xc76c51a30654be30ULL,
    0xd192e819d6ef5218ULL, 0xd69906245565a910ULL, 0xf40e35855771202aULL, 0x106aa07032bbd1b8ULL,
    0x19a4c116b8d2d0c8ULL, 0x1e376c085141ab53ULL, 0x2748774cdf8eeb99ULL, 0x34b0bcb5e19b48a8ULL,
    0x391c0cb3c5c95a63ULL, 0x4ed8aa4ae3418acbULL, 0x5b9cca4f7763e373ULL, 0x682e6ff3d6b2b8a3ULL,
    0x748f82ee5defb2fcULL, 0x78a5636f43172f60ULL, 0x84c87814a1f0ab72ULL, 0x8cc702081a6439ecULL,
    0x90befffa23631e28ULL, 0xa4506cebde82bde9ULL, 0xbef9a3f7b2c67915ULL, 0xc67178f2e372532bULL,
    0xca273eceea26619cULL, 0xd186b8c721c0c207ULL, 0xeada7dd6cde0eb1eULL, 0xf57d4f7fee6ed178ULL,
    0x06f067aa72176fbaULL, 0x0a637dc5a2c898a6ULL, 0x113f9804bef90daeULL, 0x1b710b35131c471bULL,
    0x28db77f523047d84ULL, 0x32caab7b40c72493ULL, 0x3c9ebe0a15c9bebcULL, 0x431d67c49c100d4cULL,
    0x4cc5d4becb3e42b6ULL, 0x597f299cfc657e2aULL, 0x5fcb6fab3ad6faecULL, 0x6c44198c4a475817ULL
};

#define ROTR64(x, n) (((x) >> (n)) | ((x) << (64 - (n))))

#define CH(x, y, z)  (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x, y, z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x)       (ROTR64(x, 28) ^ ROTR64(x, 34) ^ ROTR64(x, 39))
#define EP1(x)       (ROTR64(x, 14) ^ ROTR64(x, 18) ^ ROTR64(x, 41))
#define SIG0(x)      (ROTR64(x, 1) ^ ROTR64(x, 8) ^ ((x) >> 7))
#define SIG1(x)      (ROTR64(x, 19) ^ ROTR64(x, 61) ^ ((x) >> 6))

static void sha512_transform(sha512_state_t *st, const uint8_t *block) {
    uint64_t a, b, c, d, e, f, g, h, t1, t2, W[80];
    int i;

    for (i = 0; i < 16; i++) {
        W[i] = ((uint64_t)block[i * 8 + 0] << 56) |
               ((uint64_t)block[i * 8 + 1] << 48) |
               ((uint64_t)block[i * 8 + 2] << 40) |
               ((uint64_t)block[i * 8 + 3] << 32) |
               ((uint64_t)block[i * 8 + 4] << 24) |
               ((uint64_t)block[i * 8 + 5] << 16) |
               ((uint64_t)block[i * 8 + 6] << 8) |
               ((uint64_t)block[i * 8 + 7]);
    }
    for (i = 16; i < 80; i++) {
        W[i] = SIG1(W[i - 2]) + W[i - 7] + SIG0(W[i - 15]) + W[i - 16];
    }

    a = st->state[0];
    b = st->state[1];
    c = st->state[2];
    d = st->state[3];
    e = st->state[4];
    f = st->state[5];
    g = st->state[6];
    h = st->state[7];

    for (i = 0; i < 80; i++) {
        t1 = h + EP1(e) + CH(e, f, g) + K[i] + W[i];
        t2 = EP0(a) + MAJ(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + t1;
        d = c;
        c = b;
        b = a;
        a = t1 + t2;
    }

    st->state[0] += a;
    st->state[1] += b;
    st->state[2] += c;
    st->state[3] += d;
    st->state[4] += e;
    st->state[5] += f;
    st->state[6] += g;
    st->state[7] += h;
}

void cfx_sha512_init(cfx_sha512_ctx_t *ctx) {
    sha512_state_t *st = CTX(ctx);
    st->state[0] = 0x6a09e667f3bcc908ULL;
    st->state[1] = 0xbb67ae8584caa73bULL;
    st->state[2] = 0x3c6ef372fe94f82bULL;
    st->state[3] = 0xa54ff53a5f1d36f1ULL;
    st->state[4] = 0x510e527fade682d1ULL;
    st->state[5] = 0x9b05688c2b3e6c1fULL;
    st->state[6] = 0x1f83d9abfb41bd6bULL;
    st->state[7] = 0x5be0cd19137e2179ULL;
    st->count[0] = 0;
    st->count[1] = 0;
}

void cfx_sha512_update(cfx_sha512_ctx_t *ctx, const uint8_t *data, size_t len) {
    sha512_state_t *st = CTX(ctx);
    size_t i, index, partLen;

    index = (size_t)((st->count[0] >> 3) & 0x7f);
    st->count[0] += ((uint64_t)len) << 3;
    if (st->count[0] < (((uint64_t)len) << 3)) {
        st->count[1]++;
    }
    st->count[1] += ((uint64_t)len) >> 61;

    partLen = 128 - index;

    if (len >= partLen) {
        memcpy(&st->buffer[index], data, partLen);
        sha512_transform(st, st->buffer);
        for (i = partLen; i + 127 < len; i += 128) {
            sha512_transform(st, &data[i]);
        }
        index = 0;
    } else {
        i = 0;
    }

    memcpy(&st->buffer[index], &data[i], len - i);
}

void cfx_sha512_final(cfx_sha512_ctx_t *ctx, uint8_t out[64]) {
    sha512_state_t *st = CTX(ctx);
    uint8_t bits[16];
    size_t index, padLen;
    static const uint8_t padding[128] = { 0x80 };
    int i;

    for (i = 0; i < 8; i++) {
        bits[i] = (uint8_t)(st->count[1] >> (56 - i * 8));
        bits[i + 8] = (uint8_t)(st->count[0] >> (56 - i * 8));
    }

    index = (size_t)((st->count[0] >> 3) & 0x7f);
    padLen = (index < 112) ? (112 - index) : (240 - index);
    cfx_sha512_update(ctx, padding, padLen);
    cfx_sha512_update(ctx, bits, 16);

    for (i = 0; i < 8; i++) {
        out[i * 8 + 0] = (uint8_t)(st->state[i] >> 56);
        out[i * 8 + 1] = (uint8_t)(st->state[i] >> 48);
        out[i * 8 + 2] = (uint8_t)(st->state[i] >> 40);
        out[i * 8 + 3] = (uint8_t)(st->state[i] >> 32);
        out[i * 8 + 4] = (uint8_t)(st->state[i] >> 24);
        out[i * 8 + 5] = (uint8_t)(st->state[i] >> 16);
        out[i * 8 + 6] = (uint8_t)(st->state[i] >> 8);
        out[i * 8 + 7] = (uint8_t)(st->state[i]);
    }
}

void cfx_sha512(uint8_t out[64], const uint8_t *data, size_t len) {
    cfx_sha512_ctx_t ctx;
    cfx_sha512_init(&ctx);
    cfx_sha512_update(&ctx, data, len);
    cfx_sha512_final(&ctx, out);
}
