#include "cfx/sha256.h"
#include "cfx/macros.h"
#include "cfx/bits.h"

#include <string.h>

typedef struct {
    uint32_t state[8];     /* H0..H7 */
    uint64_t total_bits;   /* total message length in bits (mod 2^64) */
    uint8_t buffer[64];    /* partial block buffer */
    size_t buffer_len;     /* number of bytes currently in buffer [0..64) */
} cfx_sha256_state_t;

CFX_STATIC_ASSERT(
    sizeof(cfx_sha256_state_t) <= CFX_SHA256_CTX_SIZE,
    CFX_SHA256_CTX_SIZE_too_small_for_cfx_sha256_state
    );

#define rotr32 cfx_rotr32

static inline uint32_t Sigma0(uint32_t x) {
    return rotr32(x, 2) ^ rotr32(x, 13) ^ rotr32(x, 22);
}

static inline uint32_t Sigma1(uint32_t x) {
    return rotr32(x, 6) ^ rotr32(x, 11) ^ rotr32(x, 25);
}

static inline uint32_t sigma0(uint32_t x) {
    return rotr32(x, 7) ^ rotr32(x, 18) ^ (x >> 3);
}

static inline uint32_t sigma1(uint32_t x) {
    return rotr32(x, 17) ^ rotr32(x, 19) ^ (x >> 10);
}

static inline uint32_t Ch(uint32_t x, uint32_t y, uint32_t z) {
    return (x & y) ^ (~x & z);
}

static inline uint32_t Maj(uint32_t x, uint32_t y, uint32_t z) {
    return (x & y) ^ (x & z) ^ (y & z);
}

static inline uint32_t load_be32(const uint8_t *p) {
    return ((uint32_t)p[0] << 24)
           | ((uint32_t)p[1] << 16)
           | ((uint32_t)p[2] <<  8)
           | ((uint32_t)p[3]);
}

static inline void store_be32(uint8_t *p, uint32_t x) {
    p[0] = (uint8_t)(x >> 24);
    p[1] = (uint8_t)(x >> 16);
    p[2] = (uint8_t)(x >> 8);
    p[3] = (uint8_t)(x);
}

static const uint32_t K[64] = {
    0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u,
    0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
    0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
    0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
    0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu,
    0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
    0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u,
    0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
    0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
    0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
    0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u,
    0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
    0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u,
    0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
    0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
    0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u
};

static void cfx_sha256_compress(cfx_sha256_state_t *st, const uint8_t block[64]) {
    uint32_t W[64];
    uint32_t a, b, c, d, e, f, g, h;

    for (unsigned i = 0; i < 16; ++i) {
        W[i] = load_be32(block + 4u * i);
    }

    for (unsigned t = 16; t < 64; ++t) {
        uint32_t s0 = sigma0(W[t - 15]);
        uint32_t s1 = sigma1(W[t - 2]);
        W[t] = W[t - 16] + s0 + W[t - 7] + s1;
    }

    a = st->state[0];
    b = st->state[1];
    c = st->state[2];
    d = st->state[3];
    e = st->state[4];
    f = st->state[5];
    g = st->state[6];
    h = st->state[7];

    for (unsigned t = 0; t < 64; ++t) {
        uint32_t T1 = h + Sigma1(e) + Ch(e, f, g) + K[t] + W[t];
        uint32_t T2 = Sigma0(a) + Maj(a, b, c);

        h = g;
        g = f;
        f = e;
        e = d + T1;
        d = c;
        c = b;
        b = a;
        a = T1 + T2;
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



void cfx_sha256_init(cfx_sha256_ctx *ctx) {
    cfx_sha256_state_t *st = (cfx_sha256_state_t *)ctx->opaque;

    st->state[0] = 0x6a09e667u;
    st->state[1] = 0xbb67ae85u;
    st->state[2] = 0x3c6ef372u;
    st->state[3] = 0xa54ff53au;
    st->state[4] = 0x510e527fu;
    st->state[5] = 0x9b05688cu;
    st->state[6] = 0x1f83d9abu;
    st->state[7] = 0x5be0cd19u;

    st->total_bits = 0;
    st->buffer_len = 0;
}

void cfx_sha256_update(cfx_sha256_ctx *ctx, const uint8_t *data, size_t len) {
    cfx_sha256_state_t *st = (cfx_sha256_state_t *)ctx->opaque;

    if (len == 0) {
        return;
    }

    st->total_bits += (uint64_t)len * 8u;

    size_t buffer_len = st->buffer_len;

    if (buffer_len > 0) {
        size_t to_fill = 64u - buffer_len;
        if (len < to_fill) {
            memcpy(st->buffer + buffer_len, data, len);
            st->buffer_len = buffer_len + len;
            return;
        }

        memcpy(st->buffer + buffer_len, data, to_fill);
        cfx_sha256_compress(st, st->buffer);

        data        += to_fill;
        len         -= to_fill;
        buffer_len   = 0;
    }

    while (len >= 64u) {
        cfx_sha256_compress(st, data);
        data += 64u;
        len  -= 64u;
    }

    if (len > 0) {
        memcpy(st->buffer + buffer_len, data, len);
        st->buffer_len = buffer_len + len;
    } else {
        st->buffer_len = 0;
    }
}

/* add padding, length, and output hash */
void cfx_sha256_final(cfx_sha256_ctx *ctx, uint8_t out[32]) {
    cfx_sha256_state_t *st = (cfx_sha256_state_t *)ctx->opaque;
    uint8_t *buf = st->buffer;
    size_t i   = st->buffer_len;

    /* append the '1' bit as 0x80 */
    buf[i++] = 0x80;

    /* If there isn't enough room for length (8 bytes) in this block,
       pad with zeros, compress, and start a new block. */
    if (i > 56u) {
        while (i < 64u) {
            buf[i++] = 0x00;
        }
        cfx_sha256_compress(st, buf);
        i = 0;
    }


    while (i < 56u) {
        buf[i++] = 0x00;
    }

    uint64_t bits = st->total_bits;
    buf[56] = (uint8_t)(bits >> 56);
    buf[57] = (uint8_t)(bits >> 48);
    buf[58] = (uint8_t)(bits >> 40);
    buf[59] = (uint8_t)(bits >> 32);
    buf[60] = (uint8_t)(bits >> 24);
    buf[61] = (uint8_t)(bits >> 16);
    buf[62] = (uint8_t)(bits >>  8);
    buf[63] = (uint8_t)(bits      );

    cfx_sha256_compress(st, buf);

    /* "digest" */
    for (unsigned j = 0; j < 8; ++j) {
        store_be32(out + 4u * j, st->state[j]);
    }

    memset(st, 0, sizeof(*st));
}

void cfx_sha256(uint8_t out[32], const uint8_t *data, size_t len) {
    cfx_sha256_ctx ctx;
    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, data, len);
    cfx_sha256_final(&ctx, out);
}