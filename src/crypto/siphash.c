#include "cfx/siphash.h"
#include "cfx/macros.h"
#include "cfx/memory.h"
#include "cfx/bits.h"
#include <string.h>

typedef struct {
    uint64_t v[4];
    uint64_t total;
    uint8_t buf[8];
    uint8_t buflen;
    uint8_t c_rounds;
    uint8_t d_rounds;
} cfx_siphash_state_t;

CFX_STATIC_ASSERT(
    sizeof(cfx_siphash_state_t) <= CFX_SIPHASH_CTX_SIZE,
    CFX_SIPHASH_CTX_SIZE_too_small
    );

#define rotl64 cfx_rotl64

#define SIPROUND(v0, v1, v2, v3) do { \
            v0 += v1; v1 = rotl64(v1, 13); v1 ^= v0; v0 = rotl64(v0, 32); \
            v2 += v3; v3 = rotl64(v3, 16); v3 ^= v2; \
            v0 += v3; v3 = rotl64(v3, 21); v3 ^= v0; \
            v2 += v1; v1 = rotl64(v1, 17); v1 ^= v2; v2 = rotl64(v2, 32); \
} while (0)

uint64_t cfx_siphash_cd(const uint8_t *data, size_t len, const uint8_t key[16], int c_rounds, int d_rounds) {
    if (c_rounds < 1 || d_rounds < 1) return 0;

    uint64_t k0 = cfx_load64_le(key);
    uint64_t k1 = cfx_load64_le(key + 8);

    uint64_t v0 = k0 ^ 0x736f6d6570736575ULL;
    uint64_t v1 = k1 ^ 0x646f72616e646f6dULL;
    uint64_t v2 = k0 ^ 0x6c7967656e657261ULL;
    uint64_t v3 = k1 ^ 0x7465646279746573ULL;

    const uint8_t *end = data + (len & ~7ULL);
    const size_t left = len & 7;

    /* process full 8-byte blocks */
    while (data < end) {
        uint64_t m = cfx_load64_le(data);
        v3 ^= m;
        for (int i = 0; i < c_rounds; i++) SIPROUND(v0, v1, v2, v3);
        v0 ^= m;
        data += 8;
    }

    /* final block with length byte */
    uint64_t b = (uint64_t)len << 56;
    switch (left) {
    case 7: b |= (uint64_t)data[6] << 48;     /* fall through */
    case 6: b |= (uint64_t)data[5] << 40;     /* fall through */
    case 5: b |= (uint64_t)data[4] << 32;     /* fall through */
    case 4: b |= (uint64_t)data[3] << 24;     /* fall through */
    case 3: b |= (uint64_t)data[2] << 16;     /* fall through */
    case 2: b |= (uint64_t)data[1] << 8;      /* fall through */
    case 1: b |= (uint64_t)data[0];           /* fall through */
    case 0: break;
    }

    v3 ^= b;
    for (int i = 0; i < c_rounds; i++) SIPROUND(v0, v1, v2, v3);
    v0 ^= b;

    v2 ^= 0xff;
    for (int i = 0; i < d_rounds; i++) SIPROUND(v0, v1, v2, v3);

    return v0 ^ v1 ^ v2 ^ v3;
}

uint64_t cfx_siphash(const uint8_t *data, size_t len, const uint8_t key[16]) {
    return cfx_siphash_cd(data, len, key, 2, 4);
}

void cfx_siphash128(uint8_t out[16], const uint8_t *data, size_t len,
    const uint8_t key[16]) {
    uint64_t k0 = cfx_load64_le(key);
    uint64_t k1 = cfx_load64_le(key + 8);

    uint64_t v0 = k0 ^ 0x736f6d6570736575ULL;
    uint64_t v1 = k1 ^ 0x646f72616e646f6dULL;
    uint64_t v2 = k0 ^ 0x6c7967656e657261ULL;
    uint64_t v3 = k1 ^ 0x7465646279746573ULL;

    v1 ^= 0xee;  /* 128-bit output mode */

    const uint8_t *end = data + (len & ~7ULL);
    const size_t left = len & 7;

    while (data < end) {
        uint64_t m = cfx_load64_le(data);
        v3 ^= m;
        SIPROUND(v0, v1, v2, v3);
        SIPROUND(v0, v1, v2, v3);
        v0 ^= m;
        data += 8;
    }

    uint64_t b = (uint64_t)len << 56;
    switch (left) {
    case 7: b |= (uint64_t)data[6] << 48;     /* fall through */
    case 6: b |= (uint64_t)data[5] << 40;     /* fall through */
    case 5: b |= (uint64_t)data[4] << 32;     /* fall through */
    case 4: b |= (uint64_t)data[3] << 24;     /* fall through */
    case 3: b |= (uint64_t)data[2] << 16;     /* fall through */
    case 2: b |= (uint64_t)data[1] << 8;      /* fall through */
    case 1: b |= (uint64_t)data[0];           /* fall through */
    case 0: break;
    }

    v3 ^= b;
    SIPROUND(v0, v1, v2, v3);
    SIPROUND(v0, v1, v2, v3);
    v0 ^= b;

    v2 ^= 0xee;
    SIPROUND(v0, v1, v2, v3);
    SIPROUND(v0, v1, v2, v3);
    SIPROUND(v0, v1, v2, v3);
    SIPROUND(v0, v1, v2, v3);

    uint64_t h0 = v0 ^ v1 ^ v2 ^ v3;

    v1 ^= 0xdd;
    SIPROUND(v0, v1, v2, v3);
    SIPROUND(v0, v1, v2, v3);
    SIPROUND(v0, v1, v2, v3);
    SIPROUND(v0, v1, v2, v3);

    uint64_t h1 = v0 ^ v1 ^ v2 ^ v3;

    /* store little-endian */
    for (int i = 0; i < 8; i++) out[i]     = (uint8_t)(h0 >> (8 * i));
    for (int i = 0; i < 8; i++) out[8 + i] = (uint8_t)(h1 >> (8 * i));
}

/* streaming interface */

void cfx_siphash_init_cd(cfx_siphash_ctx_t *ctx, const uint8_t key[16], int c_rounds, int d_rounds) {
    if (c_rounds < 1 || d_rounds < 1) return;

    cfx_siphash_state_t *s = (cfx_siphash_state_t *)ctx;
    uint64_t k0 = cfx_load64_le(key);
    uint64_t k1 = cfx_load64_le(key + 8);

    s->v[0] = k0 ^ 0x736f6d6570736575ULL;
    s->v[1] = k1 ^ 0x646f72616e646f6dULL;
    s->v[2] = k0 ^ 0x6c7967656e657261ULL;
    s->v[3] = k1 ^ 0x7465646279746573ULL;
    s->total = 0;
    s->buflen = 0;
    s->c_rounds = (uint8_t)c_rounds;
    s->d_rounds = (uint8_t)d_rounds;
}

void cfx_siphash_init(cfx_siphash_ctx_t *ctx, const uint8_t key[16]) {
    cfx_siphash_init_cd(ctx, key, 2, 4);
}

void cfx_siphash_update(cfx_siphash_ctx_t *ctx, const uint8_t *data, size_t len) {
    cfx_siphash_state_t *s = (cfx_siphash_state_t *)ctx;
    s->total += len;

    /* fill buffer if we have partial data */
    if (s->buflen > 0) {
        size_t need = 8 - s->buflen;
        if (len < need) {
            memcpy(s->buf + s->buflen, data, len);
            s->buflen += (uint8_t)len;
            return;
        }
        memcpy(s->buf + s->buflen, data, need);
        uint64_t m = cfx_load64_le(s->buf);
        s->v[3] ^= m;
        for (int i = 0; i < s->c_rounds; i++)
            SIPROUND(s->v[0], s->v[1], s->v[2], s->v[3]);
        s->v[0] ^= m;
        data += need;
        len -= need;
        s->buflen = 0;
    }

    /* process full blocks */
    while (len >= 8) {
        uint64_t m = cfx_load64_le(data);
        s->v[3] ^= m;
        for (int i = 0; i < s->c_rounds; i++)
            SIPROUND(s->v[0], s->v[1], s->v[2], s->v[3]);
        s->v[0] ^= m;
        data += 8;
        len -= 8;
    }

    /* save remaining */
    if (len > 0) {
        memcpy(s->buf, data, len);
        s->buflen = (uint8_t)len;
    }
}

uint64_t cfx_siphash_final(cfx_siphash_ctx_t *ctx) {
    cfx_siphash_state_t *s = (cfx_siphash_state_t *)ctx;
    uint64_t b = s->total << 56;

    switch (s->buflen) {
    case 7: b |= (uint64_t)s->buf[6] << 48;     /* fall through */
    case 6: b |= (uint64_t)s->buf[5] << 40;     /* fall through */
    case 5: b |= (uint64_t)s->buf[4] << 32;     /* fall through */
    case 4: b |= (uint64_t)s->buf[3] << 24;     /* fall through */
    case 3: b |= (uint64_t)s->buf[2] << 16;     /* fall through */
    case 2: b |= (uint64_t)s->buf[1] << 8;      /* fall through */
    case 1: b |= (uint64_t)s->buf[0];           /* fall through */
    case 0: break;
    }

    s->v[3] ^= b;
    for (int i = 0; i < s->c_rounds; i++)
        SIPROUND(s->v[0], s->v[1], s->v[2], s->v[3]);
    s->v[0] ^= b;

    s->v[2] ^= 0xff;
    for (int i = 0; i < s->d_rounds; i++)
        SIPROUND(s->v[0], s->v[1], s->v[2], s->v[3]);

    return s->v[0] ^ s->v[1] ^ s->v[2] ^ s->v[3];
}
