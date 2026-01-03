/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * blake2.c - BLAKE2b and BLAKE2s implementation (RFC 7693)
 */

#include "cfx/blake2.h"
#include "cfx/macros.h"
#include "cfx/arch.h"
#include <string.h>

#if CFX_HAVE_AVX2
#include <immintrin.h>
#endif

/* internal state types - must fit in opaque ctx */
typedef struct {
    uint64_t h[8];
    uint64_t t[2];
    uint64_t f[2];
    uint8_t  buf[128];
    size_t   buflen;
    size_t   outlen;
} cfx_blake2b_state_t;

typedef struct {
    uint32_t h[8];
    uint32_t t[2];
    uint32_t f[2];
    uint8_t  buf[64];
    size_t   buflen;
    size_t   outlen;
} cfx_blake2s_state_t;

CFX_STATIC_ASSERT(sizeof(cfx_blake2b_state_t) <= sizeof(cfx_blake2b_ctx_t), blake2b_ctx_too_small);
CFX_STATIC_ASSERT(sizeof(cfx_blake2s_state_t) <= sizeof(cfx_blake2s_ctx_t), blake2s_ctx_too_small);

/* BLAKE2b IV (same as SHA-512 IV) */
static const uint64_t blake2b_IV[8] = {
    0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL,
    0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
    0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
    0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL
};

/* BLAKE2s IV (same as SHA-256 IV) */
static const uint32_t blake2s_IV[8] = {
    0x6a09e667UL, 0xbb67ae85UL, 0x3c6ef372UL, 0xa54ff53aUL,
    0x510e527fUL, 0x9b05688cUL, 0x1f83d9abUL, 0x5be0cd19UL
};

/* sigma permutations for message schedule */
static const uint8_t sigma[12][16] = {
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 },
    { 11,  8, 12,  0,  5,  2, 15, 13, 10, 14,  3,  6,  7,  1,  9,  4 },
    {  7,  9,  3,  1, 13, 12, 11, 14,  2,  6,  5, 10,  4,  0, 15,  8 },
    {  9,  0,  5,  7,  2,  4, 10, 15, 14,  1, 11, 12,  6,  8,  3, 13 },
    {  2, 12,  6, 10,  0, 11,  8,  3,  4, 13,  7,  5, 15, 14,  1,  9 },
    { 12,  5,  1, 15, 14, 13,  4, 10,  0,  7,  6,  3,  9,  2,  8, 11 },
    { 13, 11,  7, 14, 12,  1,  3,  9,  5,  0, 15,  4,  8,  6,  2, 10 },
    {  6, 15, 14,  9, 11,  3,  0,  8, 12,  2, 13,  7,  1,  4, 10,  5 },
    { 10,  2,  8,  4,  7,  6,  1,  5, 15, 11,  9, 14,  3, 12, 13,  0 },
    {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15 },
    { 14, 10,  4,  8,  9, 15, 13,  6,  1, 12,  0,  2, 11,  7,  5,  3 }
};

static inline uint64_t rotr64(uint64_t x, int n) {
    return (x >> n) | (x << (64 - n));
}

static inline uint32_t rotr32(uint32_t x, int n) {
    return (x >> n) | (x << (32 - n));
}

static inline uint64_t load64_le(const void* src) {
    const uint8_t* p = (const uint8_t*)src;
    return ((uint64_t)p[0])       | ((uint64_t)p[1] << 8)  |
           ((uint64_t)p[2] << 16) | ((uint64_t)p[3] << 24) |
           ((uint64_t)p[4] << 32) | ((uint64_t)p[5] << 40) |
           ((uint64_t)p[6] << 48) | ((uint64_t)p[7] << 56);
}

static inline uint32_t load32_le(const void* src) {
    const uint8_t* p = (const uint8_t*)src;
    return ((uint32_t)p[0])       | ((uint32_t)p[1] << 8) |
           ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}

static inline void store64_le(void* dst, uint64_t x) {
    uint8_t* p = (uint8_t*)dst;
    p[0] = (uint8_t)x; p[1] = (uint8_t)(x >> 8);
    p[2] = (uint8_t)(x >> 16); p[3] = (uint8_t)(x >> 24);
    p[4] = (uint8_t)(x >> 32); p[5] = (uint8_t)(x >> 40);
    p[6] = (uint8_t)(x >> 48); p[7] = (uint8_t)(x >> 56);
}

static inline void store32_le(void* dst, uint32_t x) {
    uint8_t* p = (uint8_t*)dst;
    p[0] = (uint8_t)x; p[1] = (uint8_t)(x >> 8);
    p[2] = (uint8_t)(x >> 16); p[3] = (uint8_t)(x >> 24);
}

/* ========================================================================== */
/* BLAKE2b */
/* ========================================================================== */

#define B2B_G(a, b, c, d, x, y) do { \
    a = a + b + x; \
    d = rotr64(d ^ a, 32); \
    c = c + d; \
    b = rotr64(b ^ c, 24); \
    a = a + b + y; \
    d = rotr64(d ^ a, 16); \
    c = c + d; \
    b = rotr64(b ^ c, 63); \
} while(0)

CFX_MAYBE_UNUSED
static void blake2b_compress(cfx_blake2b_state_t* S, const uint8_t block[128]) {
    uint64_t m[16], v[16];

    for (int i = 0; i < 16; i++) {
        m[i] = load64_le(block + i * 8);
    }

    for (int i = 0; i < 8; i++) {
        v[i] = S->h[i];
        v[i + 8] = blake2b_IV[i];
    }

    v[12] ^= S->t[0];
    v[13] ^= S->t[1];
    v[14] ^= S->f[0];
    v[15] ^= S->f[1];

    for (int r = 0; r < 12; r++) {
        const uint8_t* s = sigma[r];
        B2B_G(v[0], v[4], v[ 8], v[12], m[s[ 0]], m[s[ 1]]);
        B2B_G(v[1], v[5], v[ 9], v[13], m[s[ 2]], m[s[ 3]]);
        B2B_G(v[2], v[6], v[10], v[14], m[s[ 4]], m[s[ 5]]);
        B2B_G(v[3], v[7], v[11], v[15], m[s[ 6]], m[s[ 7]]);
        B2B_G(v[0], v[5], v[10], v[15], m[s[ 8]], m[s[ 9]]);
        B2B_G(v[1], v[6], v[11], v[12], m[s[10]], m[s[11]]);
        B2B_G(v[2], v[7], v[ 8], v[13], m[s[12]], m[s[13]]);
        B2B_G(v[3], v[4], v[ 9], v[14], m[s[14]], m[s[15]]);
    }

    for (int i = 0; i < 8; i++) {
        S->h[i] ^= v[i] ^ v[i + 8];
    }
}

static void blake2b_increment_counter(cfx_blake2b_state_t* S, uint64_t inc) {
    S->t[0] += inc;
    S->t[1] += (S->t[0] < inc) ? 1 : 0;
}

#if CFX_HAVE_AVX2

static inline __m256i rotr64_32(__m256i x) {
    return _mm256_shuffle_epi32(x, _MM_SHUFFLE(2, 3, 0, 1));
}

static inline __m256i rotr64_24(__m256i x) {
    const __m256i mask = _mm256_setr_epi8(
        3, 4, 5, 6, 7, 0, 1, 2, 11, 12, 13, 14, 15, 8, 9, 10,
        3, 4, 5, 6, 7, 0, 1, 2, 11, 12, 13, 14, 15, 8, 9, 10
    );
    return _mm256_shuffle_epi8(x, mask);
}

static inline __m256i rotr64_16(__m256i x) {
    const __m256i mask = _mm256_setr_epi8(
        2, 3, 4, 5, 6, 7, 0, 1, 10, 11, 12, 13, 14, 15, 8, 9,
        2, 3, 4, 5, 6, 7, 0, 1, 10, 11, 12, 13, 14, 15, 8, 9
    );
    return _mm256_shuffle_epi8(x, mask);
}

static inline __m256i rotr64_63(__m256i x) {
    return _mm256_or_si256(_mm256_srli_epi64(x, 63), _mm256_slli_epi64(x, 1));
}

#define B2B_G_AVX2(row0, row1, row2, row3, m_x, m_y) do { \
    row0 = _mm256_add_epi64(_mm256_add_epi64(row0, row1), m_x); \
    row3 = rotr64_32(_mm256_xor_si256(row3, row0)); \
    row2 = _mm256_add_epi64(row2, row3); \
    row1 = rotr64_24(_mm256_xor_si256(row1, row2)); \
    row0 = _mm256_add_epi64(_mm256_add_epi64(row0, row1), m_y); \
    row3 = rotr64_16(_mm256_xor_si256(row3, row0)); \
    row2 = _mm256_add_epi64(row2, row3); \
    row1 = rotr64_63(_mm256_xor_si256(row1, row2)); \
} while(0)

/* row1 rotates by 1, row2 by 2, row3 by 3 for diagonal round */
#define DIAG_PERM(row1, row2, row3) do { \
    row1 = _mm256_permute4x64_epi64(row1, _MM_SHUFFLE(0, 3, 2, 1)); \
    row2 = _mm256_permute4x64_epi64(row2, _MM_SHUFFLE(1, 0, 3, 2)); \
    row3 = _mm256_permute4x64_epi64(row3, _MM_SHUFFLE(2, 1, 0, 3)); \
} while(0)

#define UNDIAG_PERM(row1, row2, row3) do { \
    row1 = _mm256_permute4x64_epi64(row1, _MM_SHUFFLE(2, 1, 0, 3)); \
    row2 = _mm256_permute4x64_epi64(row2, _MM_SHUFFLE(1, 0, 3, 2)); \
    row3 = _mm256_permute4x64_epi64(row3, _MM_SHUFFLE(0, 3, 2, 1)); \
} while(0)

static void blake2b_compress_avx2(cfx_blake2b_state_t* S, const uint8_t block[128]) {
    __m256i row0, row1, row2, row3;

    row0 = _mm256_loadu_si256((const __m256i*)&S->h[0]);
    row1 = _mm256_loadu_si256((const __m256i*)&S->h[4]);
    row2 = _mm256_loadu_si256((const __m256i*)&blake2b_IV[0]);
    row3 = _mm256_loadu_si256((const __m256i*)&blake2b_IV[4]);

    row3 = _mm256_xor_si256(row3, _mm256_setr_epi64x(
        (int64_t)S->t[0], (int64_t)S->t[1], (int64_t)S->f[0], (int64_t)S->f[1]
    ));

    for (int r = 0; r < 12; r++) {
        const uint8_t* s = sigma[r];

        /* column round */
        __m256i mx0 = _mm256_setr_epi64x(
            (int64_t)load64_le(block + s[0]*8), (int64_t)load64_le(block + s[2]*8),
            (int64_t)load64_le(block + s[4]*8), (int64_t)load64_le(block + s[6]*8));
        __m256i my0 = _mm256_setr_epi64x(
            (int64_t)load64_le(block + s[1]*8), (int64_t)load64_le(block + s[3]*8),
            (int64_t)load64_le(block + s[5]*8), (int64_t)load64_le(block + s[7]*8));
        B2B_G_AVX2(row0, row1, row2, row3, mx0, my0);

        DIAG_PERM(row1, row2, row3);

        /* diagonal round */
        __m256i mx1 = _mm256_setr_epi64x(
            (int64_t)load64_le(block + s[8]*8), (int64_t)load64_le(block + s[10]*8),
            (int64_t)load64_le(block + s[12]*8), (int64_t)load64_le(block + s[14]*8));
        __m256i my1 = _mm256_setr_epi64x(
            (int64_t)load64_le(block + s[9]*8), (int64_t)load64_le(block + s[11]*8),
            (int64_t)load64_le(block + s[13]*8), (int64_t)load64_le(block + s[15]*8));
        B2B_G_AVX2(row0, row1, row2, row3, mx1, my1);

        UNDIAG_PERM(row1, row2, row3);
    }

    __m256i h0 = _mm256_loadu_si256((const __m256i*)&S->h[0]);
    __m256i h1 = _mm256_loadu_si256((const __m256i*)&S->h[4]);
    h0 = _mm256_xor_si256(h0, _mm256_xor_si256(row0, row2));
    h1 = _mm256_xor_si256(h1, _mm256_xor_si256(row1, row3));
    _mm256_storeu_si256((__m256i*)&S->h[0], h0);
    _mm256_storeu_si256((__m256i*)&S->h[4], h1);
}

#undef B2B_G_AVX2
#undef DIAG_PERM
#undef UNDIAG_PERM

#define blake2b_compress_fn blake2b_compress_avx2
#else
#define blake2b_compress_fn blake2b_compress
#endif /* CFX_HAVE_AVX2 */

int cfx_blake2b_init(cfx_blake2b_ctx_t* ctx, size_t outlen) {
    if (outlen == 0 || outlen > CFX_BLAKE2B_OUTBYTES) return -1;

    cfx_blake2b_state_t* S = (cfx_blake2b_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    for (int i = 0; i < 8; i++) {
        S->h[i] = blake2b_IV[i];
    }

    S->h[0] ^= 0x01010000 ^ outlen;
    S->outlen = outlen;

    return 0;
}

int cfx_blake2b_init_key(cfx_blake2b_ctx_t* ctx, size_t outlen,
                         const void* key, size_t keylen) {
    if (outlen == 0 || outlen > CFX_BLAKE2B_OUTBYTES) return -1;
    if (keylen > CFX_BLAKE2B_KEYBYTES) return -1;

    cfx_blake2b_state_t* S = (cfx_blake2b_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    for (int i = 0; i < 8; i++) {
        S->h[i] = blake2b_IV[i];
    }

    S->h[0] ^= 0x01010000 ^ (keylen << 8) ^ outlen;
    S->outlen = outlen;

    if (keylen > 0) {
        uint8_t block[CFX_BLAKE2B_BLOCKBYTES] = {0};
        memcpy(block, key, keylen);
        cfx_blake2b_update(ctx, block, CFX_BLAKE2B_BLOCKBYTES);
        memset(block, 0, sizeof(block));
    }

    return 0;
}

int cfx_blake2b_update(cfx_blake2b_ctx_t* ctx, const void* in, size_t inlen) {
    cfx_blake2b_state_t* S = (cfx_blake2b_state_t*)ctx;
    const uint8_t* pin = (const uint8_t*)in;

    if (inlen == 0) return 0;

    size_t left = S->buflen;
    size_t fill = CFX_BLAKE2B_BLOCKBYTES - left;

    if (inlen > fill) {
        S->buflen = 0;
        memcpy(S->buf + left, pin, fill);
        blake2b_increment_counter(S, CFX_BLAKE2B_BLOCKBYTES);
        blake2b_compress_fn(S, S->buf);
        pin += fill;
        inlen -= fill;

        while (inlen > CFX_BLAKE2B_BLOCKBYTES) {
            blake2b_increment_counter(S, CFX_BLAKE2B_BLOCKBYTES);
            blake2b_compress_fn(S, pin);
            pin += CFX_BLAKE2B_BLOCKBYTES;
            inlen -= CFX_BLAKE2B_BLOCKBYTES;
        }
    }

    memcpy(S->buf + S->buflen, pin, inlen);
    S->buflen += inlen;

    return 0;
}

int cfx_blake2b_final(cfx_blake2b_ctx_t* ctx, void* out) {
    cfx_blake2b_state_t* S = (cfx_blake2b_state_t*)ctx;
    uint8_t buffer[CFX_BLAKE2B_OUTBYTES];

    blake2b_increment_counter(S, S->buflen);
    S->f[0] = ~(uint64_t)0;

    memset(S->buf + S->buflen, 0, CFX_BLAKE2B_BLOCKBYTES - S->buflen);
    blake2b_compress_fn(S, S->buf);

    for (int i = 0; i < 8; i++) {
        store64_le(buffer + i * 8, S->h[i]);
    }

    memcpy(out, buffer, S->outlen);
    memset(S, 0, sizeof(*S));

    return 0;
}

int cfx_blake2b(void* out, size_t outlen,
                const void* in, size_t inlen,
                const void* key, size_t keylen) {
    cfx_blake2b_ctx_t ctx;
    int ret;

    if (key && keylen > 0) {
        ret = cfx_blake2b_init_key(&ctx, outlen, key, keylen);
    } else {
        ret = cfx_blake2b_init(&ctx, outlen);
    }
    if (ret != 0) return ret;

    ret = cfx_blake2b_update(&ctx, in, inlen);
    if (ret != 0) return ret;

    return cfx_blake2b_final(&ctx, out);
}

/* ========================================================================== */
/* BLAKE2s */
/* ========================================================================== */

#define B2S_G(a, b, c, d, x, y) do { \
    a = a + b + x; \
    d = rotr32(d ^ a, 16); \
    c = c + d; \
    b = rotr32(b ^ c, 12); \
    a = a + b + y; \
    d = rotr32(d ^ a, 8); \
    c = c + d; \
    b = rotr32(b ^ c, 7); \
} while(0)

static void blake2s_compress(cfx_blake2s_state_t* S, const uint8_t block[64]) {
    uint32_t m[16], v[16];

    for (int i = 0; i < 16; i++) {
        m[i] = load32_le(block + i * 4);
    }

    for (int i = 0; i < 8; i++) {
        v[i] = S->h[i];
        v[i + 8] = blake2s_IV[i];
    }

    v[12] ^= S->t[0];
    v[13] ^= S->t[1];
    v[14] ^= S->f[0];
    v[15] ^= S->f[1];

    for (int r = 0; r < 10; r++) {
        const uint8_t* s = sigma[r];
        B2S_G(v[0], v[4], v[ 8], v[12], m[s[ 0]], m[s[ 1]]);
        B2S_G(v[1], v[5], v[ 9], v[13], m[s[ 2]], m[s[ 3]]);
        B2S_G(v[2], v[6], v[10], v[14], m[s[ 4]], m[s[ 5]]);
        B2S_G(v[3], v[7], v[11], v[15], m[s[ 6]], m[s[ 7]]);
        B2S_G(v[0], v[5], v[10], v[15], m[s[ 8]], m[s[ 9]]);
        B2S_G(v[1], v[6], v[11], v[12], m[s[10]], m[s[11]]);
        B2S_G(v[2], v[7], v[ 8], v[13], m[s[12]], m[s[13]]);
        B2S_G(v[3], v[4], v[ 9], v[14], m[s[14]], m[s[15]]);
    }

    for (int i = 0; i < 8; i++) {
        S->h[i] ^= v[i] ^ v[i + 8];
    }
}

static void blake2s_increment_counter(cfx_blake2s_state_t* S, uint32_t inc) {
    S->t[0] += inc;
    S->t[1] += (S->t[0] < inc) ? 1 : 0;
}

int cfx_blake2s_init(cfx_blake2s_ctx_t* ctx, size_t outlen) {
    if (outlen == 0 || outlen > CFX_BLAKE2S_OUTBYTES) return -1;

    cfx_blake2s_state_t* S = (cfx_blake2s_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    for (int i = 0; i < 8; i++) {
        S->h[i] = blake2s_IV[i];
    }

    S->h[0] ^= 0x01010000 ^ (uint32_t)outlen;
    S->outlen = outlen;

    return 0;
}

int cfx_blake2s_init_key(cfx_blake2s_ctx_t* ctx, size_t outlen,
                         const void* key, size_t keylen) {
    if (outlen == 0 || outlen > CFX_BLAKE2S_OUTBYTES) return -1;
    if (keylen > CFX_BLAKE2S_KEYBYTES) return -1;

    cfx_blake2s_state_t* S = (cfx_blake2s_state_t*)ctx;
    memset(S, 0, sizeof(*S));
    for (int i = 0; i < 8; i++) {
        S->h[i] = blake2s_IV[i];
    }

    S->h[0] ^= 0x01010000 ^ ((uint32_t)keylen << 8) ^ (uint32_t)outlen;
    S->outlen = outlen;

    if (keylen > 0) {
        uint8_t block[CFX_BLAKE2S_BLOCKBYTES] = {0};
        memcpy(block, key, keylen);
        cfx_blake2s_update(ctx, block, CFX_BLAKE2S_BLOCKBYTES);
        memset(block, 0, sizeof(block));
    }

    return 0;
}

int cfx_blake2s_update(cfx_blake2s_ctx_t* ctx, const void* in, size_t inlen) {
    cfx_blake2s_state_t* S = (cfx_blake2s_state_t*)ctx;
    const uint8_t* pin = (const uint8_t*)in;

    if (inlen == 0) return 0;

    size_t left = S->buflen;
    size_t fill = CFX_BLAKE2S_BLOCKBYTES - left;

    if (inlen > fill) {
        S->buflen = 0;
        memcpy(S->buf + left, pin, fill);
        blake2s_increment_counter(S, CFX_BLAKE2S_BLOCKBYTES);
        blake2s_compress(S, S->buf);
        pin += fill;
        inlen -= fill;

        while (inlen > CFX_BLAKE2S_BLOCKBYTES) {
            blake2s_increment_counter(S, CFX_BLAKE2S_BLOCKBYTES);
            blake2s_compress(S, pin);
            pin += CFX_BLAKE2S_BLOCKBYTES;
            inlen -= CFX_BLAKE2S_BLOCKBYTES;
        }
    }

    memcpy(S->buf + S->buflen, pin, inlen);
    S->buflen += inlen;

    return 0;
}

int cfx_blake2s_final(cfx_blake2s_ctx_t* ctx, void* out) {
    cfx_blake2s_state_t* S = (cfx_blake2s_state_t*)ctx;
    uint8_t buffer[CFX_BLAKE2S_OUTBYTES];

    blake2s_increment_counter(S, (uint32_t)S->buflen);
    S->f[0] = ~(uint32_t)0;

    memset(S->buf + S->buflen, 0, CFX_BLAKE2S_BLOCKBYTES - S->buflen);
    blake2s_compress(S, S->buf);

    for (int i = 0; i < 8; i++) {
        store32_le(buffer + i * 4, S->h[i]);
    }

    memcpy(out, buffer, S->outlen);
    memset(S, 0, sizeof(*S));

    return 0;
}

int cfx_blake2s(void* out, size_t outlen,
                const void* in, size_t inlen,
                const void* key, size_t keylen) {
    cfx_blake2s_ctx_t ctx;
    int ret;

    if (key && keylen > 0) {
        ret = cfx_blake2s_init_key(&ctx, outlen, key, keylen);
    } else {
        ret = cfx_blake2s_init(&ctx, outlen);
    }
    if (ret != 0) return ret;

    ret = cfx_blake2s_update(&ctx, in, inlen);
    if (ret != 0) return ret;

    return cfx_blake2s_final(&ctx, out);
}
