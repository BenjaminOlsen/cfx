#include "cfx/chacha20.h"
#include "cfx/memory.h"

#include <string.h>

#define ROTL32(x,n) ((uint32_t)(((x) << (n)) | ((x) >> (32-(n)))))

#define WRITE_TO_PTR 0

void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], uint8_t out[64]) {
    static const uint32_t C[4] = {0x61707865u, 0x3320646eu, 0x79622d32u, 0x6b206574u};  /* "expa" "nd 3" "2-by" "te k" */
    uint32_t s[16], w[16];

    s[0]  = C[0];
    s[1]  = C[1];
    s[2]  = C[2];
    s[3]  = C[3];

    s[4]  = CFX_LOAD32_LE(key + 0);
    s[5]  = CFX_LOAD32_LE(key + 4);
    s[6]  = CFX_LOAD32_LE(key + 8);
    s[7]  = CFX_LOAD32_LE(key + 12);
    s[8]  = CFX_LOAD32_LE(key + 16);
    s[9]  = CFX_LOAD32_LE(key + 20);
    s[10] = CFX_LOAD32_LE(key + 24);
    s[11] = CFX_LOAD32_LE(key + 28);

    s[12] = counter;                            /* 32-bit block counter */

    s[13] = CFX_LOAD32_LE(nonce + 0);
    s[14] = CFX_LOAD32_LE(nonce + 4);
    s[15] = CFX_LOAD32_LE(nonce + 8);           /* 96-bit nonce */

    /* for (size_t i = 0; i < 16; ++i) w[i] = s[i]; */
    memcpy(w, s, sizeof(w));

/* Quarter Round */
#define QR(a, b, c, d) \
    a += b; d ^= a; d = ROTL32(d,16); \
    c += d; b ^= c; b = ROTL32(b,12); \
    a += b; d ^= a; d = ROTL32(d, 8); \
    c += d; b ^= c; b = ROTL32(b, 7);

    for (size_t i = 0; i < 10; ++i){
        /* column rounds */
        QR(w[0], w[4], w[8], w[12])
        QR(w[1], w[5], w[9], w[13])
        QR(w[2], w[6], w[10], w[14])
        QR(w[3], w[7], w[11], w[15])
        /* diagonal rounds */
        QR(w[0], w[5], w[10], w[15])
        QR(w[1], w[6], w[11], w[12])
        QR(w[2], w[7], w[8], w[13])
        QR(w[3], w[4], w[9], w[14])
    }
#undef QR

    for (size_t i = 0; i < 16; ++i) w[i] += s[i];

    for (size_t i = 0; i < 16; ++i) CFX_STORE32_LE(out + 4 * i, w[i]);
}


void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t *pt, size_t pt_len, uint8_t *ct) {
    uint8_t ks[64];

    while (pt_len) {
        cfx_chacha20_block_rfc8439(key, counter, nonce, ks);
        ++counter;  /* increment per 64-byte block */

        size_t take = pt_len < 64 ? pt_len : 64;
        for (size_t i = 0; i < take; ++i) {
            ct[i] = pt[i] ^ ks[i];
        }

        pt += take;
        ct += take;
        pt_len  -= take;
    }

    CFX_MEMZERO_S(ks, sizeof(ks));
}

/* ---------------------------------------------------------------------------------------------- */
/* Here be SIMD */

typedef struct {
    uint32_t v[4];
} vec4_u32;

static inline vec4_u32 v_add(vec4_u32 a, vec4_u32 b) {
    vec4_u32 r;
    r.v[0] = a.v[0] + b.v[0];
    r.v[1] = a.v[1] + b.v[1];
    r.v[2] = a.v[2] + b.v[2];
    r.v[3] = a.v[3] + b.v[3];
    return r;
}


static inline vec4_u32 v_xor(vec4_u32 a, vec4_u32 b) {
    vec4_u32 r;
    r.v[0] = a.v[0] ^ b.v[0];
    r.v[1] = a.v[1] ^ b.v[1];
    r.v[2] = a.v[2] ^ b.v[2];
    r.v[3] = a.v[3] ^ b.v[3];
    return r;
}

static inline vec4_u32 v_rotl(vec4_u32 x, unsigned n) {
    vec4_u32 r;
    r.v[0] = ROTL32(x.v[0], n);
    r.v[1] = ROTL32(x.v[1], n);
    r.v[2] = ROTL32(x.v[2], n);
    r.v[3] = ROTL32(x.v[3], n);
    return r;
}

static inline vec4_u32 v_set1(uint32_t x) {
    vec4_u32 r = {{ x, x, x, x }};
    return r;
}

static inline vec4_u32 v_set4(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    vec4_u32 r = {{ a, b, c, d }};
    return r;
}

static inline vec4_u32 v_load32_le_1(const void* x) {
    uint32_t w = CFX_LOAD32_LE(x);
    return v_set1(w);
}


static inline vec4_u32 v_load32_le_4(const void* a, const void* b, const void* c, const void* d) {
    vec4_u32 r;
    r.v[0] = CFX_LOAD32_LE(a);
    r.v[1] = CFX_LOAD32_LE(b);
    r.v[2] = CFX_LOAD32_LE(c);
    r.v[3] = CFX_LOAD32_LE(d);
    return r;
}

typedef struct {
    vec4_u32 x[16];  // x[0] is word0 for 4 blocks, etc.
} chacha_state4;


void cfx_chacha20_block4_simd(const uint8_t key[32], const uint32_t counter[4],
                              const uint8_t nonce[4][12], uint8_t out[4][64]) {
    static const uint32_t C[4] = {0x61707865u, 0x3320646eu, 0x79622d32u, 0x6b206574u};  /* "expa" "nd 3" "2-by" "te k" */
    vec4_u32 s[16], w[16];

    s[0]  = v_set1(C[0]);
    s[1]  = v_set1(C[1]);
    s[2]  = v_set1(C[2]);
    s[3]  = v_set1(C[3]);

    s[4]  = v_load32_le_1(key + 0);
    s[5]  = v_load32_le_1(key + 4);
    s[6]  = v_load32_le_1(key + 8);
    s[7]  = v_load32_le_1(key + 12);
    s[8]  = v_load32_le_1(key + 16);
    s[9]  = v_load32_le_1(key + 20);
    s[10] = v_load32_le_1(key + 24);
    s[11] = v_load32_le_1(key + 28);

    s[12] = v_set4(counter[0], counter[1], counter[2], counter[3]);

    s[13] = v_load32_le_4(nonce[0] + 0, nonce[1] + 0, nonce[2] + 0, nonce[3] + 0);
    s[14] = v_load32_le_4(nonce[0] + 4, nonce[1] + 4, nonce[2] + 4, nonce[3] + 4);
    s[15] = v_load32_le_4(nonce[0] + 8, nonce[1] + 8, nonce[2] + 8, nonce[3] + 8);
    

    /* for (size_t i = 0; i < 16; ++i) w[i] = s[i]; */
    memcpy(w, s, sizeof(w));


#define VQR(a,b,c,d)  \
    a = v_add(a, b); d = v_xor(d, a); d = v_rotl(d,16); \
    c = v_add(c, d); b = v_xor(b, c); b = v_rotl(b,12); \
    a = v_add(a, b); d = v_xor(d, a); d = v_rotl(d, 8); \
    c = v_add(c, d); b = v_xor(b, c); b = v_rotl(b, 7);

    for (size_t i = 0; i < 10; ++i){
        /* column rounds */
        VQR(w[0], w[4], w[8], w[12])
        VQR(w[1], w[5], w[9], w[13])
        VQR(w[2], w[6], w[10], w[14])
        VQR(w[3], w[7], w[11], w[15])
        /* diagonal rounds */
        VQR(w[0], w[5], w[10], w[15])
        VQR(w[1], w[6], w[11], w[12])
        VQR(w[2], w[7], w[8], w[13])
        VQR(w[3], w[4], w[9], w[14])
    }
#undef VQR

    for (size_t i = 0; i < 16; ++i) {
        w[i] = v_add(w[i], s[i]);
    }

    for (size_t i = 0; i < 16; ++i) {
        CFX_STORE32_LE(out[0] + 4*i, w[i].v[0]);
        CFX_STORE32_LE(out[1] + 4*i, w[i].v[1]);
        CFX_STORE32_LE(out[2] + 4*i, w[i].v[2]);
        CFX_STORE32_LE(out[3] + 4*i, w[i].v[3]);
    }
}

