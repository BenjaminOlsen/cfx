#ifndef CFX_CHACHA20_H
#define CFX_CHACHA20_H

#include <stdint.h>
#include <stddef.h>

#include "cfx/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------------------------------------------------ */
/* ChaCha20
 * refs
 *  https://datatracker.ietf.org/doc/html/rfc8439
 *  https://loup-vaillant.fr/tutorials/chacha20-design
 */

#define CFX_CHACHA20_ROTL32(x,n) ((uint32_t)(((x) << (n)) | ((x) >> (32-(n)))))

/* Quarter Round */
#define CFX_CHACHA20_QR(a, b, c, d) \
    a += b; d ^= a; d = ROTL32(d,16); \
    c += d; b ^= c; b = ROTL32(b,12); \
    a += b; d ^= a; d = ROTL32(d, 8); \
    c += d; b ^= c; b = ROTL32(b, 7);

typedef struct {
    uint32_t s[16];
} cfx_chacha_state_t;

/* ---------------------------------------------------------------------------------------------- */
/* SIMD */
#if CFX_SIMD

typedef uint32_t vec4_u32 __attribute__((vector_size(16)));

static inline vec4_u32 v_add(vec4_u32 a, vec4_u32 b) {
    return a + b;
}

static inline vec4_u32 v_xor(vec4_u32 a, vec4_u32 b) {
    return a ^ b;
}

static inline vec4_u32 v_rotl(vec4_u32 x, unsigned n) {
    return (x << n) | (x >> (32 - n));
}

static inline vec4_u32 v_set1(uint32_t x) {
    return (vec4_u32){ x, x, x, x };
}

static inline vec4_u32 v_set4(uint32_t a, uint32_t b, uint32_t c, uint32_t d) {
    return (vec4_u32){ a, b, c, d };
}

static inline vec4_u32 v_load32_le_1(const void* p) {
    uint32_t w = CFX_LOAD32_LE(p);
    return v_set1(w);
}

static inline vec4_u32 v_load32_le_4(const void* a, const void* b,
                                     const void* c, const void* d) {
    return (vec4_u32){
        CFX_LOAD32_LE(a),
        CFX_LOAD32_LE(b),
        CFX_LOAD32_LE(c),
        CFX_LOAD32_LE(d)
    };
}

static inline void v_store32_le_4(uint8_t *out0, uint8_t *out1,
                                  uint8_t *out2, uint8_t *out3, vec4_u32 x) {
    CFX_STORE32_LE(out0, x[0]);
    CFX_STORE32_LE(out1, x[1]);
    CFX_STORE32_LE(out2, x[2]);
    CFX_STORE32_LE(out3, x[3]);
}

#else

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


static inline void v_store32_le_4(uint8_t *out0, uint8_t *out1,
                                  uint8_t *out2, uint8_t *out3, vec4_u32 x) {
    CFX_STORE32_LE(out0, x.v[0]);
    CFX_STORE32_LE(out1, x.v[1]);
    CFX_STORE32_LE(out2, x.v[2]);
    CFX_STORE32_LE(out3, x.v[3]);
}

#endif

typedef struct {
    vec4_u32 s[16];  /*s[0] is word0 for 4 blocks, etc. */
} cfx_chacha_state4_t;

#define CFX_CHACHA20_VQR(a,b,c,d)  \
    a = v_add(a, b); d = v_xor(d, a); d = v_rotl(d,16); \
    c = v_add(c, d); b = v_xor(b, c); b = v_rotl(b,12); \
    a = v_add(a, b); d = v_xor(d, a); d = v_rotl(d, 8); \
    c = v_add(c, d); b = v_xor(b, c); b = v_rotl(b, 7);



void cfx_chacha20_state_init(cfx_chacha_state_t* ctx, const uint8_t key[32], const uint8_t nonce[12]);
void cfx_chacha20_block(cfx_chacha_state_t* ctx, const uint32_t counter, uint8_t out[64]);

void cfx_chacha20_state_init4(cfx_chacha_state4_t* ctx, const uint8_t key[32], const uint8_t nonce[4][12]);
void cfx_chacha20_block4(cfx_chacha_state4_t* ctx, const uint32_t counter[4], uint8_t out[4][64]);

void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                                uint8_t out[64]);

void cfx_chacha20_block_rfc8439_2(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                                uint8_t out[64]);

void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t *pt, size_t pt_len, uint8_t *ct);

void cfx_chacha20_encrypt_bytes(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t *pt, size_t pt_len, uint8_t *ct);

void cfx_chacha20_block4_simd(const uint8_t key[32], const uint32_t counter[4],
                              const uint8_t nonce[4][12], uint8_t out[4][64]);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_CHACHA20_H */
