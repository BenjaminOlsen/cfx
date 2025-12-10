#include "cfx/chacha20.h"
#include "cfx/arch.h"

#if CFX_HAVE_AVX2
#include <immintrin.h>
#endif

#include <string.h>


typedef struct {
    uint32_t s[16];
} cfx_chacha20_state_t;


typedef struct {
    vec4_u32 s[16];  /*s[0] is word0 for 4 blocks, etc. */
} cfx_chacha20_state4_t;

#define ROTL32(x, n) CFX_CHACHA20_ROTL32(x, n)
#define QR(a, b, c, d) CFX_CHACHA20_QR(a, b, c, d)

#define _EXPA 0x61707865u
#define _ND_3 0x3320646eu
#define _2_BY 0x79622d32u
#define _TE_K 0x6b206574u

void cfx_chacha20_ctx_init(cfx_chacha20_ctx_t* ctx, const uint8_t key[32], const uint8_t nonce[12]) {
    cfx_chacha20_state_t* st = (cfx_chacha20_state_t*)ctx->opaque;

    st->s[0]  = _EXPA;
    st->s[1]  = _ND_3;
    st->s[2]  = _2_BY;
    st->s[3]  = _TE_K;
    st->s[4]  = CFX_LOAD32_LE(key + 0);
    st->s[5]  = CFX_LOAD32_LE(key + 4);
    st->s[6]  = CFX_LOAD32_LE(key + 8);
    st->s[7]  = CFX_LOAD32_LE(key + 12);
    st->s[8]  = CFX_LOAD32_LE(key + 16);
    st->s[9]  = CFX_LOAD32_LE(key + 20);
    st->s[10] = CFX_LOAD32_LE(key + 24);
    st->s[11] = CFX_LOAD32_LE(key + 28);
    st->s[12] = 0;     /* counter */
    st->s[13] = CFX_LOAD32_LE(nonce + 0);
    st->s[14] = CFX_LOAD32_LE(nonce + 4);
    st->s[15] = CFX_LOAD32_LE(nonce + 8);
}

void cfx_chacha20_block(cfx_chacha20_ctx_t* ctx, uint32_t counter, uint8_t out[64]) {
    cfx_chacha20_state_t* st = (cfx_chacha20_state_t*)ctx->opaque;
    uint32_t w[16];
    st->s[12] = counter;
    memcpy(w, st->s, sizeof(w));

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

    for (size_t i = 0; i < 16; ++i) w[i] += st->s[i];
    for (size_t i = 0; i < 16; ++i) CFX_STORE32_LE(out + 4 * i, w[i]);
}

void cfx_chacha20_encrypt_ctx(cfx_chacha20_ctx_t* ctx, uint32_t* const counter, const uint8_t* pt, size_t pt_len, uint8_t* ct) {
    uint8_t ks[64];

    while (pt_len) {
        cfx_chacha20_block(ctx, *counter, ks);
        ++*counter;

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


void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], uint8_t out[64]) {
    uint32_t s[16], w[16];

    s[0]  = 0x61707865u;    /* "expa" */
    s[1]  = 0x3320646eu;    /* "nd 3" */
    s[2]  = 0x79622d32u;    /* "2-by" */
    s[3]  = 0x6b206574u;    /* "te k" */

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

    for (size_t i = 0; i < 16; ++i) w[i] += s[i];
    for (size_t i = 0; i < 16; ++i) CFX_STORE32_LE(out + 4 * i, w[i]);
}


void cfx_chacha20_block_rfc8439_2(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], uint8_t out[64]) {

    uint32_t s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15;
    uint32_t w0,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10,w11,w12,w13,w14,w15;

    s0  = 0x61707865u;  /* "expa" */
    s1  = 0x3320646eu;  /* "nd 3" */
    s2  = 0x79622d32u;  /* "2-by" */
    s3  = 0x6b206574u;  /* "te k" */

    s4  = CFX_LOAD32_LE(key + 0);
    s5  = CFX_LOAD32_LE(key + 4);
    s6  = CFX_LOAD32_LE(key + 8);
    s7  = CFX_LOAD32_LE(key + 12);
    s8  = CFX_LOAD32_LE(key + 16);
    s9  = CFX_LOAD32_LE(key + 20);
    s10 = CFX_LOAD32_LE(key + 24);
    s11 = CFX_LOAD32_LE(key + 28);

    s12 = counter;                            /* 32-bit block counter */

    s13 = CFX_LOAD32_LE(nonce + 0);
    s14 = CFX_LOAD32_LE(nonce + 4);
    s15 = CFX_LOAD32_LE(nonce + 8);           /* 96-bit nonce */

    w0  = s0 ;
    w1  = s1 ;
    w2  = s2 ;
    w3  = s3 ;
    w4  = s4 ;
    w5  = s5 ;
    w6  = s6 ;
    w7  = s7 ;
    w8  = s8 ;
    w9  = s9 ;
    w10 = s10;
    w11 = s11;
    w12 = s12;
    w13 = s13;
    w14 = s14;
    w15 = s15;

    for (size_t i = 20; i; i-=2){
        /* column rounds */
        QR(w0, w4, w8, w12)
        QR(w1, w5, w9, w13)
        QR(w2, w6, w10, w14)
        QR(w3, w7, w11, w15)
        /* diagonal rounds */
        QR(w0, w5, w10, w15)
        QR(w1, w6, w11, w12)
        QR(w2, w7, w8, w13)
        QR(w3, w4, w9, w14)
    }

    w0  += s0 ;
    w1  += s1 ;
    w2  += s2 ;
    w3  += s3 ;
    w4  += s4 ;
    w5  += s5 ;
    w6  += s6 ;
    w7  += s7 ;
    w8  += s8 ;
    w9  += s9 ;
    w10 += s10;
    w11 += s11;
    w12 += s12;
    w13 += s13;
    w14 += s14;
    w15 += s15;

    CFX_STORE32_LE(out + 0, w0);
    CFX_STORE32_LE(out + 4, w1);
    CFX_STORE32_LE(out + 8, w2);
    CFX_STORE32_LE(out + 12, w3);
    CFX_STORE32_LE(out + 16, w4);
    CFX_STORE32_LE(out + 20, w5);
    CFX_STORE32_LE(out + 24, w6);
    CFX_STORE32_LE(out + 28, w7);
    CFX_STORE32_LE(out + 32, w8);
    CFX_STORE32_LE(out + 36, w9);
    CFX_STORE32_LE(out + 40, w10);
    CFX_STORE32_LE(out + 44, w11);
    CFX_STORE32_LE(out + 48, w12);
    CFX_STORE32_LE(out + 52, w13);
    CFX_STORE32_LE(out + 56, w14);
    CFX_STORE32_LE(out + 60, w15);
}

void cfx_chacha20_block_rfc8439_3(const uint8_t *CFX_RESTRICT key,
                                  uint32_t counter,
                                  const uint8_t *CFX_RESTRICT nonce,
                                  uint8_t *CFX_RESTRICT out) {

    uint32_t x0  = _EXPA;
    uint32_t x1  = _ND_3;
    uint32_t x2  = _2_BY;
    uint32_t x3  = _TE_K;
    uint32_t x4  = CFX_LOAD32_LE(key + 0);
    uint32_t x5  = CFX_LOAD32_LE(key + 4);
    uint32_t x6  = CFX_LOAD32_LE(key + 8);
    uint32_t x7  = CFX_LOAD32_LE(key + 12);
    uint32_t x8  = CFX_LOAD32_LE(key + 16);
    uint32_t x9  = CFX_LOAD32_LE(key + 20);
    uint32_t x10 = CFX_LOAD32_LE(key + 24);
    uint32_t x11 = CFX_LOAD32_LE(key + 28);
    uint32_t x12 = counter;
    uint32_t x13 = CFX_LOAD32_LE(nonce + 0);
    uint32_t x14 = CFX_LOAD32_LE(nonce + 4);
    uint32_t x15 = CFX_LOAD32_LE(nonce + 8);

    uint32_t o0 = x0,  o1 = x1,  o2 = x2,  o3 = x3;
    uint32_t o4 = x4,  o5 = x5,  o6 = x6,  o7 = x7;
    uint32_t o8 = x8,  o9 = x9,  o10 = x10, o11 = x11;
    uint32_t o12 = x12, o13 = x13, o14 = x14, o15 = x15;

    for (unsigned i = 20; i; i -= 2) {
        QR(x0, x4, x8,  x12);
        QR(x1, x5, x9,  x13);
        QR(x2, x6, x10, x14);
        QR(x3, x7, x11, x15);
        QR(x0, x5, x10, x15);
        QR(x1, x6, x11, x12);
        QR(x2, x7, x8,  x13);
        QR(x3, x4, x9,  x14);
    }

    x0  += o0;  x1  += o1;  x2  += o2;  x3  += o3;
    x4  += o4;  x5  += o5;  x6  += o6;  x7  += o7;
    x8  += o8;  x9  += o9;  x10 += o10; x11 += o11;
    x12 += o12; x13 += o13; x14 += o14; x15 += o15;

    CFX_STORE32_LE(out +  0, x0);
    CFX_STORE32_LE(out +  4, x1);
    CFX_STORE32_LE(out +  8, x2);
    CFX_STORE32_LE(out + 12, x3);
    CFX_STORE32_LE(out + 16, x4);
    CFX_STORE32_LE(out + 20, x5);
    CFX_STORE32_LE(out + 24, x6);
    CFX_STORE32_LE(out + 28, x7);
    CFX_STORE32_LE(out + 32, x8);
    CFX_STORE32_LE(out + 36, x9);
    CFX_STORE32_LE(out + 40, x10);
    CFX_STORE32_LE(out + 44, x11);
    CFX_STORE32_LE(out + 48, x12);
    CFX_STORE32_LE(out + 52, x13);
    CFX_STORE32_LE(out + 56, x14);
    CFX_STORE32_LE(out + 60, x15);
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

void cfx_chacha20_encrypt_bytes(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t *pt, size_t pt_len, uint8_t *ct) {
    uint8_t ks[64];

    cfx_chacha20_ctx_t ctx;
    cfx_chacha20_ctx_init(&ctx, key, nonce);

    while (pt_len) {
        cfx_chacha20_block(&ctx, counter, ks);
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

/* ---------------------------------------------------------------------------------------------- */
/* block4 helpers */
#if CFX_SIMD

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

#define VQR(a,b,c,d) CFX_CHACHA20_VQR(a,b,c,d)


void cfx_chacha20_state_init4(cfx_chacha20_ctx4_t* ctx, const uint8_t key[32], const uint8_t nonce[4][12]) {
    cfx_chacha20_state4_t*st = (cfx_chacha20_state4_t*)ctx;
    st->s[0]  = v_set1(_EXPA);
    st->s[1]  = v_set1(_ND_3);
    st->s[2]  = v_set1(_2_BY);
    st->s[3]  = v_set1(_TE_K);

    st->s[4]  = v_load32_le_1(key + 0);
    st->s[5]  = v_load32_le_1(key + 4);
    st->s[6]  = v_load32_le_1(key + 8);
    st->s[7]  = v_load32_le_1(key + 12);
    st->s[8]  = v_load32_le_1(key + 16);
    st->s[9]  = v_load32_le_1(key + 20);
    st->s[10] = v_load32_le_1(key + 24);
    st->s[11] = v_load32_le_1(key + 28);

    st->s[12] = v_set1(0);

    st->s[13] = v_load32_le_4(nonce[0] + 0, nonce[1] + 0, nonce[2] + 0, nonce[3] + 0);
    st->s[14] = v_load32_le_4(nonce[0] + 4, nonce[1] + 4, nonce[2] + 4, nonce[3] + 4);
    st->s[15] = v_load32_le_4(nonce[0] + 8, nonce[1] + 8, nonce[2] + 8, nonce[3] + 8);
}

void cfx_chacha20_block4(cfx_chacha20_ctx4_t* ctx, const uint32_t counter[4], uint8_t out[4][64]) {
    cfx_chacha20_state4_t*st = (cfx_chacha20_state4_t*)ctx;
    vec4_u32 w[16];

    st->s[12] = v_set4(counter[0], counter[1], counter[2], counter[3]);

    memcpy(w, st->s, sizeof(w));

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

    for (size_t i = 0; i < 16; ++i) {
        w[i] = v_add(w[i], st->s[i]);
    }

    for (size_t i = 0; i < 16; ++i) {
        v_store32_le_4(out[0] + 4*i,
                       out[1] + 4*i,
                       out[2] + 4*i,
                       out[3] + 4*i,
                       w[i]);
    }
}


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

    for (size_t i = 0; i < 16; ++i) {
        w[i] = v_add(w[i], s[i]);
    }

    for (size_t i = 0; i < 16; ++i) {
        v_store32_le_4(out[0] + 4*i,
                       out[1] + 4*i,
                       out[2] + 4*i,
                       out[3] + 4*i,
                       w[i]);
    }
}

#if CFX_HAVE_AVX2

static inline __m256i cfx_mm256_rotl32(__m256i x, int n) {
    return _mm256_or_si256(_mm256_slli_epi32(x, n),
                           _mm256_srli_epi32(x, 32 - n));
}

#define MM256_QR(a,b,c,d) do {         \
    a = _mm256_add_epi32(a, b);     \
    d = _mm256_xor_si256(d, a);     \
    d = cfx_mm256_rotl32(d, 16);      \
    c = _mm256_add_epi32(c, d);     \
    b = _mm256_xor_si256(b, c);     \
    b = cfx_mm256_rotl32(b, 12);      \
    a = _mm256_add_epi32(a, b);     \
    d = _mm256_xor_si256(d, a);     \
    d = cfx_mm256_rotl32(d, 8);       \
    c = _mm256_add_epi32(c, d);     \
    b = _mm256_xor_si256(b, c);     \
    b = cfx_mm256_rotl32(b, 7);       \
} while (0)

/* This function takes 16 __m256i vectors x[0..15], each representing:
 *
 *   x[w] = [in[0][w], in[1][w], ..., in[7][w]]
 *
 * and reconstructs per-block AoS layout into out[8][16]:
 *
 *   out[b][w] = original in[b][w]
 *
 * using AVX2 + SSE transposes.
 */
#define LANE32(v, idx) (uint32_t)_mm_extract_epi32(v, idx)

/* x[w] = [ in[0][w], in[1][w], ..., in[7][w] ] */
static inline void transpose_16x8_to_blocks(const __m256i x[16], uint32_t out[8][16]) {
    __m128i lo[16], hi[16];

    for (int w = 0; w < 16; ++w) {
        /* lo[w] = _mm256_castsi256_si128(x[w]);     */
        lo[w] = _mm256_extracti128_si256(x[w], 0); /* lanes 0..3 */
        hi[w] = _mm256_extracti128_si256(x[w], 1);  /* lanes 4..7 */
    }

    uint32_t *dst;

#define EXTRACT_BLOCK(ARR, lane_idx, out_idx)                           \
    do {                                                                \
        dst = out[(out_idx)];                                           \
        for (int w = 0; w < 16; ++w) {                                  \
            dst[w] = LANE32((ARR)[w], (lane_idx));                      \
        }                                                               \
    } while (0)

    /* Blocks 0..3 from lo */
    EXTRACT_BLOCK(lo, 0, 0);
    EXTRACT_BLOCK(lo, 1, 1);
    EXTRACT_BLOCK(lo, 2, 2);
    EXTRACT_BLOCK(lo, 3, 3);

    /* Blocks 4..7 from hi */
    EXTRACT_BLOCK(hi, 0, 4);
    EXTRACT_BLOCK(hi, 1, 5);
    EXTRACT_BLOCK(hi, 2, 6);
    EXTRACT_BLOCK(hi, 3, 7);

#undef EXTRACT_BLOCK
#undef LANE32
}

void cfx_chacha20_block8_avx2(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], uint8_t out[8][64]) {

    /* Expand key/nonce to u32 */
    uint32_t k[8];
    uint32_t n[3];

    for (size_t i = 0; i < 8; ++i) {
        k[i] = CFX_LOAD32_LE(key + 4*i);
    }
    for (size_t i = 0; i < 3; ++i) {
        n[i] = CFX_LOAD32_LE(nonce + 4*i);
    }

    /* lane offsets  {0,1,2,3,4,5,6,7} */
    const __m256i lane = _mm256_setr_epi32(0,1,2,3,4,5,6,7);

    __m256i x[16], o[16];

    x[0]  = _mm256_set1_epi32(_EXPA);
    x[1]  = _mm256_set1_epi32(_ND_3);
    x[2]  = _mm256_set1_epi32(_2_BY);
    x[3]  = _mm256_set1_epi32(_TE_K);
    x[4]  = _mm256_set1_epi32((int32_t)k[0]);
    x[5]  = _mm256_set1_epi32((int32_t)k[1]);
    x[6]  = _mm256_set1_epi32((int32_t)k[2]);
    x[7]  = _mm256_set1_epi32((int32_t)k[3]);
    x[8]  = _mm256_set1_epi32((int32_t)k[4]);
    x[9]  = _mm256_set1_epi32((int32_t)k[5]);
    x[10] = _mm256_set1_epi32((int32_t)k[6]);
    x[11] = _mm256_set1_epi32((int32_t)k[7]);

    /* lane 0..7 has counter value: {counter, counter+1, ... counter+7} */
    x[12] = _mm256_add_epi32(_mm256_set1_epi32((int32_t)counter), lane);

    x[13] = _mm256_set1_epi32((int32_t)n[0]);
    x[14] = _mm256_set1_epi32((int32_t)n[1]);
    x[15] = _mm256_set1_epi32((int32_t)n[2]);

    /* originals */
    for (size_t i = 0; i < 16; ++i) {
        o[i] = x[i];
    }

    for (size_t i = 0; i < 10; ++i) {
        MM256_QR(x[0], x[4], x[8],  x[12]);
        MM256_QR(x[1], x[5], x[9],  x[13]);
        MM256_QR(x[2], x[6], x[10], x[14]);
        MM256_QR(x[3], x[7], x[11], x[15]);
        MM256_QR(x[0], x[5], x[10], x[15]);
        MM256_QR(x[1], x[6], x[11], x[12]);
        MM256_QR(x[2], x[7], x[8],  x[13]);
        MM256_QR(x[3], x[4], x[9],  x[14]);
    }

    for (size_t i = 0; i < 16; ++i) {
        x[i] = _mm256_add_epi32(x[i], o[i]);
    }

    uint32_t (*out_words)[16] = (uint32_t (*)[16])out;
    transpose_16x8_to_blocks(x, out_words);
}

#undef MM256_QR

#endif /* CFX_HAVE_AVX2 */



#undef VQR
#undef QR
