#include "cfx/crypto.h"
#include "cfx/memory.h"

#define ROTL32(x,n) ((uint32_t)(((x) << (n)) | ((x) >> (32-(n)))))

static uint32_t load32_le(const void *p) {
    const unsigned char *b = (const unsigned char*)p;
    return ((uint32_t)b[0]) | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

static void store32_le(void *p, uint32_t x) {
    unsigned char *b = (unsigned char*)p;
    b[0] = (unsigned char)(x      );
    b[1] = (unsigned char)(x >>  8);
    b[2] = (unsigned char)(x >> 16);
    b[3] = (unsigned char)(x >> 24);
}

void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], uint8_t out[64]) {
    static const uint32_t C[4] = {0x61707865u, 0x3320646eu, 0x79622d32u, 0x6b206574u}; /* "expa" "nd 3" "2-by" "te k" */
    uint32_t s[16], w[16];

    s[0] = C[0];
    s[1] = C[1];
    s[2] = C[2];
    s[3] = C[3];
    s[4] = load32_le(key + 0);
    s[5] = load32_le(key + 4);
    s[6] = load32_le(key + 8);
    s[7] = load32_le(key + 12);
    s[8] = load32_le(key + 16);
    s[9] = load32_le(key + 20);
    s[10] = load32_le(key + 24);
    s[11] = load32_le(key + 28);

    s[12] = counter;                          /* 32-bit block counter */
    s[13] = load32_le(nonce+0);
    s[14] = load32_le(nonce+4);
    s[15] = load32_le(nonce+8);               /* 96-bit nonce */

    for (int i=0;i<16;++i) w[i]=s[i];

#define QR(a, b, c, d) \
    a += b; d ^= a; d = ROTL32(d,16); \
    c += d; b ^= c; b = ROTL32(b,12); \
    a += b; d ^= a; d = ROTL32(d, 8); \
    c += d; b ^= c; b = ROTL32(b, 7);

    for (int i=0;i<10;++i){
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

    for (int i = 0; i < 16; ++i) w[i] += s[i];
    for (int i = 0; i < 16; ++i) store32_le(out + 4 * i, w[i]);
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

    cfx_memzero_s(ks, sizeof ks);
}

