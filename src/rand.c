#include "cfx/rand.h"

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>


typedef struct {
    uint8_t  key[32];
    uint8_t  nonce[12];
    uint64_t counter;
    uint8_t  buf[64];
    unsigned idx;       /* next unread byte in buf (0..64) */
    int      seeded;
} cfx_rand_state_t;

static cfx_rand_state_t G = {0};

cfx_limb_t cfx_xorshift64(cfx_limb_t* s) {
    cfx_limb_t x = *s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *s = x;
    return x;
}

cfx_limb_t cfx_rand_limb(void) {
    return (cfx_limb_t)rand();
}

static int os_getrandom(void *out, size_t len) {
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd < 0) return -1;
    size_t got = 0;
    while (got < len) {
        ssize_t r = read(fd, (unsigned char*)out + got, len - got);
        if (r < 0) {
            if (errno == EINTR) continue;
            close(fd);
            return -1;
        }
        if (r == 0) { close(fd); return -1; }
        got += (size_t)r;
    }
    close(fd);
    return 0;
}

/* make it impossible to optimize away a memory clear to avoid dead-store elimination,
so we dont leave anything hanging around in RAM */
static void memzero_s(void* p, size_t n) {
    volatile unsigned char* v = (unsigned char*)p;
    while (n--) *v++ = 0;
}

/* Fill state from OS RNG (key, nonce, counter) */
void cfx_srand_os(void) {
    uint8_t tmp[32 + 12 + 8];
    if (os_getrandom(tmp, sizeof tmp) != 0) {
        /* TODO: behavior if OS RNG unavailable */
        memset(tmp, 0, sizeof tmp);
    }
    memcpy(G.key,   tmp,       32);
    memcpy(G.nonce, tmp + 32,  12);
    G.counter = 0;
    for (int i = 0; i < 8; ++i) G.counter |= (uint64_t)tmp[44 + i] << (8*i);
    memzero_s(tmp, sizeof tmp);
    G.idx = 64; /* force refill on first use */
    G.seeded = 1;
}


/* SplitMix64 for deterministic seed expansion */
static uint64_t splitmix64_next(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

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
static void store64_le(void *p, uint64_t x) {
    unsigned char *b = (unsigned char*)p;
    b[0] = (unsigned char)(x      );
    b[1] = (unsigned char)(x >>  8);
    b[2] = (unsigned char)(x >> 16);
    b[3] = (unsigned char)(x >> 24);
    b[4] = (unsigned char)(x >> 32);
    b[5] = (unsigned char)(x >> 40);
    b[6] = (unsigned char)(x >> 48);
    b[7] = (unsigned char)(x >> 56);
}

/* Deterministic seeding from unsigned int (std srand compatibility) */
void cfx_srand(unsigned int seed) {
    uint64_t s = ((uint64_t)seed << 32) | 0xA5A5A5A5u;
    /* Fill key, nonce, counter from SplitMix64 stream */
    for (int i = 0; i < 4; ++i) {
        uint64_t x = splitmix64_next(&s);
        store64_le(G.key + 8*i, x);
    }
    uint64_t n0 = splitmix64_next(&s);
    memcpy(G.nonce, &n0, 12); /* first 12 bytes */
    G.counter = splitmix64_next(&s);
    G.idx = 64;   /* force refill */
    G.seeded = 1;
}

#define ROTL32(x,n) ((uint32_t)(((x) << (n)) | ((x) >> (32-(n)))))

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

    memzero_s(ks, sizeof ks);
}


/* Internal: ensure we have bytes buffered */
static void refill_if_needed(void) {
    if (!G.seeded) cfx_srand_os();
    if (G.idx >= 64) {
        cfx_chacha20_block_rfc8439(G.key, G.counter++, G.nonce, G.buf);
        G.idx = 0;
    }
}

/* Return 31 random bits as int in [0, 0x7fffffff] */
int cfx_rand(void) {
    /* Pull 4 bytes, assemble LE uint32, mask to 31 bits */
    uint8_t b[4];
    for (int i = 0; i < 4; ++i) {
        refill_if_needed();
        b[i] = G.buf[G.idx++];
    }
    uint32_t v = ((uint32_t)b[0]) | ((uint32_t)b[1] << 8) |
                 ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
    return (int)(v & 0x7fffffff);
}
