#include "cfx/rand.h"
#include "cfx/chacha20.h"
#include "cfx/poly1305.h"
#include "cfx/memory.h"

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>


typedef struct {
    uint8_t  buf[64];
    uint8_t  key[32];
    uint8_t  nonce[12];
    uint64_t counter;
    unsigned idx;       /* next unread byte in buf (0..64) */
    int      seeded;
} chacha_rand_state_t;

static chacha_rand_state_t G = {0};

/* ---------------------------------------------------------------------------------------------- */
/* uses splitmix64 to generate nonce, counter, key for chacha20*/
void cfx_srand(unsigned int seed) {
    uint64_t s = ((uint64_t)seed << 32) | 0xA5A5A5A5u;
    /* Fill key, nonce, counter from SplitMix64 stream */
    for (int i = 0; i < 4; ++i) {
        uint64_t x = cfx_splitmix64(&s);
        CFX_STORE64_LE(G.key + 8*i, x);
    }
    uint64_t n0 = cfx_splitmix64(&s);
    uint32_t n1 = (uint32_t) cfx_splitmix64(&s);

    CFX_STORE64_LE(G.nonce,      n0);
    CFX_STORE32_LE(G.nonce + 8,  n1);

    G.counter = cfx_splitmix64(&s);
    G.idx = 64;   /* force refill */
    G.seeded = 1;
}

static void refill_if_needed(void) {
    if (!G.seeded) cfx_srand_os();
    if (G.idx >= 64) {
        cfx_chacha20_block_rfc8439(G.key, G.counter++, G.nonce, G.buf);
        G.idx = 0;
    }
}

uint32_t cfx_urand(void) {
    uint8_t b[4];
    for (size_t i = 0; i < 4; ++i) {
        refill_if_needed();
        b[i] = G.buf[G.idx++];
    }
    uint32_t v = ((uint32_t)b[0]) | ((uint32_t)b[1] << 8) |
                 ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
    return v;
}

int cfx_rand(void) {
    return (int)(cfx_urand() & 0x7fffffff);
}

void cfx_randombytes(void* buf, size_t len) {
    uint8_t* p = (uint8_t*)buf;
    for (size_t i = 0; i < len; ++i) {
        refill_if_needed();
        p[i] = G.buf[G.idx++];
    }
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

/* Fill state from OS RNG (key, nonce, counter) */
void cfx_srand_os(void) {
    uint8_t tmp[32 + 12 + 8];
    if (os_getrandom(tmp, sizeof(tmp)) != 0) {
        /* TODO: behavior if OS RNG unavailable */
        memset(tmp, 0, sizeof(tmp));
    }
    memcpy(G.key,   tmp,       32);
    memcpy(G.nonce, tmp + 32,  12);
    G.counter = 0;
    for (int i = 0; i < 8; ++i) G.counter |= (uint64_t)tmp[44 + i] << (8*i);
    cfx_memzero_s(tmp, sizeof(tmp));
    G.idx = 64; /* force refill on first use */
    G.seeded = 1;
}

void cfx_randombytes_os(void* buf, size_t len) {
    if (os_getrandom(buf, len) != 0) {
        memset(buf, 0, len);
    }
}

cfx_limb_t cfx_rand_limb(void) {
    return (cfx_limb_t)rand();
}

/* ---------------------------------------------------------------------------------------------- */
uint32_t g_lcg_seed = 0x126;

void cfx_lcg_seed(unsigned seed) {
    g_lcg_seed = seed;
}

uint32_t cfx_lcg_gen(void) {
    return cfx_lcg(&g_lcg_seed);
}

uint32_t cfx_lcg(uint32_t* s) {
    *s *= 1664525u;
    *s += 1013904223u;
    return *s;
}

void cfx_lcg_bytes(uint32_t seed, uint8_t *data, size_t len) {
    /*
    ref: LCG - https://en.wikipedia.org/wiki/Linear_congruential_generator
    */
    uint32_t x = seed;
    uint32_t i;

    for (i = 0; i < len; i++) {
        x = 1664525u * x + 1013904223u;
        data[i] = (uint8_t)(x >> 24);
    }
}


/* ---------------------------------------------------------------------------------------------- */
/* poly1305 - this is a toy example not suitable for crypto use at all! */
static uint8_t g_poly1305_key[32];
static uint64_t g_poly1305_counter = 0;

void cfx_poly1305_seed(unsigned seed) {
    cfx_lcg_bytes(seed, g_poly1305_key, sizeof g_poly1305_key);
    g_poly1305_counter = 0;
}

uint32_t cfx_poly1305_gen(void) {
    uint8_t tag[16];
    uint8_t msg[8];

    /* Use counter as message */
    uint64_t ctr = g_poly1305_counter++;
    for (size_t i = 0; i < 8; ++i) {
        msg[i] = (uint8_t)(ctr & 0xFFu);
        ctr >>= 8;
    }
    cfx_poly1305_mac(g_poly1305_key, msg, sizeof msg, tag);

    uint32_t r =
        ((uint32_t)tag[0]      ) |
        ((uint32_t)tag[1] <<  8) |
        ((uint32_t)tag[2] << 16) |
        ((uint32_t)tag[3] << 24);

    return r;
}

/*-------------*/
static uint8_t g_chacha20_key[32];
static uint8_t g_chacha20_nonce[12];
static uint32_t g_chacha20_counter = 0;

void cfx_chacha20_seed(uint32_t seed) {

    cfx_lcg_bytes(seed, g_chacha20_key, sizeof g_chacha20_key);
    cfx_lcg_bytes(seed, g_chacha20_nonce, sizeof g_chacha20_nonce);

    g_chacha20_counter = 0;
}

uint32_t cfx_chacha20_gen(void) {
    static uint8_t buf[64];
    static size_t idx = 64;  /* force initial refill */

    if (idx >= 64) {
        /* Fill next block */
        cfx_chacha20_block_rfc8439(g_chacha20_key, g_chacha20_counter++, g_chacha20_nonce, buf);
        idx = 0;
    }

    uint32_t v =
        ((uint32_t)buf[idx + 0]      ) |
        ((uint32_t)buf[idx + 1] <<  8) |
        ((uint32_t)buf[idx + 2] << 16) |
        ((uint32_t)buf[idx + 3] << 24);
    idx += 4;
    return v;
}

/* ---------------------------------------------------------------------------------------------- */
/* xorshift */
static uint32_t g_xorshift32_seed = 0x01;

void cfx_xorshift32_seed(unsigned seed) {
    g_xorshift32_seed = seed;
}

uint32_t cfx_xorshift32_gen(void) {
    return cfx_xorshift32(&g_xorshift32_seed);
}

uint32_t cfx_xorshift32(uint32_t* s) {
    uint32_t x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    return x;
}

static uint32_t g_xorshift32_star_seed = 0x078;


void cfx_xorshift32_star_seed(unsigned seed) {
    g_xorshift32_star_seed = seed;
}

uint32_t cfx_xorshift32_star_gen(void) {
    return cfx_xorshift32_star(&g_xorshift32_star_seed);
}

uint32_t cfx_xorshift32_star(uint32_t* s) {
    uint32_t x = *s;
    /* xorshift32 (13,17,5) */
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    /* scrambler: Knuth's 32 bit golden ratio multiplier */
    return x * 0x9E3779BBu;
}

static uint64_t g_xorshift64_seed = UINT64_C(0x057);

void cfx_xorshift64_seed(unsigned seed) {
    g_xorshift64_seed = (uint64_t)seed;
}

uint32_t cfx_xorshift64_gen(void) {
    return cfx_xorshift64(&g_xorshift64_seed);
}

uint64_t cfx_xorshift64(uint64_t* s) {
    uint64_t x = *s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *s = x;
    return x;
}

static uint64_t g_xorshift_star_seed = UINT64_C(0x1238);

void cfx_xorshift64_star_seed(unsigned seed) {
    g_xorshift_star_seed = (uint64_t)seed;
}

uint32_t cfx_xorshift64_star_gen(void) {
    return (uint32_t)cfx_xorshift64_star(&g_xorshift_star_seed);
}

uint64_t cfx_xorshift64_star(uint64_t* s) {
    uint64_t x = *s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *s = x;
    return x * UINT64_C(2685821657736338717);
}


static cfx_limb_t g_xorshift_state = (cfx_limb_t)0xABC;

void cfx_xorshift_seed(unsigned seed) {
    g_xorshift_state = (cfx_limb_t)seed;
}

cfx_limb_t cfx_xorshift_gen(void) {
    return cfx_xorshift(&g_xorshift_state);
}

cfx_limb_t cfx_xorshift(cfx_limb_t* s) {
    cfx_limb_t x = *s;
    #if (CFX_LIMB_BITS == 64)
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    #elif (CFX_LIMB_BITS == 32)
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    #else
    /* ?? */
    #endif
    *s = x;
    return x;
}



/* ---------------------------------------------------------------------------------------------- */
static uint64_t g_drand_state = UINT64_C(0x123456789ABCDEF);

void cfx_drand48_seed(unsigned seed) {
    g_drand_state = (uint64_t)seed;
}

uint32_t cfx_drand48_gen(void) {
    return cfx_drand48(&g_drand_state);
}

uint32_t cfx_drand48(uint64_t* s) {
    const uint64_t pow_2_48 = UINT64_C(281474976710656);
    *s = ((*s * 25214903917) + 11 ) % pow_2_48;
    return (uint32_t)(*s >> 16);
}

/* ---------------------------------------------------------------------------------------------- */
static uint64_t g_minstd_state = UINT64_C(0x123456789ABCDEF);

void cfx_minstd_seed(unsigned seed) {
    g_minstd_state = (uint64_t)seed;
}

uint32_t cfx_minstd_gen(void) {
    return cfx_minstd(&g_minstd_state);
}

uint32_t cfx_minstd(uint64_t* s) {
    *s *= UINT64_C(16807);
    *s %= UINT64_C(2147483647);
    return (uint32_t)*s;
}


/* ---------------------------------------------------------------------------------------------- */
/* splitmix */

static uint32_t g_splitmix32_seed = 0x12678u;

void cfx_splitmix32_seed(unsigned seed) {
    g_splitmix32_seed = seed;
}

uint32_t cfx_splitmix32_gen(void) {
    return cfx_splitmix32(&g_splitmix32_seed);
}

uint32_t cfx_splitmix32(uint32_t *s) {
    uint32_t z = (*s += 0x9E3779B9u);  /* golden ratio increment*/
    z ^= z >> 15;
    z *= 0x85EBCA6Bu;
    z ^= z >> 13;
    z *= 0xC2B2AE35u;
    z ^= z >> 16;
    return z;
}

static uint64_t g_sm64_state = UINT64_C(0x123456789ABCDEF);

void cfx_splitmix64_seed(unsigned seed) {
    g_sm64_state = (uint64_t)seed;
}

uint32_t cfx_splitmix64_gen(void) {
    return (uint32_t)cfx_splitmix64(&g_sm64_state);
}

uint64_t cfx_splitmix64(uint64_t *s) {
    uint64_t z = (*s += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

/* ---------------------------------------------------------------------------------------------- */
/* PCG "permuted congruential generator"
 *   ref: https://www.pcg-random.org/index.html
 */
static uint64_t g_pcg_state = UINT64_C(0x853c49e6748fea9b);

void cfx_pcg32_seed(unsigned seed) {
    g_pcg_state = (uint64_t)seed;
}
uint32_t cfx_pcg32_gen(void) {
    return cfx_pcg32(&g_pcg_state);
}

uint32_t cfx_pcg32(uint64_t* s) {
    const uint64_t pcg_inc = UINT64_C(0xda3e39cb94b95bdb);
    uint64_t oldstate = *s;
    *s = oldstate * UINT64_C(6364136223846793005) + (pcg_inc | 1);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}
