#include "cfx/rand.h"
#include "cfx/chacha20.h"
#include "cfx/poly1305.h"
#include "cfx/memory.h"

#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>

#define RNG_ENTRY(ID) { "cfx_" #ID, cfx_##ID##_gen32, cfx_##ID##_seed, cfx_##ID##_bytes }

/* ---------------------------------------------------------------------------------------------- */
/* for testing's sake: c std library rand */
/* seed with srand (stdlib.h) */
// static uint32_t rand_gen32(void) {
//     uint32_t r = (uint32_t)(rand() & 0x7FFFFFFF);
//     r ^= 1;
//     r ^= (uint32_t)(rand() & 0x7FFFFFFF);
//     return r;
// }

const cfx_rand_desc_t g_rand_gens[] = {
    RNG_ENTRY(chacha20),
    RNG_ENTRY(poly1305),
    RNG_ENTRY(xorshift64),
    RNG_ENTRY(xorshift64star),
    RNG_ENTRY(xorshift32),
    RNG_ENTRY(xorshift32star),
    RNG_ENTRY(splitmix32),
    RNG_ENTRY(splitmix64),
    RNG_ENTRY(pcg32),
    RNG_ENTRY(drand48),
    RNG_ENTRY(minstd),
    RNG_ENTRY(xoroshiro64star),
    RNG_ENTRY(xoroshiro64starstar),
    RNG_ENTRY(xoroshiro128plus),
    RNG_ENTRY(xoroshiro128starstar),
    RNG_ENTRY(xoshiro256plus),
    RNG_ENTRY(xoshiro256starstar),
    {"cfx_rand", cfx_urand, cfx_srand, cfx_rand_bytes},
    // {"rand",     rand_gen32, srand,    rand_bytes
};

const size_t g_rand_gen_cnt = sizeof(g_rand_gens) / sizeof(g_rand_gens[0]);

/* ----------------------------------------------------------------------   */
/* This macro generates a byte streaming function for the given rng32
    using aligned writes for speed.

    For RNGS that produce uint32_t.

    Q: "WHY??? Why not use a function pointer argument to the rng32??""
    A: because we want to inline for speed, and we can't inline through a
       function pointer...

       DEFINE_BYTES32_FN(cfx_pcg32, cfx_pcg32_gen32)
    ---> generates void NAME##_bytes(void *buf, size_t len);
 */
#define DEFINE_BYTES32_FN(NAME, GEN32)                                    \
void NAME##_bytes(void *buf, size_t len)                                    \
{                                                                           \
    uint8_t *out = (uint8_t *)buf;                                          \
    const size_t word_align = CFX_ALIGNOF(uint32_t);                        \
                                                                            \
    /* 1) Head: align to uint32_t alignment */                              \
    while (len > 0 && ((uintptr_t)out % word_align) != 0) {                 \
        uint32_t w = (GEN32);                                               \
        *out++ = (uint8_t)w;                                                \
        len--;                                                              \
    }                                                                       \
                                                                            \
    /* 2) Main body: multiples of 4 bytes */                                \
    size_t main_len = len & ~(size_t)3;                                     \
    uint32_t *out4 = (uint32_t *)out;                                       \
                                                                            \
    while (main_len >= 4) {                                                 \
        *out4++ = (GEN32);                                                 \
        main_len -= 4;                                                      \
        len      -= 4;                                                      \
    }                                                                       \
                                                                            \
    out = (uint8_t *)out4;                                                  \
                                                                            \
    /* 3) Tail: 0..3 bytes */                                               \
    if (len > 0) {                                                          \
        uint32_t w = (GEN32);                                              \
        for (size_t i = 0; i < len; ++i) {                                  \
            out[i] = (uint8_t)(w >> (8 * i));                               \
        }                                                                   \
    }                                                                       \
}



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

void cfx_rand_bytes(void* buf, size_t len) {
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
    G.idx = 64;  /* force refill on first use */
    G.seeded = 1;
}

uint32_t cfx_rand_os(void) {
    uint32_t x = 0;
    os_getrandom(&x, sizeof x);
    return x;
}

void cfx_rand_bytes_os(void* buf, size_t len) {
    if (os_getrandom(buf, len) != 0) {
        memset(buf, 0, len);
    }
}

cfx_limb_t cfx_rand_limb(void) {
    return (cfx_limb_t)rand();
}

/* ---------------------------------------------------------------------------------------------- */
uint32_t g_lcg_seed = 0x126;

void cfx_lcg_seed(uint32_t seed) {
    g_lcg_seed = seed;
}

uint32_t cfx_lcg_gen32(void) {
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

void cfx_poly1305_seed(uint32_t seed) {
    cfx_lcg_bytes(seed, g_poly1305_key, sizeof g_poly1305_key);
    g_poly1305_counter = 0;
}

uint32_t cfx_poly1305_gen32(void) {
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
static chacha_rand_state_t R = {0};

void cfx_chacha20_seed(uint32_t seed) {
    cfx_lcg_bytes(seed, R.key, sizeof R.key);
    cfx_lcg_bytes(seed, R.nonce, sizeof R.nonce);

    R.counter = 0;
    R.seeded = 1;
}

uint32_t cfx_chacha20_gen32(void) {
    static uint8_t buf[64];
    static size_t idx = 64;  /* force initial refill */

    if (idx >= 64) {
        /* Fill next block */
        cfx_chacha20_block_rfc8439(R.key, R.counter++, R.nonce, buf);
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

#if 2
static inline void cfx_chacha20_refill_block(uint8_t out[64]) {
    cfx_chacha20_block_rfc8439(R.key, R.counter++, R.nonce, out);
    ++R.counter;
}

/* a more efficient _bytes function for chacha20  than the DEFINE_BYTES32_FN*/
void cfx_chacha20_bytes(void* buf, size_t len) {
    uint8_t *out = (uint8_t *)buf;

    if (!R.seeded) {
        return;
    }

    /* use any leftover keystream in R.buf */
    while (len && R.idx < 64) {
        *out++ = R.buf[R.idx++];
        len--;
    }

    if (len == 0) {
        return;
    }

    /* here, R.idx == 64 (buffer empty). */

    /* Fill whole 64-byte blocks directly */
    while (len >= 64) {
        cfx_chacha20_refill_block(out);   /* write directly into out */
        out += 64;
        len -= 64;
    }

    if (len == 0) {
        return;
    }

    /* Tail: generate one more block into R.buf and consume a prefix */
    cfx_chacha20_refill_block(R.buf);
    R.idx = (unsigned)len;  /* next unread byte in R.buf */
    for (unsigned i = 0; i < R.idx; ++i) {
        out[i] = R.buf[i];
    }
}
#endif

/* ---------------------------------------------------------------------------------------------- */
/* xorshift */
static uint32_t g_xorshift32_seed = 0x01;

void cfx_xorshift32_seed(uint32_t seed) {
    g_xorshift32_seed = seed;
}

uint32_t cfx_xorshift32_gen32(void) {
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

static uint32_t g_xorshift32star_seed = 0x078;


void cfx_xorshift32star_seed(uint32_t seed) {
    g_xorshift32star_seed = seed;
}

uint32_t cfx_xorshift32star_gen32(void) {
    return cfx_xorshift32star(&g_xorshift32star_seed);
}

uint32_t cfx_xorshift32star(uint32_t* s) {
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

void cfx_xorshift64_seed(uint32_t seed) {
    g_xorshift64_seed = (uint64_t)seed;
}

uint32_t cfx_xorshift64_gen32(void) {
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

void cfx_xorshift64star_seed(uint32_t seed) {
    g_xorshift_star_seed = (uint64_t)seed;
}

uint32_t cfx_xorshift64star_gen32(void) {
    return (uint32_t)cfx_xorshift64star(&g_xorshift_star_seed);
}

uint64_t cfx_xorshift64star(uint64_t* s) {
    uint64_t x = *s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *s = x;
    return x * UINT64_C(2685821657736338717);
}


static cfx_limb_t g_xorshift_state = (cfx_limb_t)0xABC;

void cfx_xorshift_seed(uint32_t seed) {
    g_xorshift_state = (cfx_limb_t)seed;
}

uint32_t cfx_xorshift_gen32(void) {
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

void cfx_drand48_seed(uint32_t seed) {
    g_drand_state = (uint64_t)seed;
}

uint32_t cfx_drand48_gen32(void) {
    return cfx_drand48(&g_drand_state);
}

uint32_t cfx_drand48(uint64_t* s) {
    const uint64_t pow_2_48 = UINT64_C(281474976710656);
    *s = ((*s * 25214903917) + 11 ) % pow_2_48;
    return (uint32_t)(*s >> 16);
}

/* ---------------------------------------------------------------------------------------------- */
static uint64_t g_minstd_state = UINT64_C(0x123456789ABCDEF);

void cfx_minstd_seed(uint32_t seed) {
    g_minstd_state = (uint64_t)seed;
}

uint32_t cfx_minstd_gen32(void) {
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

void cfx_splitmix32_seed(uint32_t seed) {
    g_splitmix32_seed = seed;
}

uint32_t cfx_splitmix32_gen32(void) {
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

void cfx_splitmix64_seed(uint32_t seed) {
    g_sm64_state = (uint64_t)seed;
}

uint32_t cfx_splitmix64_gen32(void) {
    return (uint32_t)cfx_splitmix64(&g_sm64_state);
}

uint64_t cfx_splitmix64(uint64_t *s) {
    uint64_t z = (*s += UINT64_C(0x9E3779B97F4A7C15));
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}

/* ---------------------------------------------------------------------------------------------- */
/** PCG "permuted congruential generator"
 *  ref: https://www.pcg-random.org/index.html
 **/
static uint64_t g_pcg_state = UINT64_C(0x853c49e6748fea9b);

void cfx_pcg32_seed(uint32_t seed) {
    g_pcg_state = (uint64_t)seed;
}
uint32_t cfx_pcg32_gen32(void) {
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

/* ---------------------------------------------------------------------------------------------- */
/* xoroshiro family */
static inline uint64_t rotl64(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint32_t rotl32(uint32_t x, int k) {
    return (x << k) | (x >> (32 - k));
}

static void seed2_32(uint32_t s[2], uint32_t seed) {
    uint32_t x = seed ? seed : 0x1u;   /* avoid trivial zero */
    s[0] = cfx_splitmix32(&x);
    s[1] = cfx_splitmix32(&x);
    if (s[0] == 0u && s[1] == 0u) {
        s[0] = 1u;  /* ensure nonzero state */
    }
}

static void seed2_64(uint64_t s[2], uint32_t seed) {
    uint64_t x = seed ? (uint64_t)seed : UINT64_C(0x1);   /* avoid trivial zero */
    s[0] = cfx_splitmix64(&x);
    s[1] = cfx_splitmix64(&x);
    if (s[0] == 0 && s[1] == 0) {
        s[0] = UINT64_C(1);  /* ensure nonzero state */
    }
}

static void seed4_64(uint64_t s[4], uint32_t seed) {
    uint64_t x = seed ? (uint64_t)seed : UINT64_C(0x1);   /* avoid trivial zero */
    s[0] = cfx_splitmix64(&x);
    s[1] = cfx_splitmix64(&x);
    s[2] = cfx_splitmix64(&x);
    s[3] = cfx_splitmix64(&x);
    if (!s[0] && !s[1] && !s[2] && !s[3]) {
        s[0] = UINT64_C(1);  /* ensure nonzero state */
    }
}


/* ................................. */
static uint32_t g_xoroshiro64star_state[2];

void cfx_xoroshiro64star_seed(uint32_t seed) {
    seed2_32(g_xoroshiro64star_state, seed);
}

uint32_t cfx_xoroshiro64star_gen32(void) {
    return cfx_xoroshiro64star(g_xoroshiro64star_state);
}

uint32_t cfx_xoroshiro64star(uint32_t s[2]) {
    const uint32_t s0 = s[0];
    uint32_t       s1 = s[1];

    const uint32_t result = s0 * 0x9E3779BBu;  /* * */

    s1 ^= s0;
    s[0] = rotl32(s0, 26) ^ s1 ^ (s1 << 9);
    s[1] = rotl32(s1, 13);

    return result;
}

/* ................................. */
static uint32_t g_xoroshiro64starstar_state[2];

void cfx_xoroshiro64starstar_seed(uint32_t seed) {
    seed2_32(g_xoroshiro64starstar_state, seed);
}

uint32_t cfx_xoroshiro64starstar_gen32(void) {
    return cfx_xoroshiro64starstar(g_xoroshiro64starstar_state);
}
uint32_t cfx_xoroshiro64starstar(uint32_t s[2]) {
    const uint32_t s0 = s[0];
    uint32_t       s1 = s[1];

    /* output function: ** */
    const uint32_t result = rotl32(s0 * 0x9E3779BBu, 5) * 5u;

    /* state transition */
    s1 ^= s0;
    s[0] = rotl32(s0, 26) ^ s1 ^ (s1 << 9);
    s[1] = rotl32(s1, 13);

    return result;
}

/* ................................. */
static uint64_t g_xoroshiro128plus_state[2];

void cfx_xoroshiro128plus_seed(uint32_t seed) {
    seed2_64(g_xoroshiro128plus_state, seed);
}

uint32_t cfx_xoroshiro128plus_gen32(void) {
    /* Q: better lower or higher bits? */
    return (uint32_t)(cfx_xoroshiro128plus(g_xoroshiro128plus_state) >> 32);
}

uint64_t cfx_xoroshiro128plus(uint64_t s[2]) {
    const uint64_t s0 = s[0];
    uint64_t s1 = s[1];
    const uint64_t result = s0 + s1;

    s1 ^= s0;
    s[0] = rotl64(s0, 55) ^ s1 ^ (s1 << 14); /* a, b */
    s[1] = rotl64(s1, 36);                   /* c */

    return result;
}

/* ................................. */
static uint64_t g_xoroshiro128starstar_state[2];

void cfx_xoroshiro128starstar_seed(uint32_t seed) {
    seed2_64(g_xoroshiro128starstar_state, seed);
}

uint32_t cfx_xoroshiro128starstar_gen32(void) {
    return (uint32_t)(cfx_xoroshiro128starstar(g_xoroshiro128starstar_state) >> 32);
}

uint64_t cfx_xoroshiro128starstar(uint64_t s[2]) {
    const uint64_t s0 = s[0];
    uint64_t s1 = s[1];
    const uint64_t result = rotl64(s0 * 5, 7) * 9;

    s1 ^= s0;
    s[0] = rotl64(s0, 55) ^ s1 ^ (s1 << 14);
    s[1] = rotl64(s1, 36);

    return result;
}

/* ................................. */
static uint64_t g_xoshiro256plus_state[4];

void cfx_xoshiro256plus_seed(uint32_t seed) {
    seed4_64(g_xoshiro256plus_state, seed);
}

uint32_t cfx_xoshiro256plus_gen32(void) {
    return (uint32_t)(cfx_xoshiro256plus(g_xoshiro256plus_state) >> 32);
}

uint64_t cfx_xoshiro256plus(uint64_t s[4]) {
    const uint64_t result = s[0] + s[3];
    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = rotl64(s[3], 45);

    return result;
}

/* ................................. */
static uint64_t g_xoshiro256starstar[4];

void cfx_xoshiro256starstar_seed(uint32_t seed) {
    seed4_64(g_xoshiro256starstar, seed);
}

uint32_t cfx_xoshiro256starstar_gen32(void) {
    return (uint32_t)(cfx_xoshiro256starstar(g_xoshiro256starstar) >> 32);
}

uint64_t cfx_xoshiro256starstar(uint64_t s[4]) {
    const uint64_t result = rotl64(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = rotl64(s[3], 45);

    return result;
}

/* DEFINE_BYTES32_FN(cfx_chacha20,               cfx_chacha20_gen32()) -- better version above*/
DEFINE_BYTES32_FN(cfx_poly1305,               cfx_poly1305_gen32())
DEFINE_BYTES32_FN(cfx_xorshift32,             cfx_xorshift32_gen32())
DEFINE_BYTES32_FN(cfx_xorshift32star,         cfx_xorshift32star_gen32())
DEFINE_BYTES32_FN(cfx_xorshift64,             cfx_xorshift64_gen32())
DEFINE_BYTES32_FN(cfx_xorshift64star,         cfx_xorshift64star_gen32())
DEFINE_BYTES32_FN(cfx_xorshift,               cfx_xorshift_gen32())
DEFINE_BYTES32_FN(cfx_splitmix32,             cfx_splitmix32_gen32())
DEFINE_BYTES32_FN(cfx_splitmix64,             cfx_splitmix64_gen32())
DEFINE_BYTES32_FN(cfx_pcg32,                  cfx_pcg32_gen32())
DEFINE_BYTES32_FN(cfx_drand48,                cfx_drand48_gen32())
DEFINE_BYTES32_FN(cfx_minstd,                 cfx_minstd_gen32())
DEFINE_BYTES32_FN(cfx_xoroshiro64star,        cfx_xoroshiro64star_gen32())
DEFINE_BYTES32_FN(cfx_xoroshiro64starstar,    cfx_xoroshiro64starstar_gen32())
DEFINE_BYTES32_FN(cfx_xoroshiro128plus,       cfx_xoroshiro128plus_gen32())
DEFINE_BYTES32_FN(cfx_xoroshiro128starstar,   cfx_xoroshiro128starstar_gen32())
DEFINE_BYTES32_FN(cfx_xoshiro256plus,         cfx_xoshiro256plus_gen32())
DEFINE_BYTES32_FN(cfx_xoshiro256starstar,     cfx_xoshiro256starstar_gen32())
