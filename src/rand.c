#include "cfx/rand.h"
#include "cfx/chacha20.h"
#include "cfx/memory.h"

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

/* ---------------------------------------------------------------------------------------------- */
/* splitmix */
static uint64_t splitmix64_next(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

static uint64_t g_sm64_state = 0x123456789ABCDEFULL;

void cfx_splitmix_seed(unsigned seed) {   
    g_sm64_state = (uint64_t)seed;
}

uint32_t cfx_splitmix32(void) {
    splitmix64_next(&g_sm64_state);
    return (uint32_t)(g_sm64_state >> 32);
}

/* ---------------------------------------------------------------------------------------------- */
/* PCG "permuted congruential generator"
 *   ref: https://www.pcg-random.org/index.html 
*/
static uint64_t g_pcg_state = 0x853c49e6748fea9bull;

void cfx_pcg_seed(unsigned seed) {
    g_pcg_state = (uint64_t)seed;
}

uint32_t cfx_pcg32(void) {
    const uint64_t pcg_inc = 0xda3e39cb94b95bdbull;
    uint64_t oldstate = g_pcg_state;
    g_pcg_state = oldstate * 6364136223846793005ULL + (pcg_inc | 1);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

void cfx_srand(unsigned int seed) {
    uint64_t s = ((uint64_t)seed << 32) | 0xA5A5A5A5u;
    /* Fill key, nonce, counter from SplitMix64 stream */
    for (int i = 0; i < 4; ++i) {
        uint64_t x = splitmix64_next(&s);
        CFX_STORE64_LE(G.key + 8*i, x);
    }
    uint64_t n0 = splitmix64_next(&s);
    memcpy(G.nonce, &n0, 12);   /* first 12 bytes */
    G.counter = splitmix64_next(&s);
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
    /* Pull 4 bytes, assemble LE uint32, mask to 31 bits */
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
        p[i]= G.buf[G.idx++];
    }
}

void cfx_randombytes_os(void* buf, size_t len) {
    if (os_getrandom(buf, len) != 0) {
        memset(buf, 0, len);
    }
}
