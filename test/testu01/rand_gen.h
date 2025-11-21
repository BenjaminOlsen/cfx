#ifndef CFX_RAND_GEN_H
#define CFX_RAND_GEN_H

#include "cfx/rand.h"
#include "cfx/poly1305.h"
#include "cfx/chacha20.h"

#include <stdint.h>
#include <stddef.h>

/* .................................................................. */
/* chacha20 */

static uint8_t g_chacha20_key[32];
static uint8_t g_chacha20_nonce[12];
static uint32_t g_chacha20_counter = 0;

static void lcg(uint32_t seed, uint8_t *data, size_t len) {
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

static void chacha20_seed(uint32_t seed) {
    
    lcg(seed, g_chacha20_key, sizeof g_chacha20_key);
    lcg(seed, g_chacha20_nonce, sizeof g_chacha20_nonce);

    g_chacha20_counter = 0;
}

static uint32_t cfx_chacha20_gen(void) {
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

/* .................................................................. */
/* poly1305 - this is a toy example not suitable for crypto use at all! */
uint8_t g_poly1305_key[32];
static uint64_t g_poly1305_counter = 0;

static void poly1305_seed(unsigned seed) {
    lcg(seed, g_poly1305_key, sizeof g_poly1305_key);
    g_poly1305_counter = 0;
}

static uint32_t cfx_poly1305_gen(void) {
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

/* .................................................................. */
/* cfx_rand */
/* seed with cfx_srand (cfx/rand.h) */
static uint32_t cfx_rand_gen(void) {
    return cfx_urand();
}

/* .................................................................. */
/* c std library rand */
/* seed with srand (stdlib.h) */
static uint32_t rand_gen(void) {
    uint32_t r = (uint32_t)(rand() & 0x7FFFFFFF);
    r ^= 1;
    r ^= (uint32_t)(rand() & 0x7FFFFFFF);
    return r;
}

/* .................................................................. */
/* Bengen */
static uint64_t g_ben_seed = 0xa5a5a5a5a5a5a5a5;


typedef uint32_t (*rand_fn)(void);
typedef void (*seed_fn)(unsigned);

typedef struct {
    const char* name;       /* name passed to TestU01 */
    rand_fn     fn;         /* RNG function */
    seed_fn     sfn;        /* fn to pass in the seed*/
} rand_desc_t;

const rand_desc_t g_rand_gens[] = {
    {"cfx_chacha20",    cfx_chacha20_gen,   chacha20_seed},
    {"cfx_poly1305",    cfx_poly1305_gen,   poly1305_seed},
    {"cfx_rand",        cfx_rand_gen,       cfx_srand},
    {"cfx_splitmix32",  cfx_splitmix32,     cfx_splitmix_seed},
    {"cfx_pcg32",       cfx_pcg32,          cfx_pcg_seed},
    {"cfx_drand48",     cfx_drand48,        cfx_drand48_seed},
    {"cfx_minstd",      cfx_minstd,         cfx_minstd_seed},
    {"rand",            rand_gen,           srand}

    /* todo later:
       { "cfx_xoshiro256ss", "--xoshiro256ss", cfx_xoshiro256ss_32 },
       ...
    */
};

static const size_t g_rand_gen_cnt = sizeof(g_rand_gens) / sizeof(g_rand_gens[0]);

#endif /* CFX_RAND_GEN_H*/
