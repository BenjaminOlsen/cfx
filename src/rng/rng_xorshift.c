#include "cfx/rand.h"
#include "rng_internal.h"

#include <string.h>

/* xorshift */
static uint32_t g_xorshift32_seed = 0x01;

void cfx_xorshift32_seed(uint32_t seed) {
    g_xorshift32_seed = seed ? seed : 1u;
}

uint32_t cfx_xorshift32_gen32(void) {
    return cfx_xorshift32(&g_xorshift32_seed);
}

uint32_t cfx_xorshift32(uint32_t *s) {
    uint32_t x = *s;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *s = x;
    return x;
}

static uint32_t g_xorshift32star_seed = 0x078;


void cfx_xorshift32star_seed(uint32_t seed) {
    g_xorshift32star_seed = seed ? seed : 1u;
}

uint32_t cfx_xorshift32star_gen32(void) {
    return cfx_xorshift32star(&g_xorshift32star_seed);
}

uint32_t cfx_xorshift32star(uint32_t *s) {
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
    g_xorshift64_seed = seed ? (uint64_t)seed : UINT64_C(1);
}

uint32_t cfx_xorshift64_gen32(void) {
    return cfx_xorshift64(&g_xorshift64_seed);
}

uint64_t cfx_xorshift64(uint64_t *s) {
    uint64_t x = *s;
    x ^= x << 13;
    x ^= x >> 7;
    x ^= x << 17;
    *s = x;
    return x;
}

static uint64_t g_xorshift_star_seed = UINT64_C(0x1238);

void cfx_xorshift64star_seed(uint32_t seed) {
    g_xorshift_star_seed = seed ? (uint64_t)seed : UINT64_C(1);
}

uint32_t cfx_xorshift64star_gen32(void) {
    return (uint32_t)cfx_xorshift64star(&g_xorshift_star_seed);
}

uint64_t cfx_xorshift64star(uint64_t *s) {
    uint64_t x = *s;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *s = x;
    return x * UINT64_C(2685821657736338717);
}


static cfx_limb_t g_xorshift_state = (cfx_limb_t)0xABC;

void cfx_xorshift_seed(uint32_t seed) {
    g_xorshift_state = seed ? (cfx_limb_t)seed : (cfx_limb_t)1;
}

uint32_t cfx_xorshift_gen32(void) {
    return cfx_xorshift(&g_xorshift_state);
}

cfx_limb_t cfx_xorshift(cfx_limb_t *s) {
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


DEFINE_BYTES32_FN(cfx_xorshift32,             cfx_xorshift32_gen32())
DEFINE_BYTES32_FN(cfx_xorshift32star,         cfx_xorshift32star_gen32())
DEFINE_BYTES32_FN(cfx_xorshift64,             cfx_xorshift64_gen32())
DEFINE_BYTES32_FN(cfx_xorshift64star,         cfx_xorshift64star_gen32())
DEFINE_BYTES32_FN(cfx_xorshift,               cfx_xorshift_gen32())
