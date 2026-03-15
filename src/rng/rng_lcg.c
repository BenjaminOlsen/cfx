#include "cfx/rng_lcg.h"
#include "rng_internal.h"

#include <string.h>

/* ---------------------------------------------------------------------------------------------- */
uint32_t g_lcg_seed = 0x126;  /* NOT THREAD-SAFE */

void cfx_lcg_seed(uint32_t seed) {
    g_lcg_seed = seed;
}

uint32_t cfx_lcg_gen32(void) {
    return cfx_lcg(&g_lcg_seed);
}

uint32_t cfx_lcg(uint32_t *s) {
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
static uint64_t g_drand_state = UINT64_C(0x123456789ABCDEF);  /* NOT THREAD-SAFE */

void cfx_drand48_seed(uint32_t seed) {
    g_drand_state = (uint64_t)seed;
}

uint32_t cfx_drand48_gen32(void) {
    return cfx_drand48(&g_drand_state);
}

uint32_t cfx_drand48(uint64_t *s) {
    *s = ((*s * UINT64_C(25214903917)) + 11) & UINT64_C(0xFFFFFFFFFFFF);
    return (uint32_t)(*s >> 16);
}

/* ---------------------------------------------------------------------------------------------- */
static uint64_t g_minstd_state = 1;  /* NOT THREAD-SAFE */

void cfx_minstd_seed(uint32_t seed) {
    g_minstd_state = seed ? (uint64_t)seed : 1;
}

uint32_t cfx_minstd_gen32(void) {
    return cfx_minstd(&g_minstd_state);
}

uint32_t cfx_minstd(uint64_t *s) {
    *s *= UINT64_C(16807);
    *s %= UINT64_C(2147483647);
    return (uint32_t)*s;
}


DEFINE_BYTES32_FN(cfx_drand48,                cfx_drand48_gen32())
DEFINE_BYTES32_FN(cfx_minstd,                 cfx_minstd_gen32())
