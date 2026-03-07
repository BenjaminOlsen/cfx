#include "cfx/rng_splitmix.h"
#include "rng_internal.h"

#include <string.h>

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

DEFINE_BYTES32_FN(cfx_splitmix32,             cfx_splitmix32_gen32())
DEFINE_BYTES32_FN(cfx_splitmix64,             cfx_splitmix64_gen32())

/* seeding helpers used by other RNG families */

void cfx_rng_seed2_32(uint32_t s[2], uint32_t seed) {
    uint32_t x = seed ? seed : 0x1u;   /* avoid trivial zero */
    s[0] = cfx_splitmix32(&x);
    s[1] = cfx_splitmix32(&x);
    if (s[0] == 0u && s[1] == 0u) {
        s[0] = 1u;  /* ensure nonzero state */
    }
}

void cfx_rng_seed2_64(uint64_t s[2], uint32_t seed) {
    uint64_t x = seed ? (uint64_t)seed : UINT64_C(0x1);   /* avoid trivial zero */
    s[0] = cfx_splitmix64(&x);
    s[1] = cfx_splitmix64(&x);
    if (s[0] == 0 && s[1] == 0) {
        s[0] = UINT64_C(1);  /* ensure nonzero state */
    }
}

void cfx_rng_seed4_64(uint64_t s[4], uint32_t seed) {
    uint64_t x = seed ? (uint64_t)seed : UINT64_C(0x1);   /* avoid trivial zero */
    s[0] = cfx_splitmix64(&x);
    s[1] = cfx_splitmix64(&x);
    s[2] = cfx_splitmix64(&x);
    s[3] = cfx_splitmix64(&x);
    if (!s[0] && !s[1] && !s[2] && !s[3]) {
        s[0] = UINT64_C(1);  /* ensure nonzero state */
    }
}
