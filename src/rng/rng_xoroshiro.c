#include "cfx/rng_xoroshiro.h"
#include "rng_internal.h"

#include <string.h>

/* xoroshiro family */

/* ................................. */
static uint32_t g_xoroshiro64star_state[2];

void cfx_xoroshiro64star_seed(uint32_t seed) {
    cfx_rng_seed2_32(g_xoroshiro64star_state, seed);
}

uint32_t cfx_xoroshiro64star_gen32(void) {
    return cfx_xoroshiro64star(g_xoroshiro64star_state);
}

uint32_t cfx_xoroshiro64star(uint32_t s[2]) {
    const uint32_t s0 = s[0];
    uint32_t s1 = s[1];

    const uint32_t result = s0 * 0x9E3779BBu;  /* * */

    s1 ^= s0;
    s[0] = cfx_rotl32(s0, 26) ^ s1 ^ (s1 << 9);
    s[1] = cfx_rotl32(s1, 13);

    return result;
}

/* ................................. */
static uint32_t g_xoroshiro64starstar_state[2];

void cfx_xoroshiro64starstar_seed(uint32_t seed) {
    cfx_rng_seed2_32(g_xoroshiro64starstar_state, seed);
}

uint32_t cfx_xoroshiro64starstar_gen32(void) {
    return cfx_xoroshiro64starstar(g_xoroshiro64starstar_state);
}
uint32_t cfx_xoroshiro64starstar(uint32_t s[2]) {
    const uint32_t s0 = s[0];
    uint32_t s1 = s[1];

    /* output function: */
    const uint32_t result = cfx_rotl32(s0 * 0x9E3779BBu, 5) * 5u;

    /* state transition */
    s1 ^= s0;
    s[0] = cfx_rotl32(s0, 26) ^ s1 ^ (s1 << 9);
    s[1] = cfx_rotl32(s1, 13);

    return result;
}

/* ................................. */
static uint64_t g_xoroshiro128plus_state[2];

void cfx_xoroshiro128plus_seed(uint32_t seed) {
    cfx_rng_seed2_64(g_xoroshiro128plus_state, seed);
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
    s[0] = cfx_rotl64(s0, 55) ^ s1 ^ (s1 << 14); /* a, b */
    s[1] = cfx_rotl64(s1, 36);                   /* c */

    return result;
}

/* ................................. */
static uint64_t g_xoroshiro128starstar_state[2];

void cfx_xoroshiro128starstar_seed(uint32_t seed) {
    cfx_rng_seed2_64(g_xoroshiro128starstar_state, seed);
}

uint32_t cfx_xoroshiro128starstar_gen32(void) {
    return (uint32_t)(cfx_xoroshiro128starstar(g_xoroshiro128starstar_state) >> 32);
}

uint64_t cfx_xoroshiro128starstar(uint64_t s[2]) {
    const uint64_t s0 = s[0];
    uint64_t s1 = s[1];
    const uint64_t result = cfx_rotl64(s0 * 5, 7) * 9;

    s1 ^= s0;
    s[0] = cfx_rotl64(s0, 55) ^ s1 ^ (s1 << 14);
    s[1] = cfx_rotl64(s1, 36);

    return result;
}

/* ................................. */
static uint64_t g_xoshiro256plus_state[4];

void cfx_xoshiro256plus_seed(uint32_t seed) {
    cfx_rng_seed4_64(g_xoshiro256plus_state, seed);
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
    s[3] = cfx_rotl64(s[3], 45);

    return result;
}

/* ................................. */
static uint64_t g_xoshiro256starstar[4];

void cfx_xoshiro256starstar_seed(uint32_t seed) {
    cfx_rng_seed4_64(g_xoshiro256starstar, seed);
}

uint32_t cfx_xoshiro256starstar_gen32(void) {
    return (uint32_t)(cfx_xoshiro256starstar(g_xoshiro256starstar) >> 32);
}

uint64_t cfx_xoshiro256starstar(uint64_t s[4]) {
    const uint64_t result = cfx_rotl64(s[1] * 5, 7) * 9;
    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;
    s[3] = cfx_rotl64(s[3], 45);

    return result;
}

DEFINE_BYTES32_FN(cfx_xoroshiro64star,        cfx_xoroshiro64star_gen32())
DEFINE_BYTES32_FN(cfx_xoroshiro64starstar,    cfx_xoroshiro64starstar_gen32())
DEFINE_BYTES32_FN(cfx_xoroshiro128plus,       cfx_xoroshiro128plus_gen32())
DEFINE_BYTES32_FN(cfx_xoroshiro128starstar,   cfx_xoroshiro128starstar_gen32())
DEFINE_BYTES32_FN(cfx_xoshiro256plus,         cfx_xoshiro256plus_gen32())
DEFINE_BYTES32_FN(cfx_xoshiro256starstar,     cfx_xoshiro256starstar_gen32())
