#include "cfx/rand.h"

#define RNG_ENTRY(ID) { "cfx_" #ID, cfx_ ## ID ## _gen32, cfx_ ## ID ## _seed, cfx_ ## ID ## _bytes }

const cfx_rand_desc_t g_rand_gens[] = {
    RNG_ENTRY(chacha20),
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
};

const size_t g_rand_gen_cnt = sizeof(g_rand_gens) / sizeof(g_rand_gens[0]);
