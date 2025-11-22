#ifndef CFX_RAND_GEN_H
#define CFX_RAND_GEN_H

#include "cfx/rand.h"
#include "cfx/poly1305.h"
#include "cfx/chacha20.h"

#include <stdint.h>
#include <stddef.h>


/* .................................................................. */
/* for testing's sake: c std library rand */
/* seed with srand (stdlib.h) */
static uint32_t rand_gen(void) {
    uint32_t r = (uint32_t)(rand() & 0x7FFFFFFF);
    r ^= 1;
    r ^= (uint32_t)(rand() & 0x7FFFFFFF);
    return r;
}


typedef uint32_t (*rand_fn)(void);
typedef void (*seed_fn)(unsigned);

typedef struct {
    const char* name;       /* name passed to TestU01 */
    rand_fn     fn;         /* RNG function */
    seed_fn     sfn;        /* fn to pass in the seed*/
} rand_desc_t;

const rand_desc_t g_rand_gens[] = {
    {"cfx_rand",                cfx_urand,                  cfx_srand},
    {"cfx_chacha20",            cfx_chacha20_gen,           cfx_chacha20_seed},
    {"cfx_poly1305",            cfx_poly1305_gen,           cfx_poly1305_seed},
    {"cfx_splitmix32",          cfx_splitmix32_gen,         cfx_splitmix32_seed},
    {"cfx_splitmix64",          cfx_splitmix64_gen,         cfx_splitmix64_seed},
    {"cfx_xorshift64",          cfx_xorshift64_gen,         cfx_xorshift64_seed},
    {"cfx_xorshift64_star",     cfx_xorshift64_star_gen,    cfx_xorshift64_seed},
    {"cfx_xorshift32",          cfx_xorshift32_gen,         cfx_xorshift32_seed},
    {"cfx_xorshift32_star",     cfx_xorshift32_star_gen,    cfx_xorshift32_star_seed},
    {"cfx_splitmix32",          cfx_splitmix32_gen,         cfx_splitmix32_seed},
    {"cfx_splitmix64",          cfx_splitmix64_gen,         cfx_splitmix64_seed},
    {"cfx_pcg32",               cfx_pcg32_gen,              cfx_pcg32_seed},
    {"cfx_drand48",             cfx_drand48_gen,            cfx_drand48_seed},
    {"cfx_minstd",              cfx_minstd_gen,             cfx_minstd_seed},
    {"rand",                    rand_gen,                   srand}

    /* todo later:
       { "cfx_xoshiro256ss", "--xoshiro256ss", cfx_xoshiro256ss_32 },
       ...
    */
};

static const size_t g_rand_gen_cnt = sizeof(g_rand_gens) / sizeof(g_rand_gens[0]);

#endif /* CFX_RAND_GEN_H*/
