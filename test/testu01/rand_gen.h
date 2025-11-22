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
static uint32_t rand_gen32(void) {
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

#define RNG_ENTRY(ID) { "cfx_" #ID, cfx_##ID##_gen32, cfx_##ID##_seed }

const rand_desc_t g_rand_gens[] = {
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
    {"cfx_rand",    cfx_urand,  cfx_srand},
    {"rand",        rand_gen32, srand}
};

static const size_t g_rand_gen_cnt = sizeof(g_rand_gens) / sizeof(g_rand_gens[0]);

#endif /* CFX_RAND_GEN_H*/
