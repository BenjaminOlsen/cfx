#include "cfx/rng_pcg.h"
#include "rng_internal.h"

#include <string.h>

/** PCG "permuted congruential generator"
 *  ref: https://www.pcg-random.org/index.html
 **/
static uint64_t g_pcg_state = UINT64_C(0x853c49e6748fea9b);  /* NOT THREAD-SAFE */

void cfx_pcg32_seed(uint32_t seed) {
    const uint64_t pcg_inc = UINT64_C(0xda3e39cb94b95bdb);
    g_pcg_state = 0;
    g_pcg_state = g_pcg_state * UINT64_C(6364136223846793005) + (pcg_inc | 1);
    g_pcg_state += (uint64_t)seed;
    g_pcg_state = g_pcg_state * UINT64_C(6364136223846793005) + (pcg_inc | 1);
}
uint32_t cfx_pcg32_gen32(void) {
    return cfx_pcg32(&g_pcg_state);
}

uint32_t cfx_pcg32(uint64_t *s) {
    const uint64_t pcg_inc = UINT64_C(0xda3e39cb94b95bdb);
    uint64_t oldstate = *s;
    *s = oldstate * UINT64_C(6364136223846793005) + (pcg_inc | 1);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((32 - rot) & 31));
}

DEFINE_BYTES32_FN(cfx_pcg32,                  cfx_pcg32_gen32())
