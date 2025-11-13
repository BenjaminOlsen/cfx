#include "cfx/rand.h"
#include <stdlib.h>

cfx_limb_t cfx_rand_xorshift64(cfx_limb_t* s) {
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

void cfx_srand(unsigned seed) {  
    srand(seed);
}
