#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/numerical.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int  cfx_rand(void);                     /* returns 0..0x7fffffff */
void cfx_srand(unsigned seed);
void cfx_rand_seed_os(void);             /* (re)seed from OS RNG */
cfx_limb_t cfx_xorshift64(cfx_limb_t* s);
cfx_limb_t cfx_rand_limb(void);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_RAND_H */
