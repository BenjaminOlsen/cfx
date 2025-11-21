#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/numerical.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

uint32_t cfx_urand(void);               /* returns 0..0xffffffff */
int  cfx_rand(void);                    /* returns 0..0x7fffffff */
void cfx_srand(unsigned seed);
void cfx_randombytes(void* buf, size_t len);  /* randon bytes seeded by cfx_srand seed */

void cfx_splitmix_seed(unsigned seed);
uint32_t cfx_splitmix32(void);

void cfx_pcg_seed(unsigned seed);
uint32_t cfx_pcg32(void);

void cfx_rand_seed_os(void);             /* (re)seed from OS RNG */
void cfx_randombytes_os(void* buf, size_t len);  /* random bytes seeded by OS RNG */

cfx_limb_t cfx_xorshift64(cfx_limb_t* s);
cfx_limb_t cfx_rand_limb(void);

void cfx_drand48_seed(unsigned seed);
uint32_t cfx_drand48(void);

void cfx_minstd_seed(unsigned seed);
uint32_t cfx_minstd(void);


#ifdef __cplusplus
}
#endif

#endif  /* CFX_RAND_H */
