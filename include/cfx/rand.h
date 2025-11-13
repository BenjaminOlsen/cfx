#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/num.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

int  cfx_rand(void);                     /* returns 0..0x7fffffff */
void cfx_srand(unsigned seed);
void cfx_rand_seed_os(void);             /* (re)seed from OS RNG */
cfx_limb_t cfx_rand_xorshift64(cfx_limb_t* s);
cfx_limb_t cfx_rand_limb(void);

void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter,const uint8_t nonce[12], uint8_t out[64]);
void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], const uint8_t *pt, size_t pt_len, uint8_t *ct);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_RAND_H */
