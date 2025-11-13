#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/num.h"

cfx_limb_t cfx_rand_xorshift64(cfx_limb_t* s);
cfx_limb_t cfx_rand_limb(void);
void cfx_srand(unsigned seed);

#endif  /* CFX_RAND_H */
