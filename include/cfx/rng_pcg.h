#ifndef CFX_RNG_PCG_H
#define CFX_RNG_PCG_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void        cfx_pcg32_seed(uint32_t seed);
uint32_t    cfx_pcg32_gen32(void);
uint32_t    cfx_pcg32(uint64_t *s);
void        cfx_pcg32_bytes(void *buf, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_RNG_PCG_H */
