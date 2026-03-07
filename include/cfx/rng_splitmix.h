#ifndef CFX_RNG_SPLITMIX_H
#define CFX_RNG_SPLITMIX_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void        cfx_splitmix32_seed(uint32_t seed);
uint32_t    cfx_splitmix32_gen32(void);
uint32_t    cfx_splitmix32(uint32_t *s);
void        cfx_splitmix32_bytes(void *buf, size_t len);

void        cfx_splitmix64_seed(uint32_t seed);
uint32_t    cfx_splitmix64_gen32(void);
uint64_t    cfx_splitmix64(uint64_t *s);
void        cfx_splitmix64_bytes(void *buf, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_RNG_SPLITMIX_H */
