#ifndef CFX_RNG_XORSHIFT_H
#define CFX_RNG_XORSHIFT_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void        cfx_xorshift32_seed(uint32_t seed);
uint32_t    cfx_xorshift32_gen32(void);
uint32_t    cfx_xorshift32(uint32_t *s);
void        cfx_xorshift32_bytes(void *buf, size_t len);

void        cfx_xorshift32star_seed(uint32_t seed);
uint32_t    cfx_xorshift32star_gen32(void);
uint32_t    cfx_xorshift32star(uint32_t *s);
void        cfx_xorshift32star_bytes(void *buf, size_t len);

void        cfx_xorshift64_seed(uint32_t seed);
uint32_t    cfx_xorshift64_gen32(void);
uint64_t    cfx_xorshift64(uint64_t *s);
void        cfx_xorshift64_bytes(void *buf, size_t len);

void        cfx_xorshift64star_seed(uint32_t seed);
uint32_t    cfx_xorshift64star_gen32(void);
uint64_t    cfx_xorshift64star(uint64_t *s);
void        cfx_xorshift64star_bytes(void *buf, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_RNG_XORSHIFT_H */
