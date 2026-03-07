#ifndef CFX_RNG_LCG_H
#define CFX_RNG_LCG_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* LCG (Numerical Recipes) */
void        cfx_lcg_seed(uint32_t seed);
uint32_t    cfx_lcg_gen32(void);
uint32_t    cfx_lcg(uint32_t *s);
void        cfx_lcg_bytes(uint32_t seed, uint8_t *data, size_t len);

/* drand48 LCG */
void        cfx_drand48_seed(uint32_t seed);
uint32_t    cfx_drand48_gen32(void);
uint32_t    cfx_drand48(uint64_t *s);
void        cfx_drand48_bytes(void *buf, size_t len);

/* MINSTD (Park-Miller) */
void        cfx_minstd_seed(uint32_t seed);
uint32_t    cfx_minstd_gen32(void);
uint32_t    cfx_minstd(uint64_t *s);
void        cfx_minstd_bytes(void *buf, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_RNG_LCG_H */
