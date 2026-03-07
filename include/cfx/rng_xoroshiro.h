#ifndef CFX_RNG_XOROSHIRO_H
#define CFX_RNG_XOROSHIRO_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* xoroshiro64* */
void        cfx_xoroshiro64star_seed(uint32_t seed);
uint32_t    cfx_xoroshiro64star_gen32(void);
uint32_t    cfx_xoroshiro64star(uint32_t s[2]);
void        cfx_xoroshiro64star_bytes(void *buf, size_t len);

/* xoroshiro64** */
void        cfx_xoroshiro64starstar_seed(uint32_t seed);
uint32_t    cfx_xoroshiro64starstar_gen32(void);
uint32_t    cfx_xoroshiro64starstar(uint32_t s[2]);
void        cfx_xoroshiro64starstar_bytes(void *buf, size_t len);

/* xoroshiro128+ */
void        cfx_xoroshiro128plus_seed(uint32_t seed);
uint32_t    cfx_xoroshiro128plus_gen32(void);
uint64_t    cfx_xoroshiro128plus(uint64_t s[2]);
void        cfx_xoroshiro128plus_bytes(void *buf, size_t len);

/* xoroshiro128** */
void        cfx_xoroshiro128starstar_seed(uint32_t seed);
uint32_t    cfx_xoroshiro128starstar_gen32(void);
uint64_t    cfx_xoroshiro128starstar(uint64_t s[2]);
void        cfx_xoroshiro128starstar_bytes(void *buf, size_t len);

/* xoshiro256+ */
void        cfx_xoshiro256plus_seed(uint32_t seed);
uint32_t    cfx_xoshiro256plus_gen32(void);
uint64_t    cfx_xoshiro256plus(uint64_t s[4]);
void        cfx_xoshiro256plus_bytes(void *buf, size_t len);

/* xoshiro256** */
void        cfx_xoshiro256starstar_seed(uint32_t seed);
uint32_t    cfx_xoshiro256starstar_gen32(void);
uint64_t    cfx_xoshiro256starstar(uint64_t s[4]);
void        cfx_xoshiro256starstar_bytes(void *buf, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_RNG_XOROSHIRO_H */
