#ifndef CFX_BITS_H
#define CFX_BITS_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

static inline uint64_t cfx_rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint32_t cfx_rotl32(uint32_t x, int k) {
    return (x << k) | (x >> (32 - k));
}

static inline uint64_t cfx_rotr64(uint64_t x, int k) {
    return (x >> k) | (x << (64 - k));
}

static inline uint32_t cfx_rotr32(uint32_t x, int k) {
    return (x >> k) | (x << (32 - k));
}

#ifdef __cplusplus
}
#endif

#endif /* CFX_BITS_H */
