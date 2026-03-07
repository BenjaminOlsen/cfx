#ifndef CFX_RNG_INTERNAL_H
#define CFX_RNG_INTERNAL_H

#include <stdint.h>
#include <stddef.h>
#include "cfx/macros.h"
#include "cfx/arch.h"

/* rotate helpers */
static inline uint64_t cfx_rotl64(uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

static inline uint32_t cfx_rotl32(uint32_t x, int k) {
    return (x << k) | (x >> (32 - k));
}

/* seeding helpers (implemented in rng_splitmix.c) */
void cfx_rng_seed2_32(uint32_t s[2], uint32_t seed);
void cfx_rng_seed2_64(uint64_t s[2], uint32_t seed);
void cfx_rng_seed4_64(uint64_t s[4], uint32_t seed);

/*
 * Generates a _bytes streaming function from a gen32 function.
 * Uses aligned writes for speed.
 *
 * The gen32 call is inlined because it appears directly in the macro
 * expansion (not behind a function pointer).
 *
 * Usage: DEFINE_BYTES32_FN(cfx_pcg32, cfx_pcg32_gen32())
 *  --->  void cfx_pcg32_bytes(void *buf, size_t len)
 */
#define DEFINE_BYTES32_FN(NAME, GEN32)                                      \
    void NAME ## _bytes(void *buf, size_t len)                              \
    {                                                                       \
        uint8_t *out = (uint8_t *)buf;                                      \
        const size_t word_size = sizeof(uint32_t);                          \
        const size_t word_align = CFX_ALIGNOF(uint32_t);                    \
                                                                            \
        if (len > 0 && ((uintptr_t)out % word_align) != 0) {               \
            uint32_t w = (GEN32);                                           \
            while (len > 0 && ((uintptr_t)out % word_align) != 0) {        \
                *out++ = (uint8_t)w;                                        \
                w >>= 8;                                                    \
                len--;                                                      \
            }                                                               \
        }                                                                   \
                                                                            \
        while (len >= word_size) {                                          \
            *(uint32_t *)out = (GEN32);                                     \
            out += word_size;                                               \
            len -= word_size;                                               \
        }                                                                   \
                                                                            \
        if (len > 0) {                                                      \
            uint32_t w = (GEN32);                                           \
            for (size_t i = 0; i < len; ++i) {                              \
                out[i] = (uint8_t)(w >> (8 * i));                           \
            }                                                               \
        }                                                                   \
    }

#endif /* CFX_RNG_INTERNAL_H */
