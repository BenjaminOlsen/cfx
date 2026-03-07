#ifndef CFX_RNG_CHACHA20_H
#define CFX_RNG_CHACHA20_H

#include <stdint.h>
#include <stddef.h>
#include "cfx/arch.h"
#include "cfx/macros.h"

#ifdef __cplusplus
extern "C" {
#endif

#if CFX_HAVE_AVX2
#define CFX_CHACHA_RNG_CTX_SIZE 1152  /* buf(512) + ctx(512) + fields */
#elif CFX_SIMD
#define CFX_CHACHA_RNG_CTX_SIZE 480
#else
#define CFX_CHACHA_RNG_CTX_SIZE 224
#endif

typedef union {
    uint8_t opaque[CFX_CHACHA_RNG_CTX_SIZE];
    CFX_ALIGNAS(32) uint64_t _force_align;
} cfx_rng_ctx_t;

void        cfx_chacha20_rng_init(cfx_rng_ctx_t *st, uint32_t seed);
int         cfx_chacha20_rng(void *ctx, uint8_t *out, size_t len);

void        cfx_chacha20_seed(uint32_t seed);
uint32_t    cfx_chacha20_gen32(void);
void        cfx_chacha20_bytes(void *buf, size_t len);

/* cfx's chosen rand - wraps chacha20 */
void        cfx_srand(uint32_t seed);
int         cfx_rand(void);
uint32_t    cfx_urand(void);
void        cfx_rand_bytes(void *buf, size_t len);
int         cfx_rng(void *ctx, uint8_t *out, size_t len);

/* OS entropy */
void        cfx_srand_os(void);
uint32_t    cfx_rand_os(void);
void        cfx_rand_bytes_os(void *buf, size_t len);

#ifdef __cplusplus
}
#endif

#endif /* CFX_RNG_CHACHA20_H */
