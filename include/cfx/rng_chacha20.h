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
#define CFX_CHACHA_RNG_CTX_SIZE  1152  /* buf(512) + ctx(512) + fields */
#define CFX_CHACHA_RNG_CTX_ALIGN 32
#elif CFX_SIMD
#define CFX_CHACHA_RNG_CTX_SIZE  480
#define CFX_CHACHA_RNG_CTX_ALIGN CFX_CTX_ALIGN
#else
#define CFX_CHACHA_RNG_CTX_SIZE  224
#define CFX_CHACHA_RNG_CTX_ALIGN CFX_CTX_ALIGN
#endif

typedef union {
    CFX_ALIGNAS(CFX_CHACHA_RNG_CTX_ALIGN) uint8_t opaque[CFX_CHACHA_RNG_CTX_SIZE];
} cfx_rng_ctx_t;

void        cfx_chacha20_rng_init(cfx_rng_ctx_t *st, uint32_t seed);
int         cfx_chacha20_rng(void *ctx, uint8_t *out, size_t len);

/* NOT THREAD-SAFE: these use global state. Use cfx_chacha20_rng_init/_rng for thread safety. */
void        cfx_chacha20_seed(uint32_t seed);
uint32_t    cfx_chacha20_gen32(void);
void        cfx_chacha20_bytes(void *buf, size_t len);

/* NOT THREAD-SAFE: these use global state. Use _ctx variants for thread safety. */
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
