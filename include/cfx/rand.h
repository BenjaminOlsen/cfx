#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/rng_chacha20.h"
#include "cfx/rng_lcg.h"
#include "cfx/rng_xorshift.h"
#include "cfx/rng_splitmix.h"
#include "cfx/rng_pcg.h"
#include "cfx/rng_xoroshiro.h"
#include "cfx/arith.h"
#include "cfx/arch.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef int (*cfx_rng_fn)(void *ctx, uint8_t *out, size_t len);
typedef uint32_t (*cfx_prng32_fn)(void);
typedef void (*cfx_prng_seed_fn)(uint32_t);
typedef void (*cfx_prng_bytes_fn)(void *buf, size_t len);

/* limb-width xorshift (cfx-specific, depends on cfx_limb_t) */
void        cfx_xorshift_seed(uint32_t seed);
uint32_t    cfx_xorshift_gen32(void);
cfx_limb_t  cfx_xorshift(cfx_limb_t *s);
void        cfx_xorshift_bytes(void *buf, size_t len);

cfx_limb_t  cfx_rand_limb(void);

/** the cfx_<rng>_gen32(void) functions rely on static global state to generate
 * their sequences, written by their corresponding cfx_<rng name>_seed() function.
 * They're NOT thread safe, because of their static shared state.
 * However, the explicit state generators (with uint32_t* argument)
 * ARE thread safe, as all state is passed in.
 **/


typedef struct {
    const char *name;
    cfx_prng32_fn rng32;                /* 32 byte RNG function */
    cfx_prng_seed_fn seed;              /* fn to pass in the seed */
    cfx_prng_bytes_fn bytes;
} cfx_rand_desc_t;

extern const cfx_rand_desc_t g_rand_gens[];
extern const size_t g_rand_gen_cnt;

#ifdef __cplusplus
}
#endif

#endif  /* CFX_RAND_H */
