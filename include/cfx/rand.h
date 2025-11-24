#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/types.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** the cfx_<rng>_gen32(void) functions rely on static global state to generate 
  * their sequences, written by their corresponding cfx_<rng name>_seed() function.
  * They're NOT thread safe, because of their static shared state.
  * However, the explicit state generators (with uint32_t* argument) 
  * ARE thread safe, as all state is passed in.
  **/

void        cfx_srand(uint32_t seed);
int         cfx_rand(void);                     /* returns 0..0x7fffffff - uses chacha20 internally */
uint32_t    cfx_urand(void);                    /* returns 0..0xffffffff - uses chacha20 internally */
void        cfx_randombytes(void* buf, size_t len);  /* random bytes seeded by cfx_srand seed */

void        cfx_srand_os(void);                 /* uses /dev/urandom */
void        cfx_randombytes_os(void* buf, size_t len);     /* random bytes seeded by OS RNG */

cfx_limb_t  cfx_rand_limb(void);

void        cfx_lcg_seed(uint32_t seed);
uint32_t    cfx_lcg_gen32(void);
uint32_t    cfx_lcg(uint32_t* s);
void        cfx_lcg_bytes(uint32_t seed, uint8_t *data, size_t len);

/* ------- poly1305 ------- */
void        cfx_poly1305_seed(uint32_t seed);
uint32_t    cfx_poly1305_gen32(void);

/* ------- chacha20 ------- */
void        cfx_chacha20_seed(uint32_t seed);
uint32_t    cfx_chacha20_gen32(void);

/* ------- xorshift ------- */
void        cfx_xorshift32_seed(uint32_t seed);
uint32_t    cfx_xorshift32_gen32(void);
uint32_t    cfx_xorshift32(uint32_t* s);

void        cfx_xorshift32star_seed(uint32_t seed);
uint32_t    cfx_xorshift32star_gen32(void);
uint32_t    cfx_xorshift32star(uint32_t* s);          /* not really good */

void        cfx_xorshift64_seed(uint32_t seed);
uint32_t    cfx_xorshift64_gen32(void);
uint64_t    cfx_xorshift64(uint64_t* s);

void        cfx_xorshift64star_seed(uint32_t seed);
uint32_t    cfx_xorshift64star_gen32(void);
uint64_t    cfx_xorshift64star(uint64_t* s);

void        cfx_xorshift_seed(uint32_t seed);
uint32_t    cfx_xorshift_gen32(void);
cfx_limb_t  cfx_xorshift(cfx_limb_t* s);

/* ------- splitmix ------- */
void        cfx_splitmix32_seed(uint32_t seed);
uint32_t    cfx_splitmix32_gen32(void);
uint32_t    cfx_splitmix32(uint32_t* s);

void        cfx_splitmix64_seed(uint32_t seed);
uint32_t    cfx_splitmix64_gen32(void);
uint64_t    cfx_splitmix64(uint64_t* s);

void        cfx_pcg32_seed(uint32_t seed);
uint32_t    cfx_pcg32_gen32(void);
uint32_t    cfx_pcg32(uint64_t* s);

/* ------- drand48 ------- */
void        cfx_drand48_seed(uint32_t seed);
uint32_t    cfx_drand48_gen32(void);
uint32_t    cfx_drand48(uint64_t* s);

/* ------- minstd ------- */
void        cfx_minstd_seed(uint32_t seed);
uint32_t    cfx_minstd_gen32(void);
uint32_t    cfx_minstd(uint64_t* s);

/* ------- xoroshiro family ------- */
void        cfx_xoroshiro64star_seed(uint32_t seed);
uint32_t    cfx_xoroshiro64star_gen32(void);
uint32_t    cfx_xoroshiro64star(uint32_t s[2]);

void        cfx_xoroshiro64starstar_seed(uint32_t seed);
uint32_t    cfx_xoroshiro64starstar_gen32(void);
uint32_t    cfx_xoroshiro64starstar(uint32_t s[2]);

void        cfx_xoroshiro128plus_seed(uint32_t seed);
uint32_t    cfx_xoroshiro128plus_gen32(void);
uint64_t    cfx_xoroshiro128plus(uint64_t s[2]);

void        cfx_xoroshiro128starstar_seed(uint32_t seed);
uint32_t    cfx_xoroshiro128starstar_gen32(void);
uint64_t    cfx_xoroshiro128starstar(uint64_t s[2]);

void        cfx_xoshiro256plus_seed(uint32_t seed);
uint32_t    cfx_xoshiro256plus_gen32(void);
uint64_t    cfx_xoshiro256plus(uint64_t s[4]);

void        cfx_xoshiro256starstar_seed(uint32_t seed);
uint32_t    cfx_xoshiro256starstar_gen32(void);
uint64_t    cfx_xoshiro256starstar(uint64_t s[4]);

/* ---------------------------------------------------------------------------------------------- */
typedef uint32_t (*cfx_rand_fn)(void);
typedef void (*cfx_seed_fn)(uint32_t);

typedef struct {
    const char*     name;           /* name passed to TestU01 */
    cfx_rand_fn     gen32;          /* RNG function */
    cfx_seed_fn     seed;           /* fn to pass in the seed*/
} cfx_rand_desc_t;

extern const cfx_rand_desc_t g_rand_gens[];
extern const size_t g_rand_gen_cnt;

#ifdef __cplusplus
}
#endif

#endif  /* CFX_RAND_H */
