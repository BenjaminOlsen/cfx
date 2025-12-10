#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/arith.h"
#include "cfx/arch.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif


typedef int         (*cfx_rng_fn)(void* ctx, uint8_t* out, size_t len);
typedef uint32_t    (*cfx_prng32_fn)(void);
typedef void        (*cfx_prng_seed_fn)(uint32_t);
typedef void        (*cfx_prng_bytes_fn)(void *buf, size_t len);

#if CFX_HAVE_AVX2
#define CFX_CHACHA_RNG_CTX_SIZE 736
#elif CFX_SIMD
#define CFX_CHACHA_RNG_CTX_SIZE 480
#else
#define CFX_CHACHA_RNG_CTX_SIZE 224
#endif

typedef union {
    uint8_t opaque[CFX_CHACHA_RNG_CTX_SIZE];
    CFX_ALIGNAS(32) uint64_t _force_align;
} cfx_rng_ctx_t;


void        cfx_chacha20_rng_init(cfx_rng_ctx_t* st, uint32_t seed);
void        cfx_chacha20_rng_bytes(cfx_rng_ctx_t* st, void* buf, size_t len);
int         cfx_chacha20_rng(void* ctx, uint8_t* out, size_t len);

void        cfx_srand(uint32_t seed);
int         cfx_rand(void);                         /* returns 0..0x7fffffff - uses chacha20 internally */
uint32_t    cfx_urand(void);                        /* returns 0..0xffffffff - uses chacha20 internally */
void        cfx_rand_bytes(void* buf, size_t len);  /* random bytes seeded by cfx_srand seed */
int         cfx_rng(void* ctx, uint8_t* out, size_t len);

void        cfx_srand_os(void);                     /* uses /dev/urandom */
uint32_t    cfx_rand_os(void);
void        cfx_rand_bytes_os(void* buf, size_t len);     /* random bytes seeded by OS RNG */

cfx_limb_t  cfx_rand_limb(void);

/** the cfx_<rng>_gen32(void) functions rely on static global state to generate
  * their sequences, written by their corresponding cfx_<rng name>_seed() function.
  * They're NOT thread safe, because of their static shared state.
  * However, the explicit state generators (with uint32_t* argument)
  * ARE thread safe, as all state is passed in.
  **/

/* -------- These are toy functions below -------- */
void        cfx_lcg_seed(uint32_t seed);
uint32_t    cfx_lcg_gen32(void);
uint32_t    cfx_lcg(uint32_t* s);
void        cfx_lcg_bytes(uint32_t seed, uint8_t *data, size_t len);

/* ------- poly1305 ------- */
void        cfx_poly1305_seed(uint32_t seed);
uint32_t    cfx_poly1305_gen32(void);
void        cfx_poly1305_bytes(void* buf, size_t len);

/* ------- chacha20 ------- */
void        cfx_chacha20_seed(uint32_t seed);
uint32_t    cfx_chacha20_gen32(void);
void        cfx_chacha20_bytes(void* buf, size_t len);

/* ------- xorshift ------- */
void        cfx_xorshift32_seed(uint32_t seed);
uint32_t    cfx_xorshift32_gen32(void);
uint32_t    cfx_xorshift32(uint32_t* s);
void        cfx_xorshift32_bytes(void* buf, size_t len);

void        cfx_xorshift32star_seed(uint32_t seed);
uint32_t    cfx_xorshift32star_gen32(void);
uint32_t    cfx_xorshift32star(uint32_t* s);          /* not really good */
void        cfx_xorshift32star_bytes(void* buf, size_t len);

void        cfx_xorshift64_seed(uint32_t seed);
uint32_t    cfx_xorshift64_gen32(void);
uint64_t    cfx_xorshift64(uint64_t* s);
void        cfx_xorshift64_bytes(void* buf, size_t len);

void        cfx_xorshift64star_seed(uint32_t seed);
uint32_t    cfx_xorshift64star_gen32(void);
uint64_t    cfx_xorshift64star(uint64_t* s);
void        cfx_xorshift64star_bytes(void* buf, size_t len);

void        cfx_xorshift_seed(uint32_t seed);
uint32_t    cfx_xorshift_gen32(void);
cfx_limb_t  cfx_xorshift(cfx_limb_t* s);
void        cfx_xorshift_bytes(void* buf, size_t len);

/* ------- splitmix ------- */
void        cfx_splitmix32_seed(uint32_t seed);
uint32_t    cfx_splitmix32_gen32(void);
uint32_t    cfx_splitmix32(uint32_t* s);
void        cfx_splitmix32_bytes(void* buf, size_t len);

void        cfx_splitmix64_seed(uint32_t seed);
uint32_t    cfx_splitmix64_gen32(void);
uint64_t    cfx_splitmix64(uint64_t* s);
void        cfx_splitmix64_bytes(void* buf, size_t len);

void        cfx_pcg32_seed(uint32_t seed);
uint32_t    cfx_pcg32_gen32(void);
uint32_t    cfx_pcg32(uint64_t* s);
void        cfx_pcg32_bytes(void* buf, size_t len);

/* ------- drand48 ------- */
void        cfx_drand48_seed(uint32_t seed);
uint32_t    cfx_drand48_gen32(void);
uint32_t    cfx_drand48(uint64_t* s);
void        cfx_drand48_bytes(void* buf, size_t len);

/* ------- minstd ------- */
void        cfx_minstd_seed(uint32_t seed);
uint32_t    cfx_minstd_gen32(void);
uint32_t    cfx_minstd(uint64_t* s);
void        cfx_minstd_bytes(void* buf, size_t len);

/* ------- xoroshiro family ------- */
void        cfx_xoroshiro64star_seed(uint32_t seed);
uint32_t    cfx_xoroshiro64star_gen32(void);
uint32_t    cfx_xoroshiro64star(uint32_t s[2]);
void        cfx_xoroshiro64star_bytes(void* buf, size_t len);

void        cfx_xoroshiro64starstar_seed(uint32_t seed);
uint32_t    cfx_xoroshiro64starstar_gen32(void);
uint32_t    cfx_xoroshiro64starstar(uint32_t s[2]);
void        cfx_xoroshiro64starstar_bytes(void* buf, size_t len);

void        cfx_xoroshiro128plus_seed(uint32_t seed);
uint32_t    cfx_xoroshiro128plus_gen32(void);
uint64_t    cfx_xoroshiro128plus(uint64_t s[2]);
void        cfx_xoroshiro128plus_bytes(void* buf, size_t len);

void        cfx_xoroshiro128starstar_seed(uint32_t seed);
uint32_t    cfx_xoroshiro128starstar_gen32(void);
uint64_t    cfx_xoroshiro128starstar(uint64_t s[2]);
void        cfx_xoroshiro128starstar_bytes(void* buf, size_t len);

void        cfx_xoshiro256plus_seed(uint32_t seed);
uint32_t    cfx_xoshiro256plus_gen32(void);
uint64_t    cfx_xoshiro256plus(uint64_t s[4]);
void        cfx_xoshiro256plus_bytes(void* buf, size_t len);

void        cfx_xoshiro256starstar_seed(uint32_t seed);
uint32_t    cfx_xoshiro256starstar_gen32(void);
uint64_t    cfx_xoshiro256starstar(uint64_t s[4]);
void        cfx_xoshiro256starstar_bytes(void* buf, size_t len);

/* ---------------------------------------------------------------------------------------------- */


typedef struct {
    const char*         name;
    cfx_prng32_fn       rng32;          /* 32 byte RNG function */
    cfx_prng_seed_fn    seed;           /* fn to pass in the seed */
    cfx_prng_bytes_fn   bytes;
} cfx_rand_desc_t;

extern const cfx_rand_desc_t g_rand_gens[];
extern const size_t g_rand_gen_cnt;

#ifdef __cplusplus
}
#endif

#endif  /* CFX_RAND_H */
