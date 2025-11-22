#ifndef CFX_RAND_H
#define CFX_RAND_H

#include "cfx/numerical.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/** the cfx_xxx_gen(void) functions rely on their corresponding _seed function to generate 
  * their sequences
*/

void        cfx_srand(unsigned seed);
int         cfx_rand(void);                     /* returns 0..0x7fffffff - uses chacha20 internally */
uint32_t    cfx_urand(void);                    /* returns 0..0xffffffff - uses chacha20 internally */
void        cfx_randombytes(void* buf, size_t len);  /* random bytes seeded by cfx_srand seed */

void        cfx_srand_os(void);                 /* uses /dev/urandom */
void        cfx_randombytes_os(void* buf, size_t len);     /* random bytes seeded by OS RNG */

cfx_limb_t  cfx_rand_limb(void);

void        cfx_lcg_seed(unsigned seed);
uint32_t    cfx_lcg_gen(void);
uint32_t    cfx_lcg(uint32_t* s);
void        cfx_lcg_bytes(uint32_t seed, uint8_t *data, size_t len);

/* ------- poly1305 ------- */
void        cfx_poly1305_seed(unsigned seed);
uint32_t    cfx_poly1305_gen(void);

/* ------- chacha20 ------- */
void        cfx_chacha20_seed(uint32_t seed);
uint32_t    cfx_chacha20_gen(void);

/* ------- xorshift ------- */
void        cfx_xorshift32_seed(unsigned seed);
uint32_t    cfx_xorshift32_gen(void);
uint32_t    cfx_xorshift32(uint32_t* s);

void        cfx_xorshift32_star_seed(unsigned seed);
uint32_t    cfx_xorshift32_star_gen(void);
uint32_t    cfx_xorshift32_star(uint32_t* s);          /* not really good */

void        cfx_xorshift64_seed(unsigned seed);
uint32_t    cfx_xorshift64_gen(void);
uint64_t    cfx_xorshift64(uint64_t* s);

void        cfx_xorshift64_star_seed(unsigned seed);
uint32_t    cfx_xorshift64_star_gen(void);
uint64_t    cfx_xorshift64_star(uint64_t* s);

void        cfx_xorshift_seed(unsigned seed);
cfx_limb_t  cfx_xorshift(cfx_limb_t* s);

/* ------- splitmix ------- */
void        cfx_splitmix32_seed(unsigned seed);
uint32_t    cfx_splitmix32_gen(void);
uint32_t    cfx_splitmix32(uint32_t* s);

void        cfx_splitmix64_seed(unsigned seed);
uint32_t    cfx_splitmix64_gen(void);
uint64_t    cfx_splitmix64(uint64_t* s);

void        cfx_pcg32_seed(unsigned seed);
uint32_t    cfx_pcg32_gen(void);
uint32_t    cfx_pcg32(uint64_t* s);

/* ------- drand48 ------- */
void        cfx_drand48_seed(unsigned seed);
uint32_t    cfx_drand48_gen(void);
uint32_t    cfx_drand48(uint64_t* s);

/* ------- minstd ------- */
void        cfx_minstd_seed(unsigned seed);
uint32_t    cfx_minstd_gen(void);
uint32_t    cfx_minstd(uint64_t* s);


#ifdef __cplusplus
}
#endif

#endif  /* CFX_RAND_H */
