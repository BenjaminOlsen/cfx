/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

/*
 * these wrappers use hash functions in counter mode to generate
 * pseudorandom output for TestU01 statistical testing.
 */

#ifndef CFX_HASH_RNGS_H
#define CFX_HASH_RNGS_H

#include "cfx/rand.h"

#ifdef __cplusplus
extern "C" {
#endif

/* All used in counter mode */

void     cfx_sha256_ctr_seed(uint32_t seed);
uint32_t cfx_sha256_ctr_gen32(void);
void     cfx_sha256_ctr_bytes(void *buf, size_t len);

void     cfx_sha512_ctr_seed(uint32_t seed);
uint32_t cfx_sha512_ctr_gen32(void);
void     cfx_sha512_ctr_bytes(void *buf, size_t len);

void     cfx_blake2b_ctr_seed(uint32_t seed);
uint32_t cfx_blake2b_ctr_gen32(void);
void     cfx_blake2b_ctr_bytes(void *buf, size_t len);

void     cfx_blake2s_ctr_seed(uint32_t seed);
uint32_t cfx_blake2s_ctr_gen32(void);
void     cfx_blake2s_ctr_bytes(void *buf, size_t len);

extern const cfx_rand_desc_t g_hash_rngs[];
extern const size_t g_hash_rng_cnt;

#ifdef __cplusplus
}
#endif

#endif /* CFX_HASH_RNGS_H */
