/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "hash_rngs.h"

#include "cfx/sha256.h"
#include "cfx/sha512.h"
#include "cfx/blake2.h"
#include "cfx/memory.h"

#include <string.h>


/* SHA-256 */
static struct {
    uint64_t counter;
    uint8_t buffer[32];
    size_t pos;
} sha256_ctr_state;

void cfx_sha256_ctr_seed(uint32_t seed) {
    sha256_ctr_state.counter = seed;
    sha256_ctr_state.pos = 32;
}

static void sha256_ctr_fill(void) {
    cfx_sha256_ctx ctx;
    uint8_t input[8];

    cfx_store64_le(input, sha256_ctr_state.counter);
    sha256_ctr_state.counter++;

    cfx_sha256_init(&ctx);
    cfx_sha256_update(&ctx, input, 8);
    cfx_sha256_final(&ctx, sha256_ctr_state.buffer);
    sha256_ctr_state.pos = 0;
}

uint32_t cfx_sha256_ctr_gen32(void) {
    if (sha256_ctr_state.pos >= 32) {
        sha256_ctr_fill();
    }
    uint32_t r = cfx_load32_le(&sha256_ctr_state.buffer[sha256_ctr_state.pos]);
    sha256_ctr_state.pos += 4;
    return r;
}

void cfx_sha256_ctr_bytes(void *buf, size_t len) {
    uint8_t *p = (uint8_t *)buf;
    while (len > 0) {
        if (sha256_ctr_state.pos >= 32) {
            sha256_ctr_fill();
        }
        size_t avail = 32 - sha256_ctr_state.pos;
        size_t take = (len < avail) ? len : avail;
        memcpy(p, &sha256_ctr_state.buffer[sha256_ctr_state.pos], take);
        sha256_ctr_state.pos += take;
        p += take;
        len -= take;
    }
}

/* SHA-512 */

static struct {
    uint64_t counter;
    uint8_t buffer[64];
    size_t pos;
} sha512_ctr_state;

void cfx_sha512_ctr_seed(uint32_t seed) {
    sha512_ctr_state.counter = seed;
    sha512_ctr_state.pos = 64;
}

static void sha512_ctr_fill(void) {
    uint8_t input[8];
    cfx_store64_le(input, sha512_ctr_state.counter);
    sha512_ctr_state.counter++;
    cfx_sha512(sha512_ctr_state.buffer, input, 8);
    sha512_ctr_state.pos = 0;
}

uint32_t cfx_sha512_ctr_gen32(void) {
    if (sha512_ctr_state.pos >= 64) {
        sha512_ctr_fill();
    }
    uint32_t r = cfx_load32_le(&sha512_ctr_state.buffer[sha512_ctr_state.pos]);
    sha512_ctr_state.pos += 4;
    return r;
}

void cfx_sha512_ctr_bytes(void *buf, size_t len) {
    uint8_t *p = (uint8_t *)buf;
    while (len > 0) {
        if (sha512_ctr_state.pos >= 64) {
            sha512_ctr_fill();
        }
        size_t avail = 64 - sha512_ctr_state.pos;
        size_t take = (len < avail) ? len : avail;
        memcpy(p, &sha512_ctr_state.buffer[sha512_ctr_state.pos], take);
        sha512_ctr_state.pos += take;
        p += take;
        len -= take;
    }
}


/* BLAKE2b (64-byte output) */
static struct {
    uint64_t counter;
    uint8_t buffer[64];
    size_t pos;
} blake2b_ctr_state;

void cfx_blake2b_ctr_seed(uint32_t seed) {
    blake2b_ctr_state.counter = seed;
    blake2b_ctr_state.pos = 64;
}

static void blake2b_ctr_fill(void) {
    uint8_t input[8];
    cfx_store64_le(input, blake2b_ctr_state.counter);
    blake2b_ctr_state.counter++;
    cfx_blake2b(blake2b_ctr_state.buffer, 64, input, 8, NULL, 0);
    blake2b_ctr_state.pos = 0;
}

uint32_t cfx_blake2b_ctr_gen32(void) {
    if (blake2b_ctr_state.pos >= 64) {
        blake2b_ctr_fill();
    }
    uint32_t r = cfx_load32_le(&blake2b_ctr_state.buffer[blake2b_ctr_state.pos]);
    blake2b_ctr_state.pos += 4;
    return r;
}

void cfx_blake2b_ctr_bytes(void *buf, size_t len) {
    uint8_t *p = (uint8_t *)buf;
    while (len > 0) {
        if (blake2b_ctr_state.pos >= 64) {
            blake2b_ctr_fill();
        }
        size_t avail = 64 - blake2b_ctr_state.pos;
        size_t take = (len < avail) ? len : avail;
        memcpy(p, &blake2b_ctr_state.buffer[blake2b_ctr_state.pos], take);
        blake2b_ctr_state.pos += take;
        p += take;
        len -= take;
    }
}


/* BLAKE2s (32 byte output) */
static struct {
    uint64_t counter;
    uint8_t buffer[32];
    size_t pos;
} blake2s_ctr_state;

void cfx_blake2s_ctr_seed(uint32_t seed) {
    blake2s_ctr_state.counter = seed;
    blake2s_ctr_state.pos = 32;
}

static void blake2s_ctr_fill(void) {
    uint8_t input[8];
    cfx_store64_le(input, blake2s_ctr_state.counter);
    blake2s_ctr_state.counter++;
    cfx_blake2s(blake2s_ctr_state.buffer, 32, input, 8, NULL, 0);
    blake2s_ctr_state.pos = 0;
}

uint32_t cfx_blake2s_ctr_gen32(void) {
    if (blake2s_ctr_state.pos >= 32) {
        blake2s_ctr_fill();
    }
    uint32_t r = cfx_load32_le(&blake2s_ctr_state.buffer[blake2s_ctr_state.pos]);
    blake2s_ctr_state.pos += 4;
    return r;
}

void cfx_blake2s_ctr_bytes(void *buf, size_t len) {
    uint8_t *p = (uint8_t *)buf;
    while (len > 0) {
        if (blake2s_ctr_state.pos >= 32) {
            blake2s_ctr_fill();
        }
        size_t avail = 32 - blake2s_ctr_state.pos;
        size_t take = (len < avail) ? len : avail;
        memcpy(p, &blake2s_ctr_state.buffer[blake2s_ctr_state.pos], take);
        blake2s_ctr_state.pos += take;
        p += take;
        len -= take;
    }
}


const cfx_rand_desc_t g_hash_rngs[] = {
    {"cfx_sha256_ctr",  cfx_sha256_ctr_gen32,  cfx_sha256_ctr_seed,  cfx_sha256_ctr_bytes},
    {"cfx_sha512_ctr",  cfx_sha512_ctr_gen32,  cfx_sha512_ctr_seed,  cfx_sha512_ctr_bytes},
    {"cfx_blake2b_ctr", cfx_blake2b_ctr_gen32, cfx_blake2b_ctr_seed, cfx_blake2b_ctr_bytes},
    {"cfx_blake2s_ctr", cfx_blake2s_ctr_gen32, cfx_blake2s_ctr_seed, cfx_blake2s_ctr_bytes},
};

const size_t g_hash_rng_cnt = sizeof(g_hash_rngs) / sizeof(g_hash_rngs[0]);
