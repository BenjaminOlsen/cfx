/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/hmac.h"
#include "cfx/sha256.h"
#include "cfx/sha512.h"
#include "cfx/memory.h"
#include "cfx/macros.h"

#include <string.h>


#define SHA256_BLOCK  64
#define SHA256_DIGEST 32

typedef struct {
    cfx_sha256_ctx inner;
    uint8_t        opad[SHA256_BLOCK];
} hmac_sha256_state;

CFX_STATIC_ASSERT(
    sizeof(hmac_sha256_state) <= CFX_HMAC_SHA256_CTX_SIZE,
    CFX_HMAC_SHA256_CTX_SIZE_too_small);

void cfx_hmac_sha256_init(cfx_hmac_sha256_ctx *ctx, const uint8_t *key, size_t key_len) {

    hmac_sha256_state *st = (hmac_sha256_state *)ctx;
    uint8_t kp[SHA256_BLOCK];
    size_t i;

    /* K' = H(K) if key > block size, else K zero-padded */
    memset(kp, 0, SHA256_BLOCK);
    if (key_len > SHA256_BLOCK) {
        cfx_sha256(kp, key, key_len);
    } else {
        memcpy(kp, key, key_len);
    }

    /* Build opad for later, start inner hash with ipad */
    uint8_t ipad[SHA256_BLOCK];
    for (i = 0; i < SHA256_BLOCK; i++) {
        ipad[i]     = kp[i] ^ 0x36;
        st->opad[i] = kp[i] ^ 0x5c;
    }

    CFX_MEMZERO_S(kp, sizeof kp);

    cfx_sha256_init(&st->inner);
    cfx_sha256_update(&st->inner, ipad, SHA256_BLOCK);

    CFX_MEMZERO_S(ipad, sizeof ipad);
}

void cfx_hmac_sha256_update(cfx_hmac_sha256_ctx *ctx, const uint8_t *data, size_t len) {

    hmac_sha256_state *st = (hmac_sha256_state *)ctx;
    cfx_sha256_update(&st->inner, data, len);
}

void cfx_hmac_sha256_final(cfx_hmac_sha256_ctx *ctx, uint8_t out[32]) {

    hmac_sha256_state *st = (hmac_sha256_state *)ctx;
    uint8_t inner_digest[SHA256_DIGEST];

    cfx_sha256_final(&st->inner, inner_digest);

    /* outer = H(opad || inner_digest) */
    cfx_sha256_ctx outer;
    cfx_sha256_init(&outer);
    cfx_sha256_update(&outer, st->opad, SHA256_BLOCK);
    cfx_sha256_update(&outer, inner_digest, SHA256_DIGEST);
    cfx_sha256_final(&outer, out);

    CFX_MEMZERO_S(inner_digest, sizeof inner_digest);
    CFX_MEMZERO_S(st, sizeof *st);
}

void cfx_hmac_sha256(uint8_t out[32], const uint8_t *key, size_t key_len, const uint8_t *data, size_t data_len) {
    cfx_hmac_sha256_ctx ctx;
    cfx_hmac_sha256_init(&ctx, key, key_len);
    cfx_hmac_sha256_update(&ctx, data, data_len);
    cfx_hmac_sha256_final(&ctx, out);
}


#define SHA512_BLOCK  128
#define SHA512_DIGEST 64

typedef struct {
    cfx_sha512_ctx_t inner;
    uint8_t          opad[SHA512_BLOCK];
} hmac_sha512_state;

CFX_STATIC_ASSERT(
    sizeof(hmac_sha512_state) <= CFX_HMAC_SHA512_CTX_SIZE,
    CFX_HMAC_SHA512_CTX_SIZE_too_small);

void cfx_hmac_sha512_init(cfx_hmac_sha512_ctx *ctx, const uint8_t *key, size_t key_len) {

    hmac_sha512_state *st = (hmac_sha512_state *)ctx;
    uint8_t kp[SHA512_BLOCK];
    size_t i;

    memset(kp, 0, SHA512_BLOCK);
    if (key_len > SHA512_BLOCK) {
        cfx_sha512(kp, key, key_len);
    } else {
        memcpy(kp, key, key_len);
    }

    uint8_t ipad[SHA512_BLOCK];
    for (i = 0; i < SHA512_BLOCK; i++) {
        ipad[i]     = kp[i] ^ 0x36;
        st->opad[i] = kp[i] ^ 0x5c;
    }

    CFX_MEMZERO_S(kp, sizeof kp);

    cfx_sha512_init(&st->inner);
    cfx_sha512_update(&st->inner, ipad, SHA512_BLOCK);

    CFX_MEMZERO_S(ipad, sizeof ipad);
}

void cfx_hmac_sha512_update(cfx_hmac_sha512_ctx *ctx, const uint8_t *data, size_t len) {
    hmac_sha512_state *st = (hmac_sha512_state *)ctx;
    cfx_sha512_update(&st->inner, data, len);
}

void cfx_hmac_sha512_final(cfx_hmac_sha512_ctx *ctx, uint8_t out[64]) {
    hmac_sha512_state *st = (hmac_sha512_state *)ctx;
    uint8_t inner_digest[SHA512_DIGEST];

    cfx_sha512_final(&st->inner, inner_digest);

    cfx_sha512_ctx_t outer;
    cfx_sha512_init(&outer);
    cfx_sha512_update(&outer, st->opad, SHA512_BLOCK);
    cfx_sha512_update(&outer, inner_digest, SHA512_DIGEST);
    cfx_sha512_final(&outer, out);

    CFX_MEMZERO_S(inner_digest, sizeof inner_digest);
    CFX_MEMZERO_S(st, sizeof *st);
}

void cfx_hmac_sha512(uint8_t out[64],
                      const uint8_t *key, size_t key_len,
                      const uint8_t *data, size_t data_len) {

    cfx_hmac_sha512_ctx ctx;
    cfx_hmac_sha512_init(&ctx, key, key_len);
    cfx_hmac_sha512_update(&ctx, data, data_len);
    cfx_hmac_sha512_final(&ctx, out);
}
