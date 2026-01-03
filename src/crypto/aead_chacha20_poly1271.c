/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/aead_chacha20_poly1271.h"

#include <string.h>
#include "cfx/chacha20.h"
#include "cfx/poly1271.h"

static inline size_t pad15(size_t len) {
    size_t rem = len % 15;
    return rem ? (15 - rem) : 0;
}

static int memeq_ct(const void* a, const void* b, size_t n) {
    const volatile uint8_t* pa = a;
    const volatile uint8_t* pb = b;
    uint8_t diff = 0;
    for (size_t i = 0; i < n; i++) {
        diff |= pa[i] ^ pb[i];
    }
    return diff;
}

int cfx_chacha20_poly1271_encrypt(
    uint8_t*        ct,
    uint8_t         tag[16],
    const uint8_t*  pt,
    size_t          pt_len,
    const uint8_t*  aad,
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[12])
{
    uint8_t poly_key[64] = {0};
    uint8_t len_block[16] = {0};
    const uint8_t zeros[15] = {0};
    uint32_t counter = 0;
    cfx_chacha20_ctx_t chacha_ctx;

    cfx_chacha20_ctx_init(&chacha_ctx, key, nonce);
    cfx_chacha20_block(&chacha_ctx, counter, poly_key);

    counter = 1;
    cfx_chacha20_encrypt_ctx(&chacha_ctx, &counter, pt, pt_len, ct);

#if CFX_HAVE_AVX2
    cfx_poly1271_avx2_ctx_t poly_ctx;
    cfx_poly1271_avx2_init(&poly_ctx, poly_key);
    CFX_MEMZERO_S(poly_key, sizeof poly_key);

    cfx_poly1271_avx2_update(&poly_ctx, aad, aad_len);
    cfx_poly1271_avx2_update(&poly_ctx, zeros, pad15(aad_len));
    cfx_poly1271_avx2_update(&poly_ctx, ct, pt_len);
    cfx_poly1271_avx2_update(&poly_ctx, zeros, pad15(pt_len));

    CFX_STORE64_LE(len_block, (uint64_t)aad_len);
    CFX_STORE64_LE(len_block + 8, (uint64_t)pt_len);
    cfx_poly1271_avx2_update(&poly_ctx, len_block, 16);

    cfx_poly1271_avx2_finish(&poly_ctx, tag);
#else
    cfx_poly1271_ctx_t poly_ctx;
    cfx_poly1271_init(&poly_ctx, poly_key);
    CFX_MEMZERO_S(poly_key, sizeof poly_key);

    cfx_poly1271_update(&poly_ctx, aad, aad_len);
    cfx_poly1271_update(&poly_ctx, zeros, pad15(aad_len));
    cfx_poly1271_update(&poly_ctx, ct, pt_len);
    cfx_poly1271_update(&poly_ctx, zeros, pad15(pt_len));

    CFX_STORE64_LE(len_block, (uint64_t)aad_len);
    CFX_STORE64_LE(len_block + 8, (uint64_t)pt_len);
    cfx_poly1271_update(&poly_ctx, len_block, 16);

    cfx_poly1271_finish(&poly_ctx, tag);
#endif

    return 0;
}

int cfx_chacha20_poly1271_decrypt(
    uint8_t*        pt,
    const uint8_t*  ct,
    size_t          ct_len,
    const uint8_t*  aad,
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[12],
    const uint8_t   tag[16])
{
    uint8_t poly_key[64] = {0};
    uint8_t len_block[16] = {0};
    uint8_t computed_tag[16] = {0};
    const uint8_t zeros[15] = {0};
    uint32_t counter = 0;
    cfx_chacha20_ctx_t chacha_ctx;

    cfx_chacha20_ctx_init(&chacha_ctx, key, nonce);
    cfx_chacha20_block(&chacha_ctx, counter, poly_key);

#if CFX_HAVE_AVX2
    cfx_poly1271_avx2_ctx_t poly_ctx;
    cfx_poly1271_avx2_init(&poly_ctx, poly_key);
    CFX_MEMZERO_S(poly_key, sizeof poly_key);

    cfx_poly1271_avx2_update(&poly_ctx, aad, aad_len);
    cfx_poly1271_avx2_update(&poly_ctx, zeros, pad15(aad_len));
    cfx_poly1271_avx2_update(&poly_ctx, ct, ct_len);
    cfx_poly1271_avx2_update(&poly_ctx, zeros, pad15(ct_len));

    CFX_STORE64_LE(len_block, (uint64_t)aad_len);
    CFX_STORE64_LE(len_block + 8, (uint64_t)ct_len);
    cfx_poly1271_avx2_update(&poly_ctx, len_block, 16);

    cfx_poly1271_avx2_finish(&poly_ctx, computed_tag);
#else
    cfx_poly1271_ctx_t poly_ctx;
    cfx_poly1271_init(&poly_ctx, poly_key);
    CFX_MEMZERO_S(poly_key, sizeof poly_key);

    cfx_poly1271_update(&poly_ctx, aad, aad_len);
    cfx_poly1271_update(&poly_ctx, zeros, pad15(aad_len));
    cfx_poly1271_update(&poly_ctx, ct, ct_len);
    cfx_poly1271_update(&poly_ctx, zeros, pad15(ct_len));

    CFX_STORE64_LE(len_block, (uint64_t)aad_len);
    CFX_STORE64_LE(len_block + 8, (uint64_t)ct_len);
    cfx_poly1271_update(&poly_ctx, len_block, 16);

    cfx_poly1271_finish(&poly_ctx, computed_tag);
#endif

    if (memeq_ct(tag, computed_tag, 16) != 0) {
        CFX_MEMZERO_S(computed_tag, sizeof computed_tag);
        CFX_MEMZERO_S(&chacha_ctx, sizeof chacha_ctx);
        return -1;
    }

    counter = 1;
    cfx_chacha20_encrypt_ctx(&chacha_ctx, &counter, ct, ct_len, pt);

    CFX_MEMZERO_S(&chacha_ctx, sizeof chacha_ctx);
    CFX_MEMZERO_S(computed_tag, sizeof computed_tag);
    return 0;
}

int cfx_xchacha20_poly1271_encrypt(
    uint8_t*        ct,
    uint8_t         tag[16],
    const uint8_t*  pt,
    size_t          pt_len,
    const uint8_t*  aad,
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[24])
{
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};

    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);

    int rc = cfx_chacha20_poly1271_encrypt(ct, tag, pt, pt_len, aad, aad_len, subkey, subnonce);

    CFX_MEMZERO_S(subkey, sizeof subkey);
    return rc;
}

int cfx_xchacha20_poly1271_decrypt(
    uint8_t*        pt,
    const uint8_t*  ct,
    size_t          ct_len,
    const uint8_t*  aad,
    size_t          aad_len,
    const uint8_t   key[32],
    const uint8_t   nonce[24],
    const uint8_t   tag[16])
{
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};

    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);

    int rc = cfx_chacha20_poly1271_decrypt(pt, ct, ct_len, aad, aad_len, subkey, subnonce, tag);

    CFX_MEMZERO_S(subkey, sizeof subkey);
    return rc;
}
