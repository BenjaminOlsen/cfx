/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#include "cfx/aead_chacha20_poly1305.h"

#include <string.h>
#include "cfx/chacha20.h"
#include "cfx/poly1305.h"
#include "cfx/memory.h"
#include "cfx/arch.h"

static inline size_t ct_pad16(size_t len) {
    size_t rem = len & 0xF;            /* 0..15 */
    return (16 - rem) & 0xF;           /* yields 0 if rem==0, else 16-rem */
}

int cfx_chacha20_poly1305_encrypt(
    uint8_t *ct,                    /* cyphertext - out, length = pt_len */
    uint8_t tag[16],                /* out */
    const uint8_t *pt,              /* plaintext  - in,  length = pt_len */
    size_t pt_len,
    const uint8_t *aad,             /* in, may be NULL (iff aad_len == 0) */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[12]) {

    /* poly_key = ChaCha20(key, nonce, counter=0)[0..31] */
    uint8_t poly_key[64] = {0};  /* chacha20 writes 64 bytes, poly1305 uses least significant 32 */
    uint8_t len_block[8] = {0};
    const uint8_t zeros[16] = {0};
    uint32_t counter = 0;
    size_t pad = 0;
    cfx_chacha20_ctx_t chacha_ctx;
    cfx_poly1305_ctx_t poly_ctx;

    cfx_chacha20_ctx_init(&chacha_ctx, key, nonce);
    cfx_chacha20_block(&chacha_ctx, counter, poly_key);
    cfx_poly1305_init(&poly_ctx, poly_key);
    CFX_MEMZERO_S(&poly_key, sizeof poly_key);

    /* ciphertext = ChaCha20(key, nonce, counter=1)[pt_len] XOR plaintext */
    counter = 1;
    cfx_chacha20_encrypt_ctx(&chacha_ctx, &counter, pt, pt_len, ct);

    /* make the poly1305 tag out of the message: */
    /* AAD || padding1 || ciphertext || padding2 || len(AAD) (uint64, little-endian) || len(ciphertext) (uint64, little-endian) */
    cfx_poly1305_update(&poly_ctx, aad, aad_len); /* does nothing on aad_len == 0*/
    pad = ct_pad16(aad_len);
    cfx_poly1305_update(&poly_ctx, zeros, pad);
    cfx_poly1305_update(&poly_ctx, ct, pt_len);
    pad = ct_pad16(pt_len);
    cfx_poly1305_update(&poly_ctx, zeros, pad);
    CFX_STORE64_LE(len_block, (uint64_t)aad_len);
    cfx_poly1305_update(&poly_ctx, len_block, 8);
    CFX_STORE64_LE(len_block, (uint64_t)pt_len);
    cfx_poly1305_update(&poly_ctx, len_block, 8);

    cfx_poly1305_finish(&poly_ctx, tag);
    return 0;
}



int cfx_chacha20_poly1305_decrypt(
    uint8_t *pt,                    /* plaintext  - out, length = ct_len */
    const uint8_t *ct,              /* ciphertext - in,  length = ct_len */
    size_t ct_len,
    const uint8_t *aad,             /* in, may be NULL (iff aad_len == 0) */
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[12],
    const uint8_t tag[16]){
    uint8_t poly_key[64] = {0};
    uint8_t len_block[8] = {0};
    const uint8_t zeros[16] = {0};
    uint8_t computed_tag[16] = {0};
    uint32_t counter = 0;
    size_t pad = 0;

    cfx_chacha20_ctx_t chacha_ctx;
    cfx_poly1305_ctx_t poly_ctx;

    cfx_chacha20_ctx_init(&chacha_ctx, key, nonce);
    cfx_chacha20_block(&chacha_ctx, counter, poly_key);
    cfx_poly1305_init(&poly_ctx, poly_key);
    CFX_MEMZERO_S(&poly_key, sizeof poly_key);

    cfx_poly1305_update(&poly_ctx, aad, aad_len);
    pad = ct_pad16(aad_len);
    cfx_poly1305_update(&poly_ctx, zeros, pad);
    cfx_poly1305_update(&poly_ctx, ct, ct_len);
    pad = ct_pad16(ct_len);
    cfx_poly1305_update(&poly_ctx, zeros, pad);
    CFX_STORE64_LE(len_block, (uint64_t)aad_len);
    cfx_poly1305_update(&poly_ctx, len_block, 8);
    CFX_STORE64_LE(len_block, (uint64_t)ct_len);
    cfx_poly1305_update(&poly_ctx, len_block, 8);
    cfx_poly1305_finish(&poly_ctx, computed_tag);

    if (cfx_memeq_ct(tag, computed_tag, sizeof computed_tag) != 0) {
        CFX_MEMZERO_S(&computed_tag, sizeof computed_tag);
        CFX_MEMZERO_S(&chacha_ctx, sizeof chacha_ctx);
        CFX_MEMZERO_S(&poly_ctx, sizeof poly_ctx);
        return -1;  /* authentication error */
    }

    counter = 1;
    cfx_chacha20_encrypt_ctx(&chacha_ctx, &counter, ct, ct_len, pt);
    CFX_MEMZERO_S(&chacha_ctx, sizeof chacha_ctx);
    return 0;
}

int cfx_xchacha20_poly1305_encrypt(
    uint8_t *ct,
    uint8_t tag[16],
    const uint8_t *pt,
    size_t pt_len,
    const uint8_t *aad,
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[24]) {
        
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};

    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);

    int rc = cfx_chacha20_poly1305_encrypt(ct, tag, pt, pt_len, aad, aad_len, subkey, subnonce);

    CFX_MEMZERO_S(subkey, sizeof subkey);
    return rc;
}

int cfx_xchacha20_poly1305_decrypt(
    uint8_t *pt,
    const uint8_t *ct,
    size_t ct_len,
    const uint8_t *aad,
    size_t aad_len,
    const uint8_t key[32],
    const uint8_t nonce[24],
    const uint8_t tag[16]){
    uint8_t subkey[32];
    uint8_t subnonce[12] = {0};

    cfx_hchacha20(subkey, key, nonce);
    memcpy(subnonce + 4, nonce + 16, 8);

    int rc = cfx_chacha20_poly1305_decrypt(pt, ct, ct_len, aad, aad_len, subkey, subnonce, tag);

    CFX_MEMZERO_S(subkey, sizeof subkey);
    return rc;
}

static void stream_derive_nonce(uint8_t out[24], const uint8_t base[24],
                                uint64_t counter, int is_final) {
    memcpy(out, base, 24);
    out[20] ^= (uint8_t)(counter);
    out[21] ^= (uint8_t)(counter >> 8);
    out[22] ^= (uint8_t)(counter >> 16);
    /* use the top bit for 'is_final', the rest for the counter */
    out[23] ^= (uint8_t)((counter >> 24) & 0x7F) | (is_final ? 0x80 : 0x00);
}

int cfx_stream_xchacha20_poly1305_encrypt_chunk(
    uint8_t *ct, uint8_t tag[16],
    const uint8_t *pt, size_t pt_len,
    uint64_t chunk_counter, int is_final,
    const uint8_t key[32], const uint8_t base_nonce[24]) {

    uint8_t nonce[24];
    stream_derive_nonce(nonce, base_nonce, chunk_counter, is_final);
    int rc = cfx_xchacha20_poly1305_encrypt(ct, tag, pt, pt_len, NULL, 0, key, nonce);
    CFX_MEMZERO_S(nonce, sizeof nonce);
    return rc;
}

int cfx_stream_xchacha20_poly1305_decrypt_chunk(
    uint8_t *pt,
    const uint8_t *ct, size_t ct_len,
    const uint8_t tag[16],
    uint64_t chunk_counter, int is_final,
    const uint8_t key[32], const uint8_t base_nonce[24]) {

    uint8_t nonce[24];
    stream_derive_nonce(nonce, base_nonce, chunk_counter, is_final);
    int rc = cfx_xchacha20_poly1305_decrypt(pt, ct, ct_len, NULL, 0, key, nonce, tag);
    CFX_MEMZERO_S(nonce, sizeof nonce);
    return rc;
}
