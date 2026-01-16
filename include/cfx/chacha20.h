/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_CHACHA20_H
#define CFX_CHACHA20_H

#include <stdint.h>
#include <stddef.h>

#include "cfx/memory.h"
#include "cfx/arch.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------------------------------------------------ */
/* ChaCha20
 * refs
 *  https://datatracker.ietf.org/doc/html/rfc8439
 *  https://loup-vaillant.fr/tutorials/chacha20-design
 */

#define CFX_CHACHA20_ROTL32(x,n) ((uint32_t)(((x) << (n)) | ((x) >> (32-(n)))))

/* Quarter Round */
#define CFX_CHACHA20_QR(a, b, c, d) \
    a += b; d ^= a; d = CFX_CHACHA20_ROTL32(d,16); \
    c += d; b ^= c; b = CFX_CHACHA20_ROTL32(b,12); \
    a += b; d ^= a; d = CFX_CHACHA20_ROTL32(d, 8); \
    c += d; b ^= c; b = CFX_CHACHA20_ROTL32(b, 7);



typedef union {
    uint8_t opaque[64];
    uint64_t aligner;
} cfx_chacha20_ctx_t;


#if CFX_SIMD
typedef uint32_t vec4_u32 __attribute__((vector_size(16)));
#else
typedef struct {
    uint32_t v[4];
} vec4_u32;
#endif


typedef union {
    uint8_t opaque[256];
    uint64_t aligner;
} cfx_chacha20_ctx4_t;



#define CFX_CHACHA20_VQR(a,b,c,d)  \
    a = v_add(a, b); d = v_xor(d, a); d = v_rotl(d,16); \
    c = v_add(c, d); b = v_xor(b, c); b = v_rotl(b,12); \
    a = v_add(a, b); d = v_xor(d, a); d = v_rotl(d, 8); \
    c = v_add(c, d); b = v_xor(b, c); b = v_rotl(b, 7);



void cfx_chacha20_ctx_init(cfx_chacha20_ctx_t* ctx, const uint8_t key[32], const uint8_t nonce[12]);
void cfx_chacha20_block(cfx_chacha20_ctx_t* ctx, const uint32_t counter, uint8_t out[64]);
void cfx_chacha20_encrypt_ctx(cfx_chacha20_ctx_t* ctx, uint32_t* const counter, const uint8_t* pt, size_t pt_len, uint8_t* ct);

void cfx_chacha20_state_init4(cfx_chacha20_ctx4_t* ctx, const uint8_t key[32], const uint8_t nonce[4][12]);
void cfx_chacha20_block4(cfx_chacha20_ctx4_t* ctx, const uint32_t counter[4], uint8_t out[4][64]);

void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                                uint8_t out[64]);

/* Encrypt pt_len bytes. Max 256 GB per (key, nonce) due to 32-bit counter. */
void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t *pt, size_t pt_len, uint8_t *ct);

void cfx_chacha20_encrypt_bytes(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t *pt, size_t pt_len, uint8_t *ct);

/* HChaCha20 - derives 256-bit subkey for XChaCha20 */
void cfx_hchacha20(uint8_t out[32], const uint8_t key[32], const uint8_t nonce[16]);

/* XChaCha20 - ChaCha20 with extended 192 bit (24-byte) nonce */
void cfx_xchacha20_encrypt(const uint8_t key[32], uint32_t counter,
                           const uint8_t nonce[24],
                           const uint8_t* pt, size_t pt_len, uint8_t* ct);

void cfx_chacha20_block4_simd(const uint8_t key[32], const uint32_t counter[4],
                              const uint8_t nonce[4][12], uint8_t out[4][64]);

#if CFX_HAVE_AVX2

void cfx_chacha20_block8_avx2(const uint8_t key[32],
                              uint32_t counter,
                              const uint8_t nonce[12],
                              uint8_t out[8][64]);
#endif

#ifdef __cplusplus
}
#endif

#endif  /* CFX_CHACHA20_H */
