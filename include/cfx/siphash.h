#ifndef CFX_SIPHASH_H
#define CFX_SIPHASH_H

#include <stdint.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * SipHash - fast PRF optimized for short inputs
 *
 * Designed by Jean-Philippe Aumasson and Daniel J. Bernstein (2012)
 * Reference: https://cr.yp.to/siphash/siphash-20120918.pdf
 */

/* SipHash-2-4: standard 64-bit output */
uint64_t cfx_siphash(const uint8_t* data, size_t len, const uint8_t key[16]);

/* SipHash-c-d: configurable rounds (c=compression, d=finalization) */
uint64_t cfx_siphash_cd(const uint8_t* data, size_t len, const uint8_t key[16],
                        int c_rounds, int d_rounds);

/* SipHash-128: 128-bit output variant */
void cfx_siphash128(uint8_t out[16], const uint8_t* data, size_t len,
                    const uint8_t key[16]);

#define CFX_SIPHASH_CTX_SIZE 64u

typedef union {
    uint8_t  opaque[CFX_SIPHASH_CTX_SIZE];
    uint64_t aligner;
} cfx_siphash_ctx_t;

void cfx_siphash_init(cfx_siphash_ctx_t* ctx, const uint8_t key[16]);
void cfx_siphash_init_cd(cfx_siphash_ctx_t* ctx, const uint8_t key[16],
                         int c_rounds, int d_rounds);
void cfx_siphash_update(cfx_siphash_ctx_t* ctx, const uint8_t* data, size_t len);
uint64_t cfx_siphash_final(cfx_siphash_ctx_t* ctx);

#ifdef __cplusplus
}
#endif

#endif /* CFX_SIPHASH_H */
