#ifndef CFX_CRYPTO_H
#define CFX_CRYPTO_H

#include "cfx/numerical.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter,const uint8_t nonce[12], uint8_t out[64]);
void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12], const uint8_t *pt, size_t pt_len, uint8_t *ct);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_CRYPTO_H */
