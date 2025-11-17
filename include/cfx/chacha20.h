#ifndef CFX_CHACHA20_H
#define CFX_CHACHA20_H

#include "cfx/numerical.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------------------------------------------------ */
/* ChaCha20
 * refs
 *  https://datatracker.ietf.org/doc/html/rfc8439
 *  https://loup-vaillant.fr/tutorials/chacha20-design
 */
void cfx_chacha20_block_rfc8439(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                                uint8_t out[64]);

void cfx_chacha20_encrypt(const uint8_t key[32], uint32_t counter, const uint8_t nonce[12],
                          const uint8_t *pt, size_t pt_len, uint8_t *ct);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_CHACHA20_H */
