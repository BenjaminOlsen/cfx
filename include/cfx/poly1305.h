#ifndef CFX_POLY1305_H
#define CFX_POLY1305_H

#include "cfx/arch.h"
#include "cfx/arith.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Poly1305:
 * refs
 *  https://datatracker.ietf.org/doc/html/rfc8439
 *  https://loup-vaillant.fr/tutorials/poly1305-design
 */

#define CFX_POLY1305_CTX_MAX_SIZE 288u

typedef union {
    CFX_ALIGNAS(CFX_CTX_ALIGN) uint8_t opaque[CFX_POLY1305_CTX_MAX_SIZE];
} cfx_poly1305_ctx_t;


void cfx_poly1305_init(cfx_poly1305_ctx_t *ctx, const uint8_t key[32]);
void cfx_poly1305_update(cfx_poly1305_ctx_t *ctx, const uint8_t *msg, size_t mlen);
void cfx_poly1305_finish(cfx_poly1305_ctx_t *ctx, uint8_t tag[16]);

/* one shot */
void cfx_poly1305(uint8_t tag[16], const uint8_t *msg, size_t mlen, const uint8_t key[32]);


#ifdef __cplusplus
}
#endif

#endif  /* CFX_POLY1305_H */
