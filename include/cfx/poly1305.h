#ifndef CFX_POLY1305_H
#define CFX_POLY1305_H

#include "cfx/arith.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ------------------------------------------------------------------------------------------------------------------ */
/* Poly1305:
 * refs
 *  https://datatracker.ietf.org/doc/html/rfc8439
 *  https://loup-vaillant.fr/tutorials/poly1305-design
 */

#define CFX_POLY1305_CTX_MAX_SIZE 224u

typedef union {
    uint8_t opaque[CFX_POLY1305_CTX_MAX_SIZE];
    size_t align;
} cfx_poly1305_ctx_t;


void cfx_poly1305_init(cfx_poly1305_ctx_t* ctx, const uint8_t key[32]);
void cfx_poly1305_update(cfx_poly1305_ctx_t* ctx, const uint8_t* msg, size_t mlen);
void cfx_poly1305_finish(cfx_poly1305_ctx_t* ctx, uint8_t tag[16]);

void cfx_poly1305_mac(const uint8_t key[32], const uint8_t* msg, size_t mlen, uint8_t tag[16]);
void cfx_poly1305_mac_2(const uint8_t key[32], const uint8_t* msg, size_t mlen, uint8_t tag[16]);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_POLY1305_H */
