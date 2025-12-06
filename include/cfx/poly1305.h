#ifndef CFX_POLY1305_H
#define CFX_POLY1305_H

#include "cfx/types.h"
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

typedef struct {
    uint32_t r0, r1, r2, r3, r4;
    uint64_t s1, s2, s3, s4;
    uint32_t pad0, pad1, pad2, pad3;
    uint32_t h0, h1, h2, h3, h4;

    uint8_t  buffer[16];    /* holds any incomplete (< 16 byte ) blocks */
    uint8_t  buflen;        /* 0..15 */
    uint8_t  done;          /* After cfx_poly1305_finish(state, tag), the Poly1305 key must never be reused */
} cfx_poly1305_state_t;

void cfx_poly1305_init(cfx_poly1305_state_t* s, const uint8_t key[32]);
void cfx_poly1305_update(cfx_poly1305_state_t* s, const uint8_t* msg, size_t mlen);
void cfx_poly1305_finish(cfx_poly1305_state_t* s, uint8_t tag[16]);

void cfx_poly1305_mac(const uint8_t key[32], const uint8_t* msg, size_t mlen, uint8_t tag[16]);
void cfx_poly1305_mac_2(const uint8_t key[32], const uint8_t* msg, size_t mlen, uint8_t tag[16]);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_POLY1305_H */
