#ifndef CFX_POLY1305_H
#define CFX_POLY1305_H

#include "cfx/numerical.h"
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

void cfx_poly1305_mac(const uint8_t key[32], const uint8_t* msg, size_t mlen, uint8_t tag[16]);

#ifdef __cplusplus
}
#endif

#endif  /* CFX_POLY1305_H */
