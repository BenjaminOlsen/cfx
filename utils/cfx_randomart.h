/* cfx_randomart.h - drunken bishop randomart visualization */

#ifndef CFX_RANDOMART_H
#define CFX_RANDOMART_H

#include <stddef.h>
#include <stdint.h>

#define CFX_RA_DEFAULT_WIDTH  17
#define CFX_RA_DEFAULT_HEIGHT  9
#define CFX_RA_MAX_WIDTH      79
#define CFX_RA_MAX_HEIGHT     39

/* walk the bishop on a w x h grid over raw bytes.
   out must be char[h][w+1] (each row NUL-terminated).
   w and h must be odd. */
void cfx_randomart(char *out, int w, int h,
                   const uint8_t *data, size_t len);

/* convenience: SHA-256 hash then walk on default 17x9 grid, print framed.
   matches OpenSSH behaviour. */
void cfx_print_randomart(const uint8_t *data, size_t len,
                         const char *label);

/* print framed randomart from pre-computed art buffer */
void cfx_print_randomart_frame(const char *art, int w, int h,
                               const char *label);

#endif
