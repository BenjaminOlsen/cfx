#ifndef CFX_MISC_H
#define CFX_MISC_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

int hexval(int c);

/* parse hex string into exactly outlen bytes. returns 0 on success, -1 on error */
int cfx_parse_hex(const char* s, uint8_t* out, size_t outlen);
void cfx_print_hex(const uint8_t *buf, size_t len);

#endif  /* CFX_MISC_H */
