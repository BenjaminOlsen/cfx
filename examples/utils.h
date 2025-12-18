#ifndef CFX_EXAMPLE_UTILS_H
#define CFX_EXAMPLE_UTILS_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

int hexval(int c);

/* parse hex string into exactly outlen bytes. returns 0 on success, -1 on error */
int cfx_parse_hex(const char* s, uint8_t* out, size_t outlen);
void cfx_print_hex(const uint8_t *buf, size_t len);

int cfx_base64_decode(const char *s, uint8_t *out, size_t outlen);

#endif  /* CFX_EXAMPLE_UTILS_H */
