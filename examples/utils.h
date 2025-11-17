#ifndef CFX_EXAMPLE_UTILS_H
#define CFX_EXAMPLE_UTILS_H

#include <stdint.h>
#include <stddef.h>
#include <string.h>

int hexval(int c);

/* parse hex string into exactly outlen bytes. returns 0 on success, -1 on error */
int parse_hex(const char* s, uint8_t* out, size_t outlen);
void print_hex(const uint8_t *buf, size_t len);

#endif  /* CFX_EXAMPLE_UTILS_H */
