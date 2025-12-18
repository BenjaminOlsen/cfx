#include "utils.h"
#include <stdio.h>

int hexval(int c) {
    if ('0' <= c && c <= '9') return c - '0';
    if ('a' <= c && c <= 'f') return 10 + c - 'a';
    if ('A' <= c && c <= 'F') return 10 + c - 'A';
    return -1;
}

/* parse hex string into exactly outlen bytes. returns 0 on success, -1 on error */
int cfx_parse_hex(const char* s, uint8_t* out, size_t outlen) {
    if (s[0]=='0' && (s[1]=='x' || s[1]=='X')) s += 2;

    size_t n = strlen(s);
    if (n != outlen*2) return -1;

    for (size_t i = 0; i < outlen; ++i) {
        int hi = hexval(s[2*i]);
        int lo = hexval(s[2*i+1]);
        if (hi < 0 || lo < 0) return -1;
        out[i] = (uint8_t)((hi << 4) | lo);
    }
    return 0;
}

void cfx_print_hex(const uint8_t *buf, size_t len) {
    for (size_t i = 0; i < len; ++i) {
        printf("%02x", (unsigned)buf[i]);
    }
}
