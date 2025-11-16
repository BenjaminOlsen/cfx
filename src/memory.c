#include "cfx/memory.h"


/* make it impossible to optimize away a memory clear to avoid dead-store elimination,
so we dont leave anything hanging around in RAM */
void cfx_memzero_s(void* p, size_t n) {
    volatile unsigned char* v = (unsigned char*)p;
    while (n--) *v++ = 0;
}

uint32_t cfx_load32_le(const void *p) {
    const unsigned char *b = (const unsigned char*)p;
    return ((uint32_t)b[0]) | ((uint32_t)b[1] << 8) |
           ((uint32_t)b[2] << 16) | ((uint32_t)b[3] << 24);
}

void cfx_store32_le(void *p, uint32_t x) {
    unsigned char* b = (unsigned char*)p;
    b[0] = (unsigned char)(x      );
    b[1] = (unsigned char)(x >>  8);
    b[2] = (unsigned char)(x >> 16);
    b[3] = (unsigned char)(x >> 24);
}

void cfx_store64_le(void *p, uint64_t x) {
    unsigned char *b = (unsigned char*)p;
    b[0] = (unsigned char)(x      );
    b[1] = (unsigned char)(x >>  8);
    b[2] = (unsigned char)(x >> 16);
    b[3] = (unsigned char)(x >> 24);
    b[4] = (unsigned char)(x >> 32);
    b[5] = (unsigned char)(x >> 40);
    b[6] = (unsigned char)(x >> 48);
    b[7] = (unsigned char)(x >> 56);
}

uint64_t cfx_load64_le(const void *p) {
    const uint8_t *b = (const uint8_t*)p;
    return ((uint64_t)b[0])       | ((uint64_t)b[1] << 8)  |
           ((uint64_t)b[2] << 16) | ((uint64_t)b[3] << 24) |
           ((uint64_t)b[4] << 32) | ((uint64_t)b[5] << 40) |
           ((uint64_t)b[6] << 48) | ((uint64_t)b[7] << 56);
}
