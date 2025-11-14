#include "cfx/memory.h"

/* make it impossible to optimize away a memory clear to avoid dead-store elimination,
so we dont leave anything hanging around in RAM */
void cfx_memzero_s(void* p, size_t n) {
    volatile unsigned char* v = (unsigned char*)p;
    while (n--) *v++ = 0;
}
