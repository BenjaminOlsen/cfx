
#ifndef CFX_MEMORY_H
#define CFX_MEMORY_H

#include <stddef.h>
#include <stdint.h>
#include <memory.h>

#ifdef __cplusplus
extern "C" {
#endif

void cfx_memzero_s(void* p, size_t n);
void cfx_store32_le(void *p, uint32_t x);
void cfx_store64_le(void *p, uint64_t x);
uint32_t cfx_load32_le(const void *p);
uint64_t cfx_load64_le(const void *p);

#ifdef __cplusplus
}
#endif


#endif
