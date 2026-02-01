
#ifndef CFX_MEMORY_H
#define CFX_MEMORY_H

#include "cfx/arch.h"

#include <stddef.h>
#include <stdint.h>
#include <memory.h>


#ifdef __cplusplus
extern "C" {
#endif


/**
 * Q: Why do i keep the CFX_ALL_CAPS macros that just call the static inline functions?
 * A: To communicate that they're just like macros - they "inflate" at the call sites.
 **/


/* returns 0 if the data pointed to by a equals that pointed to by b*/
static inline int cfx_memeq_ct(const void *a, const void *b, size_t n) {
    const volatile uint8_t *pa = a;
    const volatile uint8_t *pb = b;
    uint8_t r = 0;
    size_t i;
    for (i = 0; i < n; i++) {
        r |= pa[i] ^ pb[i];
    }
    return r;
}

/* ---------------------------------------------------------------------------------------------- */
#define CFX_MEMZERO_S(ptr, size) cfx_memzero_s((ptr), (size))

/* make it impossible to optimize away a memory clear to avoid dead-store elimination,
   so we dont leave anything hanging around in RAM */
static inline void cfx_memzero_s(void *p, size_t n) {
    volatile unsigned char *v = (unsigned char *)p;
    while (n--) *v++ = 0;
}

/* ---------------------------------------------------------------------------------------------- */
#define CFX_STORE64_LE(dst, x) cfx_store64_le((dst), (x))

static inline void cfx_store64_le(void *dst, uint64_t x) {
#if CFX_LITTLE_ENDIAN
    memcpy(dst, &x, sizeof x);
#else
    uint8_t *p = (uint8_t *)dst;
    p[0] = (uint8_t)(x       );
    p[1] = (uint8_t)(x >>  8);
    p[2] = (uint8_t)(x >> 16);
    p[3] = (uint8_t)(x >> 24);
    p[4] = (uint8_t)(x >> 32);
    p[5] = (uint8_t)(x >> 40);
    p[6] = (uint8_t)(x >> 48);
    p[7] = (uint8_t)(x >> 56);
#endif
}

#define CFX_LOAD64_LE(src) cfx_load64_le((src))

static inline uint64_t cfx_load64_le(const void *src) {
#if CFX_LITTLE_ENDIAN
    uint64_t w;
    memcpy(&w, src, sizeof w);
    return w;
#else
    const uint8_t *p = (const uint8_t *)src;
    return ((uint64_t)p[0])        |
           ((uint64_t)p[1] <<  8)  |
           ((uint64_t)p[2] << 16)  |
           ((uint64_t)p[3] << 24)  |
           ((uint64_t)p[4] << 32)  |
           ((uint64_t)p[5] << 40)  |
           ((uint64_t)p[6] << 48)  |
           ((uint64_t)p[7] << 56);
#endif
}

#define CFX_LOAD32_LE_2(DST, SRC) cfx_load32_le_to((DST), (SRC))

static inline void cfx_load32_le_to(uint32_t *dst, const void *src) {
#if CFX_LITTLE_ENDIAN
    uint32_t v;
    memcpy(&v, src, sizeof v);
    *dst = v;
#else
    const uint8_t *p = (const uint8_t *)src;
    *dst = (uint32_t)p[0]       |
           (uint32_t)p[1] <<  8 |
           (uint32_t)p[2] << 16 |
           (uint32_t)p[3] << 24;
#endif
}

/* ---------------------------------------------------------------------------------------------- */
#define CFX_LOAD32_LE(SRC) cfx_load32_le(SRC)

static inline uint32_t cfx_load32_le(const void *src) {
#if CFX_LITTLE_ENDIAN
    uint32_t w;
    memcpy(&w, src, sizeof w);
    return w;
#else
    const uint8_t *p = (const uint8_t *)src;
    uint32_t w = (uint32_t)p[0];
    w |= (uint32_t)p[1] <<  8;
    w |= (uint32_t)p[2] << 16;
    w |= (uint32_t)p[3] << 24;
    return w;
#endif
}

/* ---------------------------------------------------------------------------------------------- */
#define CFX_STORE32_LE(DST, W) cfx_store32_le((DST), (W))

static inline void cfx_store32_le(void *dst, uint32_t w) {
#if CFX_LITTLE_ENDIAN
    memcpy(dst, &w, sizeof w);
#else
    uint8_t *p = (uint8_t *)dst;
    p[0] = (uint8_t)(w      );
    p[1] = (uint8_t)(w >>  8);
    p[2] = (uint8_t)(w >> 16);
    p[3] = (uint8_t)(w >> 24);
#endif
}


static inline void cfx_store16_le(void *dst, uint16_t v) {
#if CFX_LITTLE_ENDIAN
    memcpy(dst, &v, sizeof v);
#else
    dst[0] = (uint8_t)(v);
    dst[1] = (uint8_t)(v >> 8);
#endif
}

static inline uint16_t cfx_load16_le(const void *src) {
#if CFX_LITTLE_ENDIAN
    uint16_t w;
    memcpy(&w, src, sizeof w);
    return w;
#else
    const uint8_t *p = (const uint8_t *)src;
    uint16_t w = (uint16_t)p[0];
    w |= (uint16_t)p[1] << 8;
    return w;
#endif

}

#ifdef __cplusplus
}
#endif


#endif
