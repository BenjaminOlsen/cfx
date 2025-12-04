#ifndef CFX_ARCH_H
#define CFX_ARCH_H

#include <stdint.h>

#if defined(__clang__) || defined(__GNUC__)
#define CFX_SIMD 1
#endif

#if defined(__AVX2__)
#define CFX_HAVE_AVX2 1
#endif

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
    #include <stdalign.h>
    #define CFX_ALIGNOF(T) alignof(T)
    #define CFX_ALIGNAS(N) alignas(N)
#else
    #include <stddef.h>
    struct offset_struct_32 {char c; uint32_t x;};
    #define CFX_ALIGNOF(T) offsetof(struct offset_struct_32, x)
    #define CFX_ALIGNAS(N) __attribute__((aligned(N)))
#endif


/* -------------------------------------------------------- */
/* --- Endianness --- */
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define CFX_LITTLE_ENDIAN 1

/* Microsoft: always LE */
#elif defined(_MSC_VER)
#define CFX_LITTLE_ENDIAN 1

/* else use system headers if available */
#elif defined(__linux__) || defined(__ANDROID__)
#include <endian.h>
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define CFX_LITTLE_ENDIAN 1
#endif

#elif defined(__APPLE__)
#include <machine/endian.h>
#if _BYTE_ORDER == _LITTLE_ENDIAN
#define CFX_LITTLE_ENDIAN 1
#endif

#else
#define CFX_LITTLE_ENDIAN 0
#endif

#endif /* CFX_ARCH_H */
