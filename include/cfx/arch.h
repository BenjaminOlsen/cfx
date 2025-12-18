/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_ARCH_H
#define CFX_ARCH_H

#include <stdint.h>

#ifdef __cplusplus

    #if defined(__GNUC__) || defined(__clang__)
        #define CFX_RESTRICT __restrict__
    #elif defined(_MSC_VER)
        #define CFX_RESTRICT __restrict
    #else
        #define CFX_RESTRICT
    #endif

#else  /* ifndef __cplusplus */

    #if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
        #define CFX_RESTRICT restrict
    #else
        #define CFX_RESTRICT
    #endif

#endif

#if defined(__clang__) || defined(__GNUC__)
#define CFX_SIMD 1
#endif

#ifndef CFX_HAVE_AVX2
#  if defined(__AVX2__)
#    define CFX_HAVE_AVX2 1
#  else
#    define CFX_HAVE_AVX2 0
#  endif
#endif

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
    #include <stdalign.h>
    #define CFX_ALIGNOF(T) alignof(T)
    #define CFX_ALIGNAS(N) alignas(N)
#elif defined(_MSC_VER)
    #include <stddef.h>
    struct offset_struct_32 {char c; uint32_t x;};
    #define CFX_ALIGNOF(T) __alignof(T)
    #define CFX_ALIGNAS(N) __declspec(align(N))
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


/* --- x86 intrinsics include --- */
#if (defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)) \
    && (defined(__has_include) && __has_include(<immintrin.h>))
  #include <immintrin.h>
  #define CFX_USE_X86_INTRINSICS 1
#else
  #define CFX_USE_X86_INTRINSICS 0
#endif

/* -------- feature detection -------- */
#if !defined(CFX_FORCE_NO_UINT128) && defined(__SIZEOF_INT128__)
  #define CFX_HAS_UINT128 1
#endif

#if !defined(CFX_NO_FP)
  #include <float.h>
  #if defined(__STDC_IEC_559__) && DBL_MANT_DIG == 53
    #define CFX_USE_FP_ISQRT 1
  #else
    #define CFX_USE_FP_ISQRT 0
  #endif
#else
  #define CFX_USE_FP_ISQRT 0
#endif


/* -------- always_inline helper -------- */
#if defined(_MSC_VER)
  #define CFX_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
  #define CFX_INLINE __attribute__((always_inline)) static inline
#else
  #define CFX_INLINE static inline
#endif

#endif /* CFX_ARCH_H */
