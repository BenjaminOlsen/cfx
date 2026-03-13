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

/* -------- always_inline helper -------- */
#if defined(_MSC_VER)
  #define CFX_INLINE static __forceinline
#elif defined(__GNUC__) || defined(__clang__)
  #define CFX_INLINE static inline __attribute__((always_inline))
#else
  #define CFX_INLINE static inline
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

/* -------------------------------------------------------- */
/* Alignment -  every compiler does this differently in C99, so
  if you ever need to add some unlisted compiler, add it... */

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
    #include <stdalign.h>
    #define CFX_ALIGNOF(T) alignof(T)
    #define CFX_ALIGNAS(N) alignas(N)
#elif defined(_MSC_VER)
    #define CFX_ALIGNOF(T) __alignof(T)
    #define CFX_ALIGNAS(N) __declspec(align(N))
#elif defined(__GNUC__) || defined(__clang__)
    #define CFX_ALIGNOF(T) __alignof__(T)
    #define CFX_ALIGNAS(N) __attribute__((aligned(N)))
#else
    /* best-effort fallback
       ALIGNAS is a no-op — natural alignment must suffice */
    #include <stddef.h>
    #define CFX_ALIGNOF(T) offsetof(struct { char c; T x; }, x)
    #define CFX_ALIGNAS(N)
#endif

/* baseline alignment for opaque context unions  */
#if defined(__LP64__) || defined(_WIN64) || defined(__x86_64__) || defined(__aarch64__)
#  define CFX_CTX_ALIGN 16
#else
#  define CFX_CTX_ALIGN 8
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
#else
#define CFX_LITTLE_ENDIAN 0
#endif

#elif defined(__APPLE__)
#include <machine/endian.h>
#if _BYTE_ORDER == _LITTLE_ENDIAN
#define CFX_LITTLE_ENDIAN 1
#else
#define CFX_LITTLE_ENDIAN 0
#endif

#else
#define CFX_LITTLE_ENDIAN 0
#endif


/* --- x86 intrinsics include --- */
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#  if defined(__has_include)
#    if __has_include(<immintrin.h>)
#      include <immintrin.h>
#      define CFX_USE_X86_INTRINSICS 1
#    else
#      define CFX_USE_X86_INTRINSICS 0
#    endif
#  elif defined(_MSC_VER)
     /* MSVC always has immintrin.h for x86/x64 targets */
#    include <immintrin.h>
#    define CFX_USE_X86_INTRINSICS 1
#  elif defined(__GNUC__) || defined(__clang__)
     /* GCC/Clang without __has_include likely has immintrin.h */
#    include <immintrin.h>
#    define CFX_USE_X86_INTRINSICS 1
#  else
#    define CFX_USE_X86_INTRINSICS 0
#  endif
#else
#  define CFX_USE_X86_INTRINSICS 0
#endif

/* -------- 128-bit arithmetic detection -------- */

/*
 * CFX_HAS_128BIT: Set if we have fast 128-bit multiply/mod capability
 * CFX_128BIT_NATIVE: GCC/Clang __uint128_t
 * CFX_128BIT_MSVC: MSVC _umul128/_udiv128 intrinsics
 */
#if !defined(CFX_FORCE_NO_UINT128) && defined(__SIZEOF_INT128__)
  #define CFX_HAS_UINT128 1
  #define CFX_HAS_128BIT 1
  #define CFX_128BIT_NATIVE 1
#elif defined(_MSC_VER) && defined(_M_X64)
  #include <intrin.h>
  #define CFX_HAS_128BIT 1
  #define CFX_128BIT_MSVC 1
  #if _MSC_VER >= 1920
    #define CFX_MSVC_HAS_UDIV128 1
  #endif
#else
  #define CFX_HAS_128BIT 0
#endif

/* -------- 64-bit modular multiplication primitive -------- */

/*
 * cfx_mulmod64(a, b, m) - computes (a * b) % m without overflow
 * Uses the best available method for the platform.
 */
#if defined(CFX_128BIT_NATIVE)

CFX_INLINE uint64_t cfx_mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    __uint128_t p = (__uint128_t)a * (__uint128_t)b;
    return (uint64_t)(p % m);
}

#elif defined(CFX_128BIT_MSVC) && defined(CFX_MSVC_HAS_UDIV128)

CFX_INLINE uint64_t cfx_mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t hi, lo;
    lo = _umul128(a, b, &hi);
    uint64_t rem;
    (void)_udiv128(hi, lo, m, &rem);
    return rem;
}

#elif defined(CFX_128BIT_MSVC)

/* MSVC x64 without _udiv128 (pre-VS2019) */
CFX_INLINE uint64_t cfx_mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t hi, lo;
    lo = _umul128(a, b, &hi);
    if (hi == 0) return lo % m;
    /* Shift-and-subtract for 128/64 mod */
    uint64_t r = hi % m;
    for (int i = 63; i >= 0; --i) {
        r = (r < m - r) ? (r + r) : (r - (m - r));
        if ((lo >> i) & 1) {
            r = (r + 1 < m) ? (r + 1) : 0;
        }
    }
    return r;
}

#else

/* Portable fallback: Russian peasant (double-and-add) */
CFX_INLINE uint64_t cfx_mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    if (m <= 1) return 0;
    a %= m;
    b %= m;
    if (a > b) {
        uint64_t t = a; a = b; b = t;
    }
    uint64_t result = 0;
    while (a > 0) {
        if (a & 1) {
            result = (result < m - b) ? (result + b) : (result - (m - b));
        }
        a >>= 1;
        b = (b < m - b) ? (b + b) : (b - (m - b));
    }
    return result;
}

#endif

/* -------- 64-bit modular exponentiation primitive -------- */

CFX_INLINE uint64_t cfx_powmod64(uint64_t base, uint64_t exp, uint64_t m) {
    if (m <= 1) return 0;
    base %= m;
    uint64_t result = 1;
    while (exp > 0) {
        if (exp & 1) result = cfx_mulmod64(result, base, m);
        base = cfx_mulmod64(base, base, m);
        exp >>= 1;
    }
    return result;
}

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


#endif  /* CFX_ARCH_H */
