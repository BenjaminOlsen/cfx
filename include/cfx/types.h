#ifndef CFX_PRIM_H
#define CFX_PRIM_H

#include <stdint.h>

/* Capability detection (override in your build with -DCFX_FORCE_NO_UINT128 if needed) */
#if !defined(CFX_FORCE_NO_UINT128) && defined(__SIZEOF_INT128__)
  #define CFX_HAS_UINT128 1
#endif

#if defined(UINT64_MAX)
  #define CFX_HAS_UINT64 1
#endif

/* Public “canonical” names. Keep your code using these typedefs wherever possible. */
#if CFX_HAS_UINT64
  typedef uint64_t     cfx_u64_t;
#else
  /* if you ever need to emulate u64, you can change this typedef later */
  typedef uint64_t     cfx_u64_t; /* keep for now */
#endif

#if CFX_HAS_UINT128
  typedef __uint128_t  cfx_u128_t;
#else
  /* later: provide an emulated struct here */
  typedef struct { uint64_t lo, hi; } cfx_u128_t; /* placeholder for future */
#endif

/* Force-inline helper */
#if defined(_MSC_VER)
  #define CFX_ALWAYS_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
  #define CFX_ALWAYS_INLINE __attribute__((always_inline)) inline
#else
  #define CFX_ALWAYS_INLINE inline
#endif

/* =========================
 *  128-bit helpers (native)
 * ========================= */
#if CFX_HAS_UINT128

/* Constructors / accessors */
CFX_ALWAYS_INLINE cfx_u128_t cfx_u128_zero(void) { return (cfx_u128_t)0; }
CFX_ALWAYS_INLINE cfx_u128_t cfx_u128_from_u64(cfx_u64_t x) { return (cfx_u128_t)x; }
CFX_ALWAYS_INLINE cfx_u64_t  cfx_u128_lo64(cfx_u128_t x) { return (cfx_u64_t)x; }
CFX_ALWAYS_INLINE cfx_u64_t  cfx_u128_hi64(cfx_u128_t x) { return (cfx_u64_t)(x >> 64); }

/* In-place adds/subs (replace: acc += x;  acc -= x;) */
CFX_ALWAYS_INLINE void cfx_u128_add_eq(cfx_u128_t* acc, cfx_u64_t add) { *acc += (cfx_u128_t)add; }
CFX_ALWAYS_INLINE void cfx_u128_add128_eq(cfx_u128_t* acc, cfx_u128_t add) { *acc += add; }
CFX_ALWAYS_INLINE void cfx_u128_sub_eq(cfx_u128_t* acc, cfx_u64_t sub) { *acc -= (cfx_u128_t)sub; }
CFX_ALWAYS_INLINE void cfx_u128_sub128_eq(cfx_u128_t* acc, cfx_u128_t sub) { *acc -= sub; }

/* Shifts (replace: x <<= k; x >>= k;) */
CFX_ALWAYS_INLINE void cfx_u128_shl_eq(cfx_u128_t* x, unsigned k) { *x <<= k; }
CFX_ALWAYS_INLINE void cfx_u128_shr_eq(cfx_u128_t* x, unsigned k) { *x >>= k; }

/* Pure operations (replace one-liners): */
CFX_ALWAYS_INLINE cfx_u128_t cfx_u128_add(cfx_u128_t a, cfx_u128_t b) { return a + b; }
CFX_ALWAYS_INLINE cfx_u128_t cfx_u128_sub(cfx_u128_t a, cfx_u128_t b) { return a - b; }
CFX_ALWAYS_INLINE cfx_u128_t cfx_u128_shl(cfx_u128_t a, unsigned k) { return a << k; }
CFX_ALWAYS_INLINE cfx_u128_t cfx_u128_shr(cfx_u128_t a, unsigned k) { return a >> k; }

/* Multiply-accumulate (replace: acc += (u128)a*b;) */
CFX_ALWAYS_INLINE void cfx_u128_mac_u64(cfx_u128_t* acc, cfx_u64_t a, cfx_u64_t b) {
  *acc += (cfx_u128_t)a * (cfx_u128_t)b;
}

/* Wide multiply (replace: (__uint128_t)a*b with explicit hi/lo if you want) */
CFX_ALWAYS_INLINE cfx_u128_t cfx_u128_mul_u64(cfx_u64_t a, cfx_u64_t b) {
  return (cfx_u128_t)a * (cfx_u128_t)b;
}

/* Comparisons (replace: a < b, etc., in rare spots you want function form) */
CFX_ALWAYS_INLINE int cfx_u128_eq(cfx_u128_t a, cfx_u128_t b) { return a == b; }
CFX_ALWAYS_INLINE int cfx_u128_lt(cfx_u128_t a, cfx_u128_t b) { return a <  b; }
CFX_ALWAYS_INLINE int cfx_u128_le(cfx_u128_t a, cfx_u128_t b) { return a <= b; }
CFX_ALWAYS_INLINE int cfx_u128_gt(cfx_u128_t a, cfx_u128_t b) { return a >  b; }
CFX_ALWAYS_INLINE int cfx_u128_ge(cfx_u128_t a, cfx_u128_t b) { return a >= b; }

#endif /* CFX_HAS_UINT128 */

/* =========================
 *  64-bit helpers (native)
 * ========================= */
#if CFX_HAS_UINT64

/* Basic add/sub (optionally return carry/borrow when useful) */
CFX_ALWAYS_INLINE cfx_u64_t cfx_u64_add(cfx_u64_t a, cfx_u64_t b) { return a + b; }
CFX_ALWAYS_INLINE cfx_u64_t cfx_u64_sub(cfx_u64_t a, cfx_u64_t b) { return a - b; }
CFX_ALWAYS_INLINE cfx_u64_t cfx_u64_shl(cfx_u64_t a, unsigned k) { return a << k; }
CFX_ALWAYS_INLINE cfx_u64_t cfx_u64_shr(cfx_u64_t a, unsigned k) { return a >> k; }

CFX_ALWAYS_INLINE unsigned cfx_u64_add_overflow(cfx_u64_t a, cfx_u64_t b, cfx_u64_t* out) {
  cfx_u64_t s = a + b; *out = s; return s < a;
}
CFX_ALWAYS_INLINE unsigned cfx_u64_sub_borrow(cfx_u64_t a, cfx_u64_t b, cfx_u64_t* out) {
  cfx_u64_t d = a - b; *out = d; return a < b;
}

/* 64x64->128 via native __uint128_t when available (still zero-cost) */
#if CFX_HAS_UINT128
CFX_ALWAYS_INLINE void cfx_u64_mul_wide(cfx_u64_t a, cfx_u64_t b, cfx_u64_t* hi, cfx_u64_t* lo) {
  cfx_u128_t p = (cfx_u128_t)a * (cfx_u128_t)b;
  *lo = (cfx_u64_t)p;
  *hi = (cfx_u64_t)(p >> 64);
}
#endif

#endif /* CFX_HAS_UINT64 */

#endif  /* CFX_PRIM_H */
