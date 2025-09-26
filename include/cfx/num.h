/* cfx_num.h : C compile-time polymorphism for limb/accumulator
*
* Configure (optional):
*   -DCFX_FORCE_LIMB_32    or  -DCFX_FORCE_LIMB_64
*   -DCFX_FORCE_ACC_STRUCT (use {lo,hi} even if native wide exists)
*   -DCFX_FORCE_NO_UINT128 (pretend no __uint128_t)
*/
#ifndef CFX_NUM_H
#define CFX_NUM_H

#include <stdint.h>

/* -------- feature detection -------- */
#if !defined(CFX_FORCE_NO_UINT128) && defined(__SIZEOF_INT128__)
  #define CFX_HAS_UINT128 1
#endif

/* -------- always_inline helper -------- */
#if defined(_MSC_VER)
  #define CFX_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__)
  #define CFX_INLINE __attribute__((always_inline)) static inline
#else
  #define CFX_INLINE static inline
#endif

/* -------- limb selection: 32 or 64 -------- */
#if defined(CFX_FORCE_LIMB_32)
  typedef uint32_t cfx_limb_t;
  #define CFX_LIMB_BITS 32
#elif defined(CFX_FORCE_LIMB_64)
  typedef uint64_t cfx_limb_t;
  #define CFX_LIMB_BITS 64
#else
  /* Auto: prefer 64-bit limbs when native 128 exists; otherwise 32-bit. */
  #if CFX_HAS_UINT128
    typedef uint64_t cfx_limb_t;
    #define CFX_LIMB_BITS 64
  #else
    typedef uint32_t cfx_limb_t;
    #define CFX_LIMB_BITS 32
  #endif
#endif

/* -------- print formats, sqrt max etc -------- */
#if (CFX_LIMB_BITS == 64)
  #define CFX_SQRT_ACC_MAX 0xFFFFFFFFFFFFFFFFllu   /* floor(sqrt(UINT128_MAX+1)) - 1 */
  #define CFX_PRIxLIMB "%" PRIx64
  #define CFX_PRIuLIMB "%" PRIu64
  #define CFX_PRI0xLIMB "%016" PRIx64
  #define CFX_PRI0uLIMB "%016" PRIu64
#elif (CFX_LIMB_BITS == 32)
  #define CFX_SQRT_ACC_MAX 0xFFFFFFFFllu           /* floor(sqrt(UINT64_MAX+1)) - 1 */
  #define CFX_PRIxLIMB "%" PRIx32
  #define CFX_PRIuLIMB "%" PRIu32
  #define CFX_PRI0xLIMB "%08" PRIx32
  #define CFX_PRI0uLIMB "%08" PRIu32
#else
/* todo - choose defaults or error  */
#endif 


/* -------- accumulator selection: 2x limb width --------
 * If limb=64 and native 128 exists and not forced to struct -> use __uint128_t.
 * Else use a portable {lo,hi} struct with limb-size halves.
 */
#if (CFX_LIMB_BITS == 64) && CFX_HAS_UINT128 && !defined(CFX_FORCE_ACC_STRUCT)
  typedef __uint128_t cfx_acc_t;
  #define CFX_ACC_NATIVE 1
#else
  typedef struct { cfx_limb_t lo, hi; } cfx_acc_t;
  #define CFX_SQRT_ACC_MAX 0xFFFFFFFFllu           /* floor(sqrt(UINT64_MAX+1)) - 1 */
#endif

/* -------- API: constructors / accessors -------- */
CFX_INLINE cfx_acc_t cfx_acc_zero(void) {
#if defined(CFX_ACC_NATIVE)
  return (cfx_acc_t)0;
#else
  cfx_acc_t a; a.lo = 0; a.hi = 0; return a;
#endif
}

CFX_INLINE cfx_acc_t cfx_acc_from_lo(cfx_limb_t lo) {
#if defined(CFX_ACC_NATIVE)
  return (cfx_acc_t)lo;
#else
  cfx_acc_t a; a.lo = lo; a.hi = 0; return a;
#endif
}

/* Low/High extractors (limb-sized) */
CFX_INLINE cfx_limb_t cfx_acc_lo(cfx_acc_t a) {
#if defined(CFX_ACC_NATIVE)
  return (cfx_limb_t)(a);
#else
  return a.lo;
#endif
}

CFX_INLINE cfx_limb_t cfx_acc_hi(cfx_acc_t a) {
#if defined(CFX_ACC_NATIVE)
  return (cfx_limb_t)(a >> CFX_LIMB_BITS);
#else
  return a.hi;
#endif
}

/* -------- primitive: limb×limb -> {hi,lo} (each limb-sized) -------- */
CFX_INLINE void cfx_mul_wide(cfx_limb_t x, cfx_limb_t y,
                             cfx_limb_t* hi, cfx_limb_t* lo)
{
#if (CFX_LIMB_BITS == 64) && CFX_HAS_UINT128
  /* 64×64→128: use native when present (even if acc is struct) */
  __uint128_t p = ( (__uint128_t)x * (__uint128_t)y );
  *lo = (cfx_limb_t)p;
  *hi = (cfx_limb_t)(p >> CFX_LIMB_BITS);
#elif (CFX_LIMB_BITS == 64)
  /* Portable 64×64→128 via 32-bit halves */
  const uint64_t x0 = (uint32_t)x, x1 = (uint64_t)x >> 32;
  const uint64_t y0 = (uint32_t)y, y1 = (uint64_t)y >> 32;

  uint64_t p00 = x0 * y0;
  uint64_t p01 = x0 * y1;
  uint64_t p10 = x1 * y0;
  uint64_t p11 = x1 * y1;

  uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
  *lo = (cfx_limb_t)((p00 & 0xffffffffu) | (mid << 32));
  *hi = (cfx_limb_t)(p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32));
#else /* CFX_LIMB_BITS == 32 */
  /* 32×32→64 fits in uint64_t portably */
  uint64_t p = (uint64_t)x * (uint64_t)y;
  *lo = (cfx_limb_t)p;
  *hi = (cfx_limb_t)(p >> 32);
#endif
}

/* -------- primitive: acc += low limb -------- */
CFX_INLINE void cfx_acc_add_lo(cfx_acc_t* acc, cfx_limb_t add)
{
#if defined(CFX_ACC_NATIVE)
  *acc += (cfx_acc_t)add;
#else
  cfx_limb_t old = acc->lo;
  acc->lo = old + add;
  acc->hi += (acc->lo < old); /* carry */
#endif
}

/* -------- primitive: acc += high limb (i.e., add at +LIMB_BITS) -------- */
CFX_INLINE void cfx_acc_add_hi(cfx_acc_t* acc, cfx_limb_t add_hi)
{
#if defined(CFX_ACC_NATIVE)
  *acc += ( (cfx_acc_t)add_hi << CFX_LIMB_BITS );
#else
  /* adding into the high half directly */
  cfx_limb_t old = acc->hi;
  acc->hi = old + add_hi;
  (void)old;
#endif
}

/* -------- primitive: acc += x*y -------- */
CFX_INLINE void cfx_acc_mac(cfx_acc_t* acc, cfx_limb_t x, cfx_limb_t y)
{
#if defined(CFX_ACC_NATIVE) && (CFX_LIMB_BITS == 64)
  *acc += ( (__uint128_t)x * (__uint128_t)y );
#else
  cfx_limb_t hi, lo;
  cfx_mul_wide(x, y, &hi, &lo);
  cfx_acc_add_lo(acc, lo);
  cfx_acc_add_hi(acc, hi);
#endif
}

#endif  /* CFX_NUM_H */
