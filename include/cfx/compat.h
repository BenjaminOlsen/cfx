/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_COMPAT_H
#define CFX_COMPAT_H

/**
 * Platform compatibility layer for POSIX/Windows portability.
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Include inttypes.h for PRIu64/PRIu32 format specifiers */
#include <inttypes.h>

/* Fallback definitions for bare-metal toolchains that don't define these */
#ifndef PRIu64
#  define PRIu64 "llu"
#endif
#ifndef PRIu32
#  define PRIu32 "u"
#endif
#ifndef PRId64
#  define PRId64 "lld"
#endif
#ifndef PRIx64
#  define PRIx64 "llx"
#endif

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  include <intrin.h>
#elif defined(__arm__) && !defined(__linux__) && !defined(__APPLE__)
/* Bare-metal ARM (Cortex-M, etc.) - no OS, no pthread */
#  define CFX_NO_THREADS 1
#  define CFX_NO_CLOCK 1
#else
#  include <pthread.h>
#  include <unistd.h>
#  include <time.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Threading
 * ============================================================ */

#ifdef _WIN32

typedef HANDLE cfx_thread_t;

typedef struct {
    void * (*func)(void *);
    void *arg;
} cfx_thread_wrapper_t;

static DWORD WINAPI cfx_thread_wrapper(void *param) {
    cfx_thread_wrapper_t *w = (cfx_thread_wrapper_t *)param;
    void * (*func)(void *) = w->func;
    void *arg = w->arg;
    free(w);
    func(arg);
    return 0;
}

static inline int cfx_thread_create(cfx_thread_t *thread, void * (*func)(void *), void *arg) {
    cfx_thread_wrapper_t *w = (cfx_thread_wrapper_t *)malloc(sizeof(*w));
    if (!w) return -1;
    w->func = func;
    w->arg = arg;
    *thread = CreateThread(NULL, 0, cfx_thread_wrapper, w, 0, NULL);
    if (*thread == NULL) {
        free(w);
        return -1;
    }
    return 0;
}

static inline int cfx_thread_join(cfx_thread_t thread, void **retval) {
    (void)retval; /* Windows doesn't easily support return values */
    DWORD result = WaitForSingleObject(thread, INFINITE);
    CloseHandle(thread);
    return (result == WAIT_OBJECT_0) ? 0 : -1;
}

#elif defined(CFX_NO_THREADS)

/* Bare-metal / freestanding - no threading support */
typedef int cfx_thread_t;

static inline int cfx_thread_create(cfx_thread_t *thread, void * (*func)(void *), void *arg) {
    (void)thread; (void)func; (void)arg;
    return -1; /* Threading not supported */
}

static inline int cfx_thread_join(cfx_thread_t thread, void **retval) {
    (void)thread; (void)retval;
    return -1;
}

#else /* POSIX */

typedef pthread_t cfx_thread_t;

static inline int cfx_thread_create(cfx_thread_t *thread, void * (*func)(void *), void *arg) {
    return pthread_create(thread, NULL, func, arg);
}

static inline int cfx_thread_join(cfx_thread_t thread, void **retval) {
    return pthread_join(thread, retval);
}

#endif

/* ============================================================
 * Mutexes
 * ============================================================ */

#ifdef _WIN32

typedef SRWLOCK cfx_mutex_t;
#define CFX_MUTEX_INITIALIZER SRWLOCK_INIT

static inline void cfx_mutex_init(cfx_mutex_t *mtx) {
    InitializeSRWLock(mtx);
}

static inline void cfx_mutex_lock(cfx_mutex_t *mtx) {
    AcquireSRWLockExclusive(mtx);
}

static inline void cfx_mutex_unlock(cfx_mutex_t *mtx) {
    ReleaseSRWLockExclusive(mtx);
}

static inline void cfx_mutex_destroy(cfx_mutex_t *mtx) {
    (void)mtx; /* SRWLOCK doesn't need destruction */
}

#elif defined(CFX_NO_THREADS)

/* Bare-metal / freestanding - no mutex support (single-threaded) */
typedef int cfx_mutex_t;
#define CFX_MUTEX_INITIALIZER 0

static inline void cfx_mutex_init(cfx_mutex_t *mtx) {
    (void)mtx;
}
static inline void cfx_mutex_lock(cfx_mutex_t *mtx) {
    (void)mtx;
}
static inline void cfx_mutex_unlock(cfx_mutex_t *mtx) {
    (void)mtx;
}
static inline void cfx_mutex_destroy(cfx_mutex_t *mtx) {
    (void)mtx;
}

#else /* POSIX */

typedef pthread_mutex_t cfx_mutex_t;
#define CFX_MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER

static inline void cfx_mutex_init(cfx_mutex_t *mtx) {
    pthread_mutex_init(mtx, NULL);
}

static inline void cfx_mutex_lock(cfx_mutex_t *mtx) {
    pthread_mutex_lock(mtx);
}

static inline void cfx_mutex_unlock(cfx_mutex_t *mtx) {
    pthread_mutex_unlock(mtx);
}

static inline void cfx_mutex_destroy(cfx_mutex_t *mtx) {
    pthread_mutex_destroy(mtx);
}

#endif

/* ============================================================
 * Atomic Operations (for int/long)
 * ============================================================ */

#ifdef _WIN32

static inline int cfx_atomic_load(volatile long *ptr) {
    return InterlockedCompareExchange(ptr, 0, 0);
}

static inline void cfx_atomic_store(volatile long *ptr, long val) {
    InterlockedExchange(ptr, val);
}

#else /* POSIX / GCC */

static inline int cfx_atomic_load(volatile int *ptr) {
    return __atomic_load_n(ptr, __ATOMIC_SEQ_CST);
}

static inline void cfx_atomic_store(volatile int *ptr, int val) {
    __atomic_store_n(ptr, val, __ATOMIC_SEQ_CST);
}

#endif

/* Portable atomic flag type */
#ifdef _WIN32
typedef volatile long cfx_atomic_int;
#else
typedef volatile int cfx_atomic_int;
#endif

/* ============================================================
 * CPU Count
 * ============================================================ */

static inline int cfx_cpu_count(void) {
#ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return (int)sysinfo.dwNumberOfProcessors;
#elif defined(CFX_NO_THREADS)
    return 1; /* Single-core embedded */
#else
    long n = sysconf(_SC_NPROCESSORS_ONLN);
    return (n > 0) ? (int)n : 1;
#endif
}

/* ============================================================
 * High-Resolution Time (nanoseconds)
 * ============================================================ */

static inline uint64_t cfx_time_ns(void) {
#ifdef _WIN32
    LARGE_INTEGER freq, count;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&count);
    /* Convert to nanoseconds: count * 1e9 / freq */
    return (uint64_t)((count.QuadPart * 1000000000ULL) / freq.QuadPart);
#elif defined(CFX_NO_CLOCK)
    /* Bare-metal: no clock support - return 0 or use SysTick if available */
    return 0;
#else
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (uint64_t)ts.tv_sec * 1000000000ULL + (uint64_t)ts.tv_nsec;
#endif
}

/* ============================================================
 * String Utilities
 * ============================================================ */

/* strndup: duplicate at most n characters of a string */
static inline char * cfx_strndup(const char *s, size_t n) {
    size_t len = 0;
    while (len < n && s[len] != '\0') len++;
    char *dup = (char *)malloc(len + 1);
    if (dup) {
        memcpy(dup, s, len);
        dup[len] = '\0';
    }
    return dup;
}

/* ============================================================
 * Bit Operations
 * ============================================================ */

/* Count leading zeros for 32-bit unsigned integer */
static inline int cfx_clz32(uint32_t x) {
    if (x == 0) return 32;
#ifdef _WIN32
    unsigned long idx;
    _BitScanReverse(&idx, x);
    return 31 - (int)idx;
#else
    return __builtin_clz(x);
#endif
}

/* Count leading zeros for 64-bit unsigned integer */
static inline int cfx_clz64(uint64_t x) {
    if (x == 0) return 64;
#ifdef _WIN32
    unsigned long idx;
#  ifdef _WIN64
    _BitScanReverse64(&idx, x);
    return 63 - (int)idx;
#  else
    /* 32-bit Windows: check high word first */
    uint32_t hi = (uint32_t)(x >> 32);
    if (hi != 0) {
        _BitScanReverse(&idx, hi);
        return 31 - (int)idx;
    }
    _BitScanReverse(&idx, (uint32_t)x);
    return 63 - (int)idx;
#  endif
#else
    return __builtin_clzll(x);
#endif
}

/* ============================================================
 * 128-bit Arithmetic (for MSVC without __uint128_t)
 * ============================================================ */

#ifdef _MSC_VER
typedef struct { uint64_t lo, hi; } cfx_uint128_t;

static inline cfx_uint128_t cfx_mul64(uint64_t a, uint64_t b) {
    cfx_uint128_t r;
    r.lo = _umul128(a, b, &r.hi);
    return r;
}

static inline cfx_uint128_t cfx_add128(cfx_uint128_t a, cfx_uint128_t b) {
    cfx_uint128_t r;
    r.lo = a.lo + b.lo;
    r.hi = a.hi + b.hi + (r.lo < a.lo);
    return r;
}

static inline uint64_t cfx_shr128_51(cfx_uint128_t x) {
    return (x.hi << 13) | (x.lo >> 51);
}

static inline uint64_t cfx_lo128(cfx_uint128_t x) {
    return x.lo;
}

static inline cfx_uint128_t cfx_u128_from_u64(uint64_t x) {
    cfx_uint128_t r;
    r.lo = x;
    r.hi = 0;
    return r;
}
#endif

#ifdef __cplusplus
}
#endif

#endif /* CFX_COMPAT_H */
