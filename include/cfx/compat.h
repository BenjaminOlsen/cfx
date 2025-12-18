/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_COMPAT_H
#define CFX_COMPAT_H

/**
 * Platform compatibility layer for POSIX/Windows portability.
 * Provides unified APIs for:
 *   - Threading (create, join)
 *   - Mutexes
 *   - Atomic operations
 *   - CPU count detection
 *   - High-resolution time
 *   - String utilities (strndup)
 *   - Bit operations (clz)
 */

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#  define WIN32_LEAN_AND_MEAN
#  include <windows.h>
#  include <intrin.h>
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
typedef DWORD (WINAPI *cfx_thread_func_t)(void*);

static inline int cfx_thread_create(cfx_thread_t* thread, void* (*func)(void*), void* arg) {
    /* Windows thread func returns DWORD, POSIX returns void*.
     * We cast and ignore return value differences. */
    *thread = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, arg, 0, NULL);
    return (*thread != NULL) ? 0 : -1;
}

static inline int cfx_thread_join(cfx_thread_t thread, void** retval) {
    (void)retval; /* Windows doesn't easily support return values */
    DWORD result = WaitForSingleObject(thread, INFINITE);
    CloseHandle(thread);
    return (result == WAIT_OBJECT_0) ? 0 : -1;
}

#else /* POSIX */

typedef pthread_t cfx_thread_t;

static inline int cfx_thread_create(cfx_thread_t* thread, void* (*func)(void*), void* arg) {
    return pthread_create(thread, NULL, func, arg);
}

static inline int cfx_thread_join(cfx_thread_t thread, void** retval) {
    return pthread_join(thread, retval);
}

#endif

/* ============================================================
 * Mutexes
 * ============================================================ */

#ifdef _WIN32

typedef SRWLOCK cfx_mutex_t;
#define CFX_MUTEX_INITIALIZER SRWLOCK_INIT

static inline void cfx_mutex_init(cfx_mutex_t* mtx) {
    InitializeSRWLock(mtx);
}

static inline void cfx_mutex_lock(cfx_mutex_t* mtx) {
    AcquireSRWLockExclusive(mtx);
}

static inline void cfx_mutex_unlock(cfx_mutex_t* mtx) {
    ReleaseSRWLockExclusive(mtx);
}

static inline void cfx_mutex_destroy(cfx_mutex_t* mtx) {
    (void)mtx; /* SRWLOCK doesn't need destruction */
}

#else /* POSIX */

typedef pthread_mutex_t cfx_mutex_t;
#define CFX_MUTEX_INITIALIZER PTHREAD_MUTEX_INITIALIZER

static inline void cfx_mutex_init(cfx_mutex_t* mtx) {
    pthread_mutex_init(mtx, NULL);
}

static inline void cfx_mutex_lock(cfx_mutex_t* mtx) {
    pthread_mutex_lock(mtx);
}

static inline void cfx_mutex_unlock(cfx_mutex_t* mtx) {
    pthread_mutex_unlock(mtx);
}

static inline void cfx_mutex_destroy(cfx_mutex_t* mtx) {
    pthread_mutex_destroy(mtx);
}

#endif

/* ============================================================
 * Atomic Operations (for int/long)
 * ============================================================ */

#ifdef _WIN32

static inline int cfx_atomic_load(volatile long* ptr) {
    return InterlockedCompareExchange(ptr, 0, 0);
}

static inline void cfx_atomic_store(volatile long* ptr, long val) {
    InterlockedExchange(ptr, val);
}

#else /* POSIX / GCC */

static inline int cfx_atomic_load(volatile int* ptr) {
    return __atomic_load_n(ptr, __ATOMIC_SEQ_CST);
}

static inline void cfx_atomic_store(volatile int* ptr, int val) {
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
static inline char* cfx_strndup(const char* s, size_t n) {
#if defined(_WIN32) || !defined(_GNU_SOURCE)
    size_t len = 0;
    while (len < n && s[len] != '\0') len++;
    char* dup = (char*)malloc(len + 1);
    if (dup) {
        memcpy(dup, s, len);
        dup[len] = '\0';
    }
    return dup;
#else
    return strndup(s, n);
#endif
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

#ifdef __cplusplus
}
#endif

#endif /* CFX_COMPAT_H */
