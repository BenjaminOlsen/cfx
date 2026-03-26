/* SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later */

#ifndef CFX_MEM_BACKEND_H
#define CFX_MEM_BACKEND_H

/*
 * Memory Backend Interface
 *
 * This header declares the memory management functions that vary based on
 * CFX_MEMORY_MODE (dynamic, static, custom). The facade (big.c) calls these
 * functions for all memory operations.
 *
 * Memory mode is orthogonal to CPU target:
 *   - CFX_TARGET selects compute backends (x86_64_bmi2, arm_neon, etc.)
 *   - CFX_MEMORY_MODE selects memory backends (dynamic, static)
 *
 * Implementations are in:
 *   - src/big/mem/dynamic/  (uses malloc/realloc/free)
 *   - src/big/mem/static/   (uses fixed-size buffers, no heap)
 */

#include "cfx/big.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * MEMORY BACKEND INTERFACE
 *
 * The memory backend functions are declared in cfx/big.h:
 *   - cfx_big_init()    : Initialize to empty state
 *   - cfx_big_reserve() : Ensure capacity for N limbs
 *   - cfx_big_free()    : Release resources
 *
 * Implementations are provided by the selected memory backend
 * (dynamic or static) in src/big/mem/MODE/
 */

/* MEMORY MODE CONFIGURATION */

/*
 * Static mode configuration
 *
 * When CFX_MEMORY_MODE=static, each cfx_big_t has a fixed-size buffer.
 * CFX_STATIC_LIMBS controls the maximum number of limbs per big integer.
 * CFX_STATIC_POOL_SIZE controls the number of buffers in the pool.
 *
 * These defaults are sized for running the test suite. For embedded use,
 * configure with smaller values via CMake:
 *   cmake -DCFX_STATIC_LIMBS=64 -DCFX_STATIC_POOL_SIZE=32 ..
 *
 * Defaults:
 *   - 256 limbs = 16384 bits (64-bit limbs) or 8192 bits (32-bit limbs)
 *   - 1024 buffers = up to 1024 concurrent big integers
 */
#ifndef CFX_STATIC_LIMBS
#define CFX_STATIC_LIMBS 256
#endif

#ifndef CFX_STATIC_POOL_SIZE
#define CFX_STATIC_POOL_SIZE 1024
#endif

/*
 * Memory mode identification for debugging/logging
 */
#if defined(CFX_MEMORY_STATIC)
    #define CFX_MEMORY_MODE_NAME "static"
#else
    #define CFX_MEMORY_MODE_NAME "dynamic"
#endif

#ifdef __cplusplus
}
#endif

#endif /* CFX_MEM_BACKEND_H */
