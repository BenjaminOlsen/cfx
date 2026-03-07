# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later

#
# Big Integer Module
#
# This module uses a flat file structure with self-guarding backends.
# Each backend file wraps itself in #ifdef guards, so we simply compile
# all .c files and let the preprocessor select the active implementation.
#
# Backend selection is controlled by:
#   - CFX_TARGET_* macros (for compute backends: portable, x86_64_bmi2, arm_cortex_m4)
#   - CFX_MEMORY_STATIC macro (for memory backends: dynamic vs static)
#
# See doc/TARGET_ARCHITECTURE.md for design rationale.
#

set(CFX_BIG_SOURCES
    # Core logic (split from big.c)
    ${CMAKE_CURRENT_SOURCE_DIR}/big_core.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_bitops.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_arith.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_mul.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_div.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_mont.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_algos.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_fmt.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_fac.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ntt.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sbig.c

    # Compute backends (self-guarding via #ifdef)
    ${CMAKE_CURRENT_SOURCE_DIR}/big_portable.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_x86_64_bmi2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_arm_cortex_m4.c

    # Memory backends (self-guarding via #ifdef)
    ${CMAKE_CURRENT_SOURCE_DIR}/big_mem_dynamic.c
    ${CMAKE_CURRENT_SOURCE_DIR}/big_mem_static.c

    PARENT_SCOPE
)
