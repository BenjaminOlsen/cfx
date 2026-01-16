# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later

# Facade and common sources (always included)
set(BIG_COMMON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/big.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ntt.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sbig.c
)

# ============================================================================
# MEMORY BACKEND SOURCES (CFX_MEMORY_MODE axis)
# ============================================================================
message(STATUS "cfx/big: Selecting memory backend sources for mode '${CFX_MEMORY_MODE}'...")

cfx_find_mem_source("big" "init" BIG_MEM_INIT_SOURCE)
cfx_find_mem_source("big" "reserve" BIG_MEM_RESERVE_SOURCE)
cfx_find_mem_source("big" "free" BIG_MEM_FREE_SOURCE)

set(BIG_MEM_SOURCES
    ${BIG_MEM_INIT_SOURCE}
    ${BIG_MEM_RESERVE_SOURCE}
    ${BIG_MEM_FREE_SOURCE}
)

# ============================================================================
# COMPUTE BACKEND SOURCES (CFX_TARGET axis)
# ============================================================================
# The cfx_find_backend_source function walks up the target hierarchy
# to find implementations, falling back to portable/ if needed.

message(STATUS "cfx/big: Selecting compute backend sources for target '${CFX_TARGET}'...")

cfx_find_backend_source("big" "mul" BIG_MUL_SOURCE)
cfx_find_backend_source("big" "add" BIG_ADD_SOURCE)
cfx_find_backend_source("big" "mont" BIG_MONT_SOURCE)

set(BIG_COMPUTE_SOURCES
    ${BIG_MUL_SOURCE}
    ${BIG_ADD_SOURCE}
    ${BIG_MONT_SOURCE}
)

set(CFX_BIG_SOURCES
    ${BIG_COMMON_SOURCES}
    ${BIG_MEM_SOURCES}
    ${BIG_COMPUTE_SOURCES}
    PARENT_SCOPE
)
