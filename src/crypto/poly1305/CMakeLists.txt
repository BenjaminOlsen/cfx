# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later

# Poly1305 backend source selection
# Uses the cfx_find_backend_source function from cfx_target.cmake
# to walk up the target hierarchy and find implementations.

message(STATUS "cfx/crypto/poly1305: Selecting backend sources for target '${CFX_TARGET}'...")

cfx_find_backend_source("crypto/poly1305" "blocks" POLY1305_BLOCKS_SOURCE)

set(CFX_POLY1305_BACKEND_SOURCES
    ${POLY1305_BLOCKS_SOURCE}
    PARENT_SCOPE
)
