# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later

# Common crypto sources (always included)
set(CRYPTO_COMMON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/sha256.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sha512.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sha3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/blake2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/chacha20/chacha20.c
    ${CMAKE_CURRENT_SOURCE_DIR}/poly1305/poly1305.c
    ${CMAKE_CURRENT_SOURCE_DIR}/aead_chacha20_poly1305.c
    ${CMAKE_CURRENT_SOURCE_DIR}/aead_chacha20_poly1271.c
    ${CMAKE_CURRENT_SOURCE_DIR}/siphash.c
    ${CMAKE_CURRENT_SOURCE_DIR}/hmac.c
    ${CMAKE_CURRENT_SOURCE_DIR}/base64.c
)

# Curve25519/Ed25519 requires 128-bit arithmetic (64-bit platforms only)
# For 32-bit targets, this would need a 32-bit radix implementation
if(NOT CFX_TARGET MATCHES "arm_cortex_m4|arm_neon")
    list(APPEND CRYPTO_COMMON_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/fe25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/ge25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/sc25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/x25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/ed25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/argon2.c
        ${CMAKE_CURRENT_SOURCE_DIR}/poly1271/poly1271.c
        ${CMAKE_CURRENT_SOURCE_DIR}/poly1271/x86_64_avx2/poly1271_avx2.c
    )
else()
    message(STATUS "cfx/crypto: Excluding curve25519/ed25519 (requires 64-bit)")
endif()

# Select backend sources for algorithms with target-specific implementations
message(STATUS "cfx/crypto: Selecting backend sources for target '${CFX_TARGET}'...")

cfx_find_backend_source("crypto/chacha20" "block" CHACHA20_BLOCK_SOURCE)
cfx_find_backend_source("crypto/chacha20" "block4" CHACHA20_BLOCK4_SOURCE)
cfx_find_backend_source("crypto/chacha20" "block8" CHACHA20_BLOCK8_SOURCE)
cfx_find_backend_source("crypto/poly1305" "blocks" POLY1305_BLOCKS_SOURCE)

set(CRYPTO_BACKEND_SOURCES
    ${CHACHA20_BLOCK_SOURCE}
    ${CHACHA20_BLOCK4_SOURCE}
    ${CHACHA20_BLOCK8_SOURCE}
    ${POLY1305_BLOCKS_SOURCE}
)

set(CFX_CRYPTO_SOURCES
    ${CRYPTO_COMMON_SOURCES}
    ${CRYPTO_BACKEND_SOURCES}
    PARENT_SCOPE
)
