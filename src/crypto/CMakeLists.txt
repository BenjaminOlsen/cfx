# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later

#
# Cryptographic Primitives Module
#
# Uses flat file structure with self-guarding backends (like src/big/).
# Each backend file wraps itself in #ifdef guards, so we simply compile
# all .c files and let the preprocessor select the active implementation.
#
# See doc/TARGET_ARCHITECTURE.md Appendix B for design rationale.
#

# Common crypto sources (always included)
set(CRYPTO_COMMON_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/sha256.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sha512.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sha3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/blake2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/poly1305.c
    ${CMAKE_CURRENT_SOURCE_DIR}/aead_chacha20_poly1305.c
    ${CMAKE_CURRENT_SOURCE_DIR}/siphash.c
    ${CMAKE_CURRENT_SOURCE_DIR}/base64.c
)

# ChaCha20 facade and backends (self-guarding)
set(CHACHA20_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/chacha20/chacha20.c
    ${CMAKE_CURRENT_SOURCE_DIR}/chacha20/chacha20_portable.c
    ${CMAKE_CURRENT_SOURCE_DIR}/chacha20/chacha20_x86_64_avx2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/chacha20/chacha20_arm_cortex_m4.c
)

# Poly1305 backends (self-guarding)
set(POLY1305_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/poly1305/poly1305_portable.c
    ${CMAKE_CURRENT_SOURCE_DIR}/poly1305/poly1305_arm_cortex_m4.c
)

# Curve25519/Ed25519/Poly1271 require 128-bit arithmetic (64-bit platforms only)
# For 32-bit targets (arm_cortex_m4, or builds with 32-bit limbs), exclude these
if(NOT CFX_TARGET STREQUAL "arm_cortex_m4" AND NOT CFX_LIMB_BITS STREQUAL "32")
    list(APPEND CRYPTO_COMMON_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/fe25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/ge25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/sc25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/x25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/ed25519.c
        ${CMAKE_CURRENT_SOURCE_DIR}/argon2.c
        ${CMAKE_CURRENT_SOURCE_DIR}/poly1271/poly1271.c
        ${CMAKE_CURRENT_SOURCE_DIR}/poly1271/poly1271_avx2.c
        ${CMAKE_CURRENT_SOURCE_DIR}/aead_chacha20_poly1271.c
    )
else()
    message(STATUS "cfx/crypto: Excluding curve25519/ed25519/poly1271 (requires 64-bit)")
endif()

set(CFX_CRYPTO_SOURCES
    ${CRYPTO_COMMON_SOURCES}
    ${CHACHA20_SOURCES}
    ${POLY1305_SOURCES}
    PARENT_SCOPE
)
