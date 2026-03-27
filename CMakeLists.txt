cmake_minimum_required(VERSION 3.16)
project(cfx LANGUAGES C VERSION 0.0.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

option(CFX_BUILD_SHARED "Build shared library" ON)
option(CFX_BUILD_TESTS  "Build tests" ON)
option(CFX_BUILD_UTILS "Build utility programs" ON)
option(CFX_NO_FLOAT "No floating point math" OFF)
option(CFX_COVERAGE "Enable code coverage (GCC/Clang only)" OFF)
option(CFX_ED25519_PARANOID "Enable paranoid Ed25519 checks (small-subgroup rejection, cofactored verification)" OFF)

if(CFX_COVERAGE)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
        add_compile_options(--coverage -O0 -g)
        add_link_options(--coverage)
        message(STATUS "Code coverage enabled")
    else()
        message(WARNING "CFX_COVERAGE requires GCC or Clang")
    endif()
endif()

find_package(OpenSSL QUIET)

if(OpenSSL_FOUND)
    add_compile_definitions(CFX_HAVE_OPENSSL=1)
endif()

include(cmake/compile_flags.cmake)
include(cmake/cfx_target.cmake)
include(cmake/cfxQemuTest.cmake)

# Set up QEMU emulation for cross-compiled tests
cfx_setup_qemu_emulation()

# ---------------- static prime list
# specify with cmake -S . -B build -DCFX_PRIMES_COUNT=XXXX
set(CFX_PRIMES_COUNT 664579 CACHE STRING
    "Number of primes to include in cfx_primes[] (<= 664579)")

include(cmake/primes.cmake)

set(CFX_PRIMES_TXT "${CMAKE_CURRENT_SOURCE_DIR}/data/primes.txt")
set(CFX_PRIMES_C   "${CMAKE_CURRENT_BINARY_DIR}/cfx_primes.c")

cfx_generate_primes("${CFX_PRIMES_TXT}"
                    "${CFX_PRIMES_C}"
                    ${CFX_PRIMES_COUNT})

set(CFX_HEADERS
    include/cfx/aead_chacha20_poly1271.h
    include/cfx/aead_chacha20_poly1305.h
    include/cfx/algo.h
    include/cfx/arch.h
    include/cfx/argon2.h
    include/cfx/arith.h
    include/cfx/base64.h
    include/cfx/big.h
    include/cfx/bits.h
    include/cfx/blake2.h
    include/cfx/chacha20.h
    include/cfx/compat.h
    include/cfx/ecm.h
    include/cfx/ed25519.h
    include/cfx/fac.h
    include/cfx/fe25519.h
    include/cfx/ge25519.h
    include/cfx/sc25519.h
    include/cfx/macros.h
    include/cfx/memory.h
    include/cfx/mont.h
    include/cfx/ntt.h
    include/cfx/poly1271.h
    include/cfx/poly1305.h
    include/cfx/primes.h
    include/cfx/rand.h
    include/cfx/rng_chacha20.h
    include/cfx/rng_lcg.h
    include/cfx/rng_pcg.h
    include/cfx/rng_splitmix.h
    include/cfx/rng_xoroshiro.h
    include/cfx/rng_xorshift.h
    include/cfx/sha256.h
    include/cfx/sha3.h
    include/cfx/sha512.h
    include/cfx/siphash.h
    include/cfx/vec.h
    include/cfx/version.h
    include/cfx/sbig.h
    include/cfx/target.h
    include/cfx/x25519.h
)

# Source modules
add_subdirectory(src/core)
add_subdirectory(src/crypto)
add_subdirectory(src/rng)
add_subdirectory(src/big)
add_subdirectory(src/math)
set(CFX_SOURCES
    ${CFX_CORE_SOURCES}
    ${CFX_CRYPTO_SOURCES}
    ${CFX_RNG_SOURCES}
    ${CFX_BIG_SOURCES}
    ${CFX_MATH_SOURCES}
    "${CFX_PRIMES_GENERATED_SRC}"
)

add_library(cfx ${CFX_SOURCES})
add_library(cfx::cfx ALIAS cfx)

cfx_apply_compile_flags(cfx)
cfx_apply_target(cfx)
cfx_apply_memory_mode(cfx)

set_property(TARGET cfx PROPERTY C_STANDARD 99)
set_property(TARGET cfx PROPERTY C_STANDARD_REQUIRED ON)
set_property(TARGET cfx PROPERTY C_EXTENSIONS OFF)

if(CFX_BUILD_SHARED)
  set_property(TARGET cfx PROPERTY POSITION_INDEPENDENT_CODE ON)
  set_property(TARGET cfx PROPERTY WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

target_include_directories(cfx
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

include(CheckLibraryExists)
check_library_exists(m log "" HAVE_LIBM)

if (HAVE_LIBM)
    target_link_libraries(cfx PUBLIC m)
    message(STATUS "Linking to libm")
else()
    message(STATUS "Target has no libm")  # (OK for macOS or freestanding)
endif()

if (CFX_NO_FLOAT)
    target_compile_definitions(cfx PUBLIC CFX_NO_FLOAT)
    message(STATUS "cfx: No floating point math")
endif()

if (CFX_ED25519_PARANOID)
    target_compile_definitions(cfx PUBLIC CFX_ED25519_PARANOID)
    message(STATUS "cfx: Paranoid Ed25519 checks enabled (small subgroup, cofactored verification)")
endif()

# Feature test macros for POSIX extensions (clock_gettime, sysconf, etc.)
# Use _GNU_SOURCE on Linux (implies _POSIX_C_SOURCE + extensions)
# macOS doesn't need this - it exposes everything by default
if (NOT WIN32 AND NOT APPLE)
    target_compile_definitions(cfx PUBLIC _GNU_SOURCE)
endif()

# ---------------- ASAN
option(ASAN "Enable ASAN" OFF)

if (ASAN)
    message(STATUS "ASAN enabled")
    set(ASAN_FLAGS "-fsanitize=address -fno-omit-frame-pointer")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ASAN_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ASAN_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ASAN_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ASAN_FLAGS}")

    # macOS specific: use fast unwinder (much nicer stack traces)
    if (APPLE)
        set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} -fsanitize-address-use-after-return=always")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize-address-use-after-return=always")
    endif()
endif()

# ---------------- Utilities
if(CFX_BUILD_UTILS)
    message(STATUS "Building utilities...")
    add_subdirectory(utils)
endif()

# ---------------- Tests (CTest)
include(CTest)

if(CFX_BUILD_TESTS)
    message(STATUS "Building tests...")
    add_subdirectory(test)
endif()

# ---------------- Limb bit selection
# use with -DCFX_LIMB_BITS=auto, 32 or 64
# "auto" = 32 on MSVC (no __uint128_t), 64 on GCC/Clang (has __uint128_t)
# Note: Some targets (arm_cortex_m4) force 32-bit limbs via cfx_apply_target()
set(CFX_LIMB_BITS "auto" CACHE STRING "Limb size for cfx: auto, 32 or 64")
set_property(CACHE CFX_LIMB_BITS PROPERTY STRINGS auto 32 64)

# Check if target already forces limb size (e.g., arm_cortex_m4 -> 32-bit)
if(CFX_TARGET STREQUAL "arm_cortex_m4")
    # Already set by cfx_apply_target() - don't override
    message(STATUS "cfx: using 32-bit limbs (forced by target ${CFX_TARGET})")
elseif (CFX_LIMB_BITS STREQUAL "auto")
    if (MSVC)
        # MSVC has no __uint128_t, so force 32-bit limbs for 64-bit accumulator
        target_compile_definitions(cfx PUBLIC CFX_FORCE_LIMB_32)
        message(STATUS "cfx: using 32-bit limbs (auto for MSVC)")
    else()
        # GCC/Clang with __uint128_t can use 64-bit limbs
        target_compile_definitions(cfx PUBLIC CFX_FORCE_LIMB_64)
        message(STATUS "cfx: using 64-bit limbs (auto for GCC/Clang)")
    endif()
elseif (CFX_LIMB_BITS STREQUAL "32")
    target_compile_definitions(cfx PUBLIC CFX_FORCE_LIMB_32)
    message(STATUS "cfx: using 32-bit limbs")
elseif (CFX_LIMB_BITS STREQUAL "64")
    target_compile_definitions(cfx PUBLIC CFX_FORCE_LIMB_64)
    message(STATUS "cfx: using 64-bit limbs")
else()
    message(FATAL_ERROR "Invalid CFX_LIMB_BITS value: ${CFX_LIMB_BITS}")
endif()

# ---------------- Install
include(GNUInstallDirs)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX
      "${CMAKE_BINARY_DIR}/install"
      CACHE PATH "Install prefix" FORCE)
  message(STATUS "Install prefix set to default: ${CMAKE_INSTALL_PREFIX}")
else()
  message(STATUS "Install prefix set to: ${CMAKE_INSTALL_PREFIX}")
endif()

install(TARGETS cfx
  EXPORT cfxTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES ${CFX_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/cfx)

# ---------------- Package config (find_package)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/cfxConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cfx
)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/cfxConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(EXPORT cfxTargets
  NAMESPACE cfx::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cfx
)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/cfxConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/cfxConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/cfx
)

# ---------------- benchmarking
option(CFX_BUILD_BENCHMARKS "Build microbenchmarks" OFF)
if (CFX_BUILD_BENCHMARKS)
    message(STATUS "Building benchmarks...")
    add_subdirectory(benchmark)
endif()

