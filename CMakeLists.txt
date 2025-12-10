cmake_minimum_required(VERSION 3.16)
project(cfx LANGUAGES C VERSION 0.0.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

option(CFX_BUILD_SHARED "Build shared library" ON)
option(CFX_BUILD_TESTS  "Build tests" ON)
option(CFX_BUILD_EXAMPLES "Build examples" ON)
option(CFX_NO_FP "No floating point math" OFF)

set(CFX_ARCH "native" CACHE STRING "Target architecture")
set_property(CACHE CFX_ARCH PROPERTY STRINGS native x86_64 armv7m)

find_package(OpenSSL QUIET)

if(OpenSSL_FOUND)
    add_compile_definitions(CFX_HAVE_OPENSSL=1)
endif()

include(cmake/cfxCompileFlags.cmake)

# ---------------- static prime list
# specify with cmake -S . -B build -DCFX_PRIMES_COUNT=XXXX
set(CFX_PRIMES_COUNT 664579 CACHE STRING
    "Number of primes to include in cfx_primes[] (<= 664579)")

include(cmake/cfxPrimes.cmake)

set(CFX_PRIMES_TXT "${CMAKE_CURRENT_SOURCE_DIR}/data/primes.txt")
set(CFX_PRIMES_C   "${CMAKE_CURRENT_BINARY_DIR}/cfx_primes.c")

cfx_generate_primes("${CFX_PRIMES_TXT}"
                    "${CFX_PRIMES_C}"
                    ${CFX_PRIMES_COUNT})

set(CFX_HEADERS
    include/cfx/cfx.h
    include/cfx/cpu.h
    include/cfx/big.h
    include/cfx/fac.h
    include/cfx/algo.h
    include/cfx/rand.h
    include/cfx/primes.h
    include/cfx/vec.h
    include/cfx/big_vec.h
    include/cfx/macros.h
    include/cfx/memory.h
    include/cfx/version.h
    include/cfx/sha256.h
    include/cfx/poly1305.h
    include/cfx/chacha20.h
    include/cfx/arith.h
    include/cfx/aead_chacha20_poly1305.h
)

set(CFX_SOURCES
    src/big.c
    src/fac.c
    src/algo.c
    src/rand.c
    src/vec.c
    src/big_vec.c
    src/sha256.c
    src/chacha20.c
    src/poly1305.c
    src/aead_chacha20_poly1305.c
    "${CFX_PRIMES_GENERATED_SRC}"
)

add_library(cfx ${CFX_SOURCES})
add_library(cfx::cfx ALIAS cfx)

cfx_apply_compile_flags(cfx)

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

if (CFX_NO_FP)
    target_compile_definitions(cfx PUBLIC CFX_NO_FP)
    message(STATUS "cfx: No floating point math")
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

# ---------------- Examples
if(CFX_BUILD_EXAMPLES)
    message(STATUS "Building examples...")
    add_subdirectory(examples)
endif()

# ---------------- Tests (CTest)
include(CTest)

if(CFX_BUILD_TESTS)
    message(STATUS "Building tests...")
    add_subdirectory(test)
endif()

# ---------------- Limb bit selection
# use with -DCFX_LIMB_BITS=64 or 32
set(CFX_LIMB_BITS "64" CACHE STRING "Limb size for cfx: 32 or 64")
# set_property(CACHE CFX_LIMB_BITS PROPERTY STRINGS 32 64)
if (CFX_LIMB_BITS STREQUAL "32")
    add_compile_definitions(CFX_FORCE_LIMB_32)
    message(STATUS "using 32 bit limbs")
elseif (CFX_LIMB_BITS STREQUAL "64")
    add_compile_definitions(CFX_FORCE_LIMB_64)
    message(STATUS "using 64 bit limbs")
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
  cmake/cfxConfig.cmake.in
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
