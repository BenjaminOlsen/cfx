cmake_minimum_required(VERSION 3.16)
project(cfx LANGUAGES C VERSION 0.0.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

option(CFX_BUILD_SHARED "Build shared library" ON)
option(CFX_BUILD_TESTS  "Build tests" ON)
option(CFX_BUILD_EXAMPLES "Build examples" ON)

set(CFX_HEADERS
    include/cfx/cfx.h
    include/cfx/big.h
    include/cfx/fac.h
    include/cfx/algo.h
    include/cfx/rand.h
    include/cfx/error.h
    include/cfx/primes.h
    include/cfx/vector.h
    include/cfx/macros.h
    include/cfx/memory.h
    include/cfx/version.h
    include/cfx/poly1305.h
    include/cfx/chacha20.h
    include/cfx/numerical.h
)

set(CFX_SOURCES
    src/big.c
    src/fac.c
    src/algo.c
    src/rand.c
    src/error.c
    src/memory.c
    src/primes.c
    src/vector.c
    src/chacha20.c
    src/poly1305.c
)

add_library(cfx ${CFX_SOURCES})
if(CFX_BUILD_SHARED)
  set_target_properties(cfx PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
  )
endif()

target_include_directories(cfx
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

if (UNIX AND NOT APPLE)
  target_link_libraries(cfx PUBLIC m)
  message(STATUS "Linking to libm")
endif()

include(cmake/cfxCompileFlags.cmake)

add_library(cfx::cfx ALIAS cfx)

# ---------------- Examples
if(CFX_BUILD_EXAMPLES)

  set(EXAMPLE_NAMES
    ulam_spiral
    cfx_keygen
    find_prime
    prime_test
    factorial
    factor
    dcmini
    choose
    encrypt
    primes
    fread
    mac
    div
    pow
  )
  
  add_library(example_utils examples/utils.h examples/utils.c)

  foreach(example_name IN LISTS EXAMPLE_NAMES)
    add_executable(${example_name} examples/${example_name}.c)
    target_link_libraries(${example_name} PRIVATE cfx::cfx example_utils)
  endforeach()
endif()

# ---------------- Tests (CTest)
include(CTest)

if(CFX_BUILD_TESTS)
  set(TEST_NAMES
    test_big
    test_big_mont
    test_fac
    test_algo
    test_pollard_rho
    test_mul_rows
    test_poly1305
    test_chacha20
    test_mul_csa_pthreads
  )

  foreach(test_name IN LISTS TEST_NAMES)
    add_executable(${test_name} test/cfx/${test_name}.c)
    target_link_libraries(${test_name} PRIVATE cfx::cfx)
    add_test(NAME ${test_name} COMMAND ${test_name})
  endforeach()
endif()

# ---------------- TestU01 testing:

option(CFX_ENABLE_TESTU01 "Build RNG tests with TestU01" OFF)
set(TESTU01_ROOT "" CACHE PATH "Path to TestU01 installation root")

if(CFX_ENABLE_TESTU01)
    if(NOT TESTU01_ROOT)
        message(FATAL_ERROR "Set TESTU01_ROOT to TestU01 install prefix")
    endif()

    include_directories(${TESTU01_ROOT}/include)

    add_executable(test_chacha20_testu01 test/testu01/test_chacha20_testu01.c)

    target_link_libraries(test_chacha20_testu01 PRIVATE
        cfx
        ${TESTU01_ROOT}/lib/libtestu01.a
        ${TESTU01_ROOT}/lib/libprobdist.a
        ${TESTU01_ROOT}/lib/libmylib.a
        m
    )
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

# Package config (find_package)
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

# ---------------- benchmarking:::::
option(CFX_BUILD_BENCHMARKS "Build microbenchmarks" OFF)
if (CFX_BUILD_BENCHMARKS)
    add_subdirectory(benchmark)
endif()
