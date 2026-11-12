# Constant-time test (no external dependencies)
# test_ct uses ed25519/x25519 which require 128-bit types (64-bit only)
if(NOT CFX_TARGET STREQUAL "arm_cortex_m4" AND NOT CFX_LIMB_BITS STREQUAL "32")
    add_executable(test_ct test_ct.c)
    target_link_libraries(test_ct PRIVATE cfx::cfx)
    target_compile_options(test_ct PRIVATE ${CFX_CFLAGS})
    message(STATUS "added stats test 'test_ct'")
else()
    message(STATUS "skipping stats test 'test_ct' (requires 64-bit)")
endif()

# TestU01 tests (optional)
option(CFX_ENABLE_TESTU01 "Build statistical tests using TestU01" OFF)

if(CFX_ENABLE_TESTU01)
    if(NOT TESTU01_ROOT)
        # check common install locations
        find_path(_testu01_inc NAMES unif01.h
            PATHS "$ENV{HOME}/libs/TestU01/include"
                  "/usr/local/include"
                  "/usr/include"
        )
        if(_testu01_inc)
            get_filename_component(TESTU01_ROOT "${_testu01_inc}" DIRECTORY)
        else()
            message(FATAL_ERROR
                "CFX_ENABLE_TESTU01=ON but TestU01 not found. "
                "Set -DTESTU01_ROOT=/path/to/TestU01")
        endif()
    endif()

    find_path(TESTU01_INCLUDE_DIR NAMES unif01.h
        PATHS "${TESTU01_ROOT}/include" NO_DEFAULT_PATH REQUIRED)
    foreach(_lib testu01 probdist mylib)
        find_library(TESTU01_${_lib}_LIBRARY NAMES ${_lib}
            PATHS "${TESTU01_ROOT}/lib" NO_DEFAULT_PATH REQUIRED)
    endforeach()
    message(STATUS "TestU01 enabled: ${TESTU01_ROOT}")

    set(CFX_STATS_TESTS
        test_testu01
    )

    foreach(test_name IN LISTS CFX_STATS_TESTS)
        add_executable(${test_name} ${test_name}.c hash_rngs.c)
        target_include_directories(${test_name} PRIVATE "${TESTU01_INCLUDE_DIR}")

        target_link_libraries(${test_name} PRIVATE
            cfx::cfx
            ${TESTU01_testu01_LIBRARY}
            ${TESTU01_probdist_LIBRARY}
            ${TESTU01_mylib_LIBRARY}
        )
        if(NOT MSVC)
            target_link_libraries(${test_name} PRIVATE m)
        endif()

        message(STATUS "added stats test '${test_name}'")
    endforeach()

    set(TESTS_USING_OPENSSL
        test_testu01
    )

    foreach(test_name IN LISTS TESTS_USING_OPENSSL)
        if (OpenSSL_FOUND)
            message(STATUS "Enabling OpenSSL for ${test_name}")
            target_link_libraries(${test_name} PRIVATE OpenSSL::Crypto)
        else()
            message(STATUS "Not using OpenSSL for ${test_name}")
        endif()
    endforeach()
else()
    message(STATUS "TestU01 not enabled (set -DCFX_ENABLE_TESTU01=ON to enable)")
endif()
