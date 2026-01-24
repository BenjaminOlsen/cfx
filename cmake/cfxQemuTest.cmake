# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later
#
# QEMU Test Emulation Support
#
# This module provides functions to run tests under QEMU emulation
# when cross-compiling for embedded targets like Cortex-M4.
#
# Usage:
#   include(cfxQemuTest)
#   cfx_setup_qemu_emulation()
#
# This will:
#   - Detect QEMU installation
#   - Set CMAKE_CROSSCOMPILING_EMULATOR for the target
#   - Configure CTest to use QEMU for test execution

include_guard(GLOBAL)

#
# Find QEMU and set up emulation for cross-compiled tests
#
function(cfx_setup_qemu_emulation)
    if(NOT CMAKE_CROSSCOMPILING)
        message(STATUS "cfx/qemu: Native build, skipping QEMU setup")
        return()
    endif()

    # Determine which QEMU to use based on target
    if(CFX_TARGET STREQUAL "arm_cortex_m4")
        set(_qemu_binary "qemu-arm")
        set(_qemu_cpu "cortex-m4")
        set(_qemu_args "-cpu" "${_qemu_cpu}")
    elseif(CFX_TARGET MATCHES "aarch64")
        set(_qemu_binary "qemu-aarch64")
        set(_qemu_cpu "cortex-a72")
        set(_qemu_args "-cpu" "${_qemu_cpu}")
    elseif(CFX_TARGET MATCHES "riscv64")
        set(_qemu_binary "qemu-riscv64")
        set(_qemu_args "")
    else()
        message(STATUS "cfx/qemu: No QEMU mapping for target '${CFX_TARGET}'")
        return()
    endif()

    # Find QEMU binary
    find_program(QEMU_EXECUTABLE ${_qemu_binary}
        HINTS
            /usr/bin
            /usr/local/bin
            /opt/homebrew/bin
            $ENV{QEMU_PATH}
    )

    if(NOT QEMU_EXECUTABLE)
        message(WARNING "cfx/qemu: ${_qemu_binary} not found. Tests will not run.")
        message(WARNING "cfx/qemu: Install QEMU or set QEMU_PATH environment variable.")
        return()
    endif()

    message(STATUS "cfx/qemu: Found ${_qemu_binary}: ${QEMU_EXECUTABLE}")
    message(STATUS "cfx/qemu: CPU emulation: ${_qemu_cpu}")

    # Set the crosscompiling emulator for CMake/CTest
    # This makes add_test() automatically prefix test commands with QEMU
    set(CMAKE_CROSSCOMPILING_EMULATOR "${QEMU_EXECUTABLE};${_qemu_args}" PARENT_SCOPE)

    # Export for use in custom scripts
    set(CFX_QEMU_EXECUTABLE "${QEMU_EXECUTABLE}" CACHE INTERNAL "")
    set(CFX_QEMU_CPU "${_qemu_cpu}" CACHE INTERNAL "")
    set(CFX_QEMU_ARGS "${_qemu_args}" CACHE INTERNAL "")

endfunction()

#
# Create a custom target to run a specific test under QEMU
#
function(cfx_add_qemu_test test_name test_executable)
    if(NOT CFX_QEMU_EXECUTABLE)
        return()
    endif()

    add_custom_target(qemu_${test_name}
        COMMAND ${CFX_QEMU_EXECUTABLE} ${CFX_QEMU_ARGS} $<TARGET_FILE:${test_executable}>
        DEPENDS ${test_executable}
        COMMENT "Running ${test_name} under QEMU (${CFX_QEMU_CPU})"
        VERBATIM
    )
endfunction()

#
# Create a target to run all tests under QEMU
#
function(cfx_add_qemu_test_all)
    if(NOT CFX_QEMU_EXECUTABLE)
        return()
    endif()

    add_custom_target(qemu_test_all
        COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure
        DEPENDS ${ARGN}
        COMMENT "Running all tests under QEMU (${CFX_QEMU_CPU})"
        VERBATIM
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )
endfunction()
