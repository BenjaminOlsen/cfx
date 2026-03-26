# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later
#
# CFX Target Selection System
#
# Usage:
#   cmake -DCFX_TARGET=x86_64_avx2 ..
#   cmake -DCFX_TARGET=arm_cortex_m4 -DCMAKE_TOOLCHAIN_FILE=... ..
#

include_guard(GLOBAL)

# Available targets (hierarchical)
set(CFX_TARGETS
    portable            # Portable C99, no intrinsics (root of inheritance tree)
    x86_64          # x86-64 portableline
    x86_64_bmi2     # x86-64 + BMI2 (mulx, adcx, adox)
    x86_64_avx2     # x86-64 + AVX2 + BMI2
    x86_64_avx512   # x86-64 + AVX-512
    arm_cortex_m4   # ARM Cortex-M4 (32-bit, DSP, no NEON)
    arm_neon        # ARMv7 + NEON
    aarch64_neon    # AArch64 (NEON is portableline)
)

# Target option with auto-detection default
set(CFX_TARGET "auto" CACHE STRING "Target architecture for optimized implementations")
set_property(CACHE CFX_TARGET PROPERTY STRINGS auto ${CFX_TARGETS})

# Target capability flags (populated by target selection)
set(CFX_CAP_BMI2 OFF CACHE INTERNAL "")
set(CFX_CAP_AVX2 OFF CACHE INTERNAL "")
set(CFX_CAP_AVX512 OFF CACHE INTERNAL "")
set(CFX_CAP_NEON OFF CACHE INTERNAL "")
set(CFX_CAP_DSP OFF CACHE INTERNAL "")

# Target parent relationships (for inheritance / fallback resolution)
# portable has no parent (it's the root)
set(CFX_TARGET_PARENT_x86_64 "portable")
set(CFX_TARGET_PARENT_x86_64_bmi2 "x86_64")
set(CFX_TARGET_PARENT_x86_64_avx2 "x86_64_bmi2")
set(CFX_TARGET_PARENT_x86_64_avx512 "x86_64_avx2")
set(CFX_TARGET_PARENT_arm_cortex_m4 "portable")
set(CFX_TARGET_PARENT_arm_neon "portable")
set(CFX_TARGET_PARENT_aarch64_neon "arm_neon")

#
# Auto-detection function
#
function(cfx_detect_target OUT_VAR)
    # Auto-detect portabled on compiler and platform
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
        # Check for AVX2/BMI2 support
        include(CheckCCompilerFlag)
        check_c_compiler_flag("-mavx2" HAS_AVX2_FLAG)
        check_c_compiler_flag("-mbmi2" HAS_BMI2_FLAG)

        if(HAS_AVX2_FLAG AND HAS_BMI2_FLAG)
            # Verify compile-time support via compile test
            include(CheckCSourceCompiles)
            set(CMAKE_REQUIRED_FLAGS "-mavx2 -mbmi2")
            check_c_source_compiles("
                #include <immintrin.h>
                int main() {
                    unsigned long long hi;
                    (void)_mulx_u64(1ULL, 2ULL, &hi);
                    return 0;
                }
            " CFX_CAN_USE_BMI2)

            if(CFX_CAN_USE_BMI2)
                set(${OUT_VAR} "x86_64_avx2" PARENT_SCOPE)
                return()
            endif()
        endif()

        set(${OUT_VAR} "x86_64" PARENT_SCOPE)

    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
        set(${OUT_VAR} "aarch64_neon" PARENT_SCOPE)

    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "^arm")
        # Check for NEON
        include(CheckCCompilerFlag)
        check_c_compiler_flag("-mfpu=neon" HAS_NEON_FLAG)
        if(HAS_NEON_FLAG)
            set(${OUT_VAR} "arm_neon" PARENT_SCOPE)
        else()
            set(${OUT_VAR} "portable" PARENT_SCOPE)
        endif()

    else()
        set(${OUT_VAR} "portable" PARENT_SCOPE)
    endif()
endfunction()

#
# Backend source resolution with tree-portabled fallback
#
# Walks up the target hierarchy to find an implementation file.
# Example: cfx_find_backend_source("big" "mul" BIG_MUL_SOURCE)
#
function(cfx_find_backend_source ALGORITHM FUNCTION OUT_SOURCE)
    set(_target ${CFX_TARGET})
    set(_checked "")

    # Walk up the tree until we find an implementation
    while(TRUE)
        set(_path "${CMAKE_SOURCE_DIR}/src/${ALGORITHM}/${_target}/${FUNCTION}.c")
        list(APPEND _checked "${_target}")

        if(EXISTS ${_path})
            set(${OUT_SOURCE} ${_path} PARENT_SCOPE)
            message(STATUS "  cfx/${ALGORITHM}/${FUNCTION}: using ${_target}")
            return()
        endif()

        # Move to parent
        if(DEFINED CFX_TARGET_PARENT_${_target})
            set(_target ${CFX_TARGET_PARENT_${_target}})
        else()
            # Reached root (portable) and still not found
            message(FATAL_ERROR
                "cfx: No implementation of '${FUNCTION}' found for algorithm '${ALGORITHM}'.\n"
                "Checked targets: ${_checked}\n"
                "Expected file: src/${ALGORITHM}/portable/${FUNCTION}.c")
        endif()
    endwhile()
endfunction()

#
# Collect all backend sources for an algorithm
#
# Example: cfx_collect_backend_sources("big" "mul;add;mont" BIG_BACKEND_SOURCES)
#
function(cfx_collect_backend_sources ALGORITHM FUNCTIONS OUT_SOURCES)
    set(_sources "")
    foreach(_func ${FUNCTIONS})
        cfx_find_backend_source(${ALGORITHM} ${_func} _src)
        list(APPEND _sources ${_src})
    endforeach()
    set(${OUT_SOURCES} ${_sources} PARENT_SCOPE)
endfunction()

#
# Apply target-specific compile definitions and flags
#
function(cfx_apply_target target_name)
    # Define the target macro (uppercase with underscores)
    string(TOUPPER ${CFX_TARGET} _TARGET_UPPER)
    target_compile_definitions(${target_name} PRIVATE
        CFX_TARGET_${_TARGET_UPPER}=1
        CFX_TARGET_NAME="${CFX_TARGET}"
    )

    # Apply capability flags and compiler options based on target
    if(CFX_TARGET MATCHES "x86_64_bmi2|x86_64_avx2|x86_64_avx512")
        set(CFX_CAP_BMI2 ON CACHE INTERNAL "" FORCE)
        target_compile_definitions(${target_name} PRIVATE CFX_CAP_BMI2=1)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${target_name} PRIVATE -mbmi2)
        endif()
    endif()

    if(CFX_TARGET MATCHES "x86_64_avx2|x86_64_avx512")
        set(CFX_CAP_AVX2 ON CACHE INTERNAL "" FORCE)
        # PUBLIC so dependents (tests, users) see the same ctx sizes in headers
        target_compile_definitions(${target_name} PUBLIC CFX_CAP_AVX2=1)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${target_name} PRIVATE -mavx2)
        elseif(MSVC)
            target_compile_options(${target_name} PRIVATE /arch:AVX2)
        endif()
    endif()

    if(CFX_TARGET STREQUAL "x86_64_avx512")
        set(CFX_CAP_AVX512 ON CACHE INTERNAL "" FORCE)
        target_compile_definitions(${target_name} PRIVATE CFX_CAP_AVX512=1)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${target_name} PRIVATE -mavx512f -mavx512vl)
        elseif(MSVC)
            target_compile_options(${target_name} PRIVATE /arch:AVX512)
        endif()
    endif()

    if(CFX_TARGET MATCHES "arm_neon|aarch64_neon")
        set(CFX_CAP_NEON ON CACHE INTERNAL "" FORCE)
        target_compile_definitions(${target_name} PRIVATE CFX_CAP_NEON=1)
        if(CFX_TARGET STREQUAL "arm_neon")
            if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
                target_compile_options(${target_name} PRIVATE -mfpu=neon)
            endif()
        endif()
    endif()

    if(CFX_TARGET STREQUAL "arm_cortex_m4")
        set(CFX_CAP_DSP ON CACHE INTERNAL "" FORCE)
        target_compile_definitions(${target_name} PRIVATE CFX_CAP_DSP=1)
        # Force 32-bit limbs for Cortex-M4
        target_compile_definitions(${target_name} PRIVATE CFX_FORCE_LIMB_32=1)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${target_name} PRIVATE
                -mcpu=cortex-m4
                -mthumb
            )
        endif()
    endif()

endfunction()

#
# Perform auto-detection if requested
#
if(CFX_TARGET STREQUAL "auto")
    cfx_detect_target(CFX_TARGET)
    message(STATUS "cfx: Auto-detected target: ${CFX_TARGET}")
    # Update the cache with the detected value
    set(CFX_TARGET ${CFX_TARGET} CACHE STRING "Target architecture for optimized implementations" FORCE)
endif()

# Validate target
list(FIND CFX_TARGETS ${CFX_TARGET} _target_idx)
if(_target_idx EQUAL -1)
    message(FATAL_ERROR "cfx: Unknown target '${CFX_TARGET}'. Valid targets: ${CFX_TARGETS}")
endif()

message(STATUS "cfx: Building for target: ${CFX_TARGET}")

# MEMORY MODE SELECTION
#
# CFX_MEMORY_MODE controls the allocation strategy:
#   - dynamic: uses malloc/realloc/free (default)
#   - static: uses fixed-size buffer pool (for embedded/no-heap environments)
#

set(CFX_MEMORY_MODES dynamic static)

set(CFX_MEMORY_MODE "dynamic" CACHE STRING "Memory allocation strategy")
set_property(CACHE CFX_MEMORY_MODE PROPERTY STRINGS ${CFX_MEMORY_MODES})

# Validate memory mode
list(FIND CFX_MEMORY_MODES ${CFX_MEMORY_MODE} _mem_idx)
if(_mem_idx EQUAL -1)
    message(FATAL_ERROR "cfx: Unknown memory mode '${CFX_MEMORY_MODE}'. Valid modes: ${CFX_MEMORY_MODES}")
endif()

message(STATUS "cfx: Memory mode: ${CFX_MEMORY_MODE}")

#
# Memory source resolution
#
# Finds memory backend source files in src/<algorithm>/mem/<mode>/
# Example: cfx_find_mem_source("big" "init" BIG_MEM_INIT_SOURCE)
#
function(cfx_find_mem_source ALGORITHM FUNCTION OUT_SOURCE)
    set(_path "${CMAKE_SOURCE_DIR}/src/${ALGORITHM}/mem/${CFX_MEMORY_MODE}/${FUNCTION}.c")

    if(EXISTS ${_path})
        set(${OUT_SOURCE} ${_path} PARENT_SCOPE)
        message(STATUS "  cfx/${ALGORITHM}/mem/${FUNCTION}: using ${CFX_MEMORY_MODE}")
    else()
        message(FATAL_ERROR
            "cfx: No memory implementation of '${FUNCTION}' found.\n"
            "Expected file: ${_path}")
    endif()
endfunction()

#
# Apply memory mode compile definitions
#
function(cfx_apply_memory_mode target_name)
    if(CFX_MEMORY_MODE STREQUAL "static")
        # PUBLIC so tests and dependents see this definition
        target_compile_definitions(${target_name} PUBLIC CFX_MEMORY_STATIC=1)
    endif()
endfunction()
