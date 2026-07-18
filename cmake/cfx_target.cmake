# SPDX-License-Identifier: LGPL-3.0-or-later OR GPL-2.0-or-later
#
# CFX Target Selection System
#
# Usage:
#   cmake -DCFX_TARGET=x86_64_avx2 ..
#   cmake -DCFX_TARGET=arm_cortex_m4 -DCMAKE_TOOLCHAIN_FILE=... ..
#

include_guard(GLOBAL)

# Available targets 
set(CFX_TARGETS
    portable        # Portable C99, no intrinsics (root of inheritance tree)
    x86_64_bmi2     # x86-64 + BMI2 (mulx, adcx, adox)
    x86_64_avx2     # x86-64 + AVX2 + BMI2
    arm_cortex_m4   # ARM Cortex-M4 (32-bit, DSP, no NEON)
)

# Target option with auto-detection default
set(CFX_TARGET "auto" CACHE STRING "Target architecture")
set_property(CACHE CFX_TARGET PROPERTY STRINGS auto ${CFX_TARGETS})

# Target parent relationships (for inheritance / fallback resolution)
# portable has no parent (it's the root)
set(CFX_TARGET_PARENT_x86_64_bmi2 "portable")
set(CFX_TARGET_PARENT_x86_64_avx2 "x86_64_bmi2")
set(CFX_TARGET_PARENT_arm_cortex_m4 "portable")

#
# Auto-detection function
#
function(cfx_detect_target OUT_VAR)
    # Auto-detect based on compiler and platform
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
    endif()

    set(${OUT_VAR} "portable" PARENT_SCOPE)
endfunction()

#
# Backend source resolution with tree-based fallback
#
# Walks up the target hierarchy to find an implementation file.
# Example: cfx_find_backend_source("big" "mul" BIG_MUL_SOURCE)
#
function(cfx_find_backend_source ALGORITHM FUNCTION OUT_SOURCE)
    set(_target ${CFX_TARGET})
    set(_checked "")

    # Walk up the tree until we find an implementation
    while(TRUE)
        set(_path "${PROJECT_SOURCE_DIR}/src/${ALGORITHM}/${_target}/${FUNCTION}.c")
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

    if(CFX_TARGET MATCHES "x86_64_bmi2|x86_64_avx2")
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${target_name} PRIVATE -mbmi2)
        endif()
    endif()

    if(CFX_TARGET STREQUAL "x86_64_avx2")
        # PUBLIC so dependents (tests, users) see the same ctx sizes in headers
        target_compile_definitions(${target_name} PUBLIC CFX_CAP_AVX2=1)
        if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            target_compile_options(${target_name} PRIVATE -mavx2)
        elseif(MSVC)
            target_compile_options(${target_name} PRIVATE /arch:AVX2)
        endif()
    endif()

    if(CFX_TARGET STREQUAL "arm_cortex_m4")
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
    set(CFX_TARGET ${CFX_TARGET} CACHE STRING "Target architecture" FORCE)
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
# Apply memory mode compile definitions
#
function(cfx_apply_memory_mode target_name)
    if(CFX_MEMORY_MODE STREQUAL "static")
        # PUBLIC so tests and dependents see this definition
        target_compile_definitions(${target_name} PUBLIC CFX_MEMORY_STATIC=1)
        if(DEFINED CFX_STATIC_LIMBS)
            target_compile_definitions(${target_name} PUBLIC CFX_STATIC_LIMBS=${CFX_STATIC_LIMBS})
        endif()
        if(DEFINED CFX_STATIC_POOL_SIZE)
            target_compile_definitions(${target_name} PUBLIC CFX_STATIC_POOL_SIZE=${CFX_STATIC_POOL_SIZE})
        endif()
    endif()
endfunction()
