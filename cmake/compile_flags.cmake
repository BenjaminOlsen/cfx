if(NOT DEFINED CFX_ARCH)
    set(CFX_ARCH "native" CACHE STRING "Target architecture")
    set_property(CACHE CFX_ARCH PROPERTY STRINGS native x86_64 armv7m aarch64)
endif()

# Compiler-specific warning flags
if(MSVC)
    set(CFX_CFLAGS_COMMON
        /W4         # High warning level
        /WX         # Warnings as errors
        /wd4100     # unreferenced formal parameter
        /wd4244     # conversion, possible loss of data
        /wd4267     # conversion from size_t
        /wd4706     # assignment within conditional expression
        /wd4996     # deprecated functions
    )
    set(CFX_SILENCED "")
    set(CFX_DEBUG_FLAGS /Od /Zi)
    set(CFX_RELEASE_FLAGS /O2)
    set(CFX_RELWITHDEBINFO_FLAGS /O2 /Zi)
else()
    # GCC/Clang flags
    set(CFX_CFLAGS_COMMON
        -Wall
        -Wextra
        -Werror
        -pedantic
        -Wshadow
        -Wcast-qual
        -Wstrict-prototypes
        -Wmissing-prototypes
        -Wmissing-declarations
        -Wpointer-arith
        -Wold-style-definition
        -Wvla
        -Werror=vla
        -Wredundant-decls
        -Wmissing-field-initializers
        -fno-strict-aliasing
    )
    set(CFX_SILENCED
        # -Wno-gnu-zero-variadic-macro-arguments
    )
    set(CFX_DEBUG_FLAGS -g)
    set(CFX_RELEASE_FLAGS -O3)
    set(CFX_RELWITHDEBINFO_FLAGS -O2 -g)
endif()

set(CFX_ARCH_FLAGS "")
set(CFX_ARCH_DEFINES "")

if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")

    if(CFX_ARCH STREQUAL "native")
        list(APPEND CFX_ARCH_FLAGS -march=native)

    elseif(CFX_ARCH STREQUAL "x86_64")
        list(APPEND CFX_ARCH_FLAGS -march=x86-64-v2)
        list(APPEND CFX_ARCH_DEFINES CFX_ARCH_X86_64=1)

    elseif(CFX_ARCH STREQUAL "armv7m")
        # Cortex-M4
        list(APPEND CFX_ARCH_FLAGS
            -mcpu=cortex-m4
            -mthumb
            -ffreestanding
            -fno-builtin
            -fdata-sections
            -ffunction-sections
        )
        list(APPEND CFX_ARCH_DEFINES CFX_ARCH_ARMV7M=1)

    # todo::
    # elseif(CFX_ARCH STREQUAL "aarch64")
    #     list(APPEND CFX_ARCH_FLAGS -march=armv8-a)
    #     list(APPEND CFX_ARCH_DEFINES CFX_ARCH_AARCH64=1)

    else()
        message(FATAL_ERROR "Unknown CFX_ARCH='${CFX_ARCH}'")
    endif()

elseif(MSVC)
    # MSVC architecture options
    if(CFX_ARCH STREQUAL "x86_64")
        list(APPEND CFX_ARCH_DEFINES CFX_ARCH_X86_64=1)
    endif()
    # native is default for MSVC, no special flags needed

endif()

# Helper to apply flags to a target
function(cfx_apply_compile_flags target)
    target_compile_options(${target} PRIVATE
        ${CFX_CFLAGS_COMMON}
        ${CFX_SILENCED}
        $<$<CONFIG:Debug>:${CFX_DEBUG_FLAGS}>
        $<$<CONFIG:Release>:${CFX_RELEASE_FLAGS}>
        $<$<CONFIG:RelWithDebInfo>:${CFX_RELWITHDEBINFO_FLAGS}>
    )

    target_compile_definitions(${target} PRIVATE
        $<$<CONFIG:Debug>:CFX_DEBUG>
        $<$<CONFIG:Release>:NDEBUG>
        $<$<CONFIG:RelWithDebInfo>:NDEBUG;CFX_DEBUG>
    )

    if(NOT MSVC)
        target_compile_options(${target} PUBLIC ${CFX_ARCH_FLAGS})
    endif()

    if(CFX_ARCH_DEFINES)
        target_compile_definitions(${target} PRIVATE ${CFX_ARCH_DEFINES})
    endif()
endfunction()

string(TOUPPER "${CMAKE_BUILD_TYPE}" CFG)

message(STATUS "cfx: Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "cfx: CFX_ARCH = ${CFX_ARCH}")
message(STATUS "cfx: Common flags: ${CFX_CFLAGS_COMMON}")
message(STATUS "cfx: Arch flags:   ${CFX_ARCH_FLAGS}")
message(STATUS "cfx: Silenced:     ${CFX_SILENCED}")
