# Bare-metal test support library for Cortex-M4
#
# Provides startup code, semihosting syscalls, and linker script
# so test binaries can run under qemu-system-arm.

add_library(cfx_m4_support STATIC
    startup.c
    syscalls.c
)

target_include_directories(cfx_m4_support PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Custom linker script + bare-metal link options
set(LINKER_SCRIPT ${CMAKE_CURRENT_SOURCE_DIR}/linker.ld)
target_link_options(cfx_m4_support PUBLIC
    -T${LINKER_SCRIPT}
    -nostartfiles
)
