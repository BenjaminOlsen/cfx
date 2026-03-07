# ARM bare-metal toolchain for Cortex-M series
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR armv7-m)
set(CMAKE_C_COMPILER   arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

# Cortex-M4 is Thumb-only — all code must be compiled in Thumb mode
# -ffunction-sections/-fdata-sections allow --gc-sections to discard unused code
set(CMAKE_C_FLAGS_INIT "-mcpu=cortex-m4 -mthumb -ffunction-sections -fdata-sections -DCFX_BAREMETAL=1")
set(CMAKE_CXX_FLAGS_INIT "-mcpu=cortex-m4 -mthumb -ffunction-sections -fdata-sections -DCFX_BAREMETAL=1")

# Use nano.specs for smaller newlib (reduced printf, etc.)
# Semihosting syscall stubs are provided by test/support/cortex-m4/syscalls.c
# --gc-sections strips unreferenced functions/data to fit in constrained flash
set(CMAKE_EXE_LINKER_FLAGS_INIT "--specs=nano.specs -Wl,--gc-sections")

# Avoid finding host libraries
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
