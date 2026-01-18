# ARM bare-metal toolchain for Cortex-M series
set(CMAKE_SYSTEM_NAME Generic)
set(CMAKE_SYSTEM_PROCESSOR armv7-m)
set(CMAKE_C_COMPILER   arm-none-eabi-gcc)
set(CMAKE_CXX_COMPILER arm-none-eabi-g++)
set(CMAKE_TRY_COMPILE_TARGET_TYPE STATIC_LIBRARY)

# Use nosys.specs for minimal syscall stubs (printf, exit, etc.)
# This provides stub implementations for _exit, _sbrk, _write, etc.
set(CMAKE_EXE_LINKER_FLAGS_INIT "--specs=nosys.specs --specs=nano.specs")

# Avoid finding host libraries
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
