# List of simple (single-file) utility names
set(UTIL_NAMES
    choose
    dc
    div
    chacha20
    factor
    factorial
    fread
    hash
    hmac
    keygen
    primegen
    mac
    modexp
    pow
    isprime
    primes
    primes_near_pow2
    rand
    randart
    ulam_spiral
    pi
    x25519
    ed25519
    sign
    verify
    aead
    argon2
    xgcd
)

# Helper library with misc utilities (hex parsing, etc.)
add_library(cfx_misc common.h common.c cfx_keyfile.h cfx_keyfile.c cfx_randomart.h cfx_randomart.c)
target_link_libraries(cfx_misc PUBLIC cfx::cfx)

# Multi-file utilities in subdirectories
add_subdirectory(bge)

# Library containing all utility logic (without main functions)
set(UTIL_SOURCES)
foreach(util_name IN LISTS UTIL_NAMES)
    list(APPEND UTIL_SOURCES cfx_${util_name}.c)
    message(STATUS "added utility '${util_name}'")
endforeach()

add_library(cfx_util_lib STATIC ${UTIL_SOURCES})
target_link_libraries(cfx_util_lib PUBLIC cfx::cfx cfx_misc cfx_bge_lib)
target_include_directories(cfx_util_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Build individual standalone binaries (cfx_factor, cfx_dc, etc.)
foreach(util_name IN LISTS UTIL_NAMES)
    add_executable(cfx_${util_name} cfx_${util_name}_main.c)
    target_link_libraries(cfx_${util_name} PRIVATE cfx_util_lib)
endforeach()

# Build unified dispatcher binary
# Note: Named cfx_tool to avoid collision with the cfx library target
add_executable(cfx_tool cfx.c)
target_link_libraries(cfx_tool PRIVATE cfx_util_lib)
set_target_properties(cfx_tool PROPERTIES OUTPUT_NAME "cfx")

# Collect all utility executable names for install
set(UTIL_BINARIES cfx_tool)
foreach(util_name IN LISTS UTIL_NAMES)
    list(APPEND UTIL_BINARIES cfx_${util_name})
endforeach()

# Install all binaries
install(TARGETS ${UTIL_BINARIES}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
