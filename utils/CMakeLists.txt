set(UTIL_NAMES
    cfx_choose
    cfx_dc
    cfx_div
    cfx_chacha20
    cfx_factor
    cfx_factorial
    cfx_fread
    cfx_hash
    cfx_hmac
    cfx_keygen
    cfx_primegen
    cfx_mac
    cfx_modexp
    cfx_pow
    cfx_isprime
    cfx_primes
    cfx_primes_near_pow2
    cfx_rand
    cfx_randart
    cfx_ulam_spiral
    cfx_pi
    cfx_x25519
    cfx_ed25519
    cfx_sign
    cfx_verify
    cfx_aead
    cfx_argon2
    cfx_xgcd
)

# Helper library with misc utilities (hex parsing, etc.)
add_library(cfx_misc cfx_utils_common.h cfx_utils_common.c cfx_keyfile.h cfx_keyfile.c cfx_randomart.h cfx_randomart.c)
target_link_libraries(cfx_misc PUBLIC cfx::cfx)

if(CFX_BUILD_BGE)
    add_subdirectory(bge)
endif()

if(CFX_BUILD_STORE)
    add_subdirectory(store)
endif()

# Library containing all utility logic (without main functions)
if(CFX_BUILD_UTILS)
set(UTIL_SOURCES)
foreach(util_name IN LISTS UTIL_NAMES)
    list(APPEND UTIL_SOURCES ${util_name}.c)
    message(STATUS "added utility '${util_name}'")
endforeach()

add_library(cfx_util_lib STATIC ${UTIL_SOURCES})
target_link_libraries(cfx_util_lib PUBLIC cfx::cfx cfx_misc)
if(CFX_BUILD_BGE)
    target_link_libraries(cfx_util_lib PUBLIC cfx_bge_command)
    target_compile_definitions(cfx_util_lib PUBLIC CFX_HAVE_BGE)
endif()
if(CFX_BUILD_STORE)
    target_link_libraries(cfx_util_lib PUBLIC cfx_store_command)
    target_compile_definitions(cfx_util_lib PUBLIC CFX_HAVE_STORE)
endif()
target_include_directories(cfx_util_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Build individual standalone binaries (cfx_factor, cfx_dc, etc.)
foreach(util_name IN LISTS UTIL_NAMES)
    add_executable(${util_name} ${util_name}_main.c)
    target_link_libraries(${util_name} PRIVATE cfx_util_lib)
endforeach()

# Build dispatcher binary
# Note: Named cfx_tool to avoid collision with the cfx library target
add_executable(cfx_tool cfx.c)
target_link_libraries(cfx_tool PRIVATE cfx_util_lib)
set_target_properties(cfx_tool PROPERTIES OUTPUT_NAME "cfx")

# Collect all utility executable names for install
set(UTIL_BINARIES cfx_tool)
foreach(util_name IN LISTS UTIL_NAMES)
    list(APPEND UTIL_BINARIES ${util_name})
endforeach()

# Install all binaries
install(TARGETS ${UTIL_BINARIES}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
endif()
