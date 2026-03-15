# BGE / Store — encrypted secrets store and file encryption
#
# Multi-file utility: crypto core, store manipulation, file encryption,
# auth grace period, and two entry points (cfx_bge, cfx_store).

add_library(cfx_bge_lib STATIC
    cfx_bge.c
    cfx_store.c
    cfx_bge_crypto.c
    cfx_bge_store.c
    cfx_bge_file.c
    cfx_bge_common.c
    cfx_bge_grace.c
)

target_link_libraries(cfx_bge_lib PUBLIC cfx::cfx cfx_misc)
target_include_directories(cfx_bge_lib PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}       # cfx_bge_internal.h
    ${CMAKE_CURRENT_SOURCE_DIR}/..    # misc.h, cfx_cmd.h, cfx_keyfile.h
)

# Standalone binaries
add_executable(cfx_bge cfx_bge_main.c)
target_link_libraries(cfx_bge PRIVATE cfx_bge_lib)

add_executable(cfx_store cfx_store_main.c)
target_link_libraries(cfx_store PRIVATE cfx_bge_lib)

install(TARGETS cfx_bge cfx_store
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
