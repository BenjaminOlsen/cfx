# Encrypted secret-store command component. Store has no public library API.

add_library(cfx_store_command STATIC
    cfx_store.c
    cfx_store_crypto.c
    cfx_store_format.c
    cfx_store_grace.c
    cfx_store_helpers.c
)

target_link_libraries(cfx_store_command PUBLIC cfx::cfx cfx_misc)
target_include_directories(cfx_store_command PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

if(CFX_BUILD_STORE_TOOL)
    add_executable(cfx_store_exe cfx_store_main.c)
    target_link_libraries(cfx_store_exe PRIVATE cfx_store_command)
    set_target_properties(cfx_store_exe PROPERTIES OUTPUT_NAME cfx_store)

    install(TARGETS cfx_store_exe
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
