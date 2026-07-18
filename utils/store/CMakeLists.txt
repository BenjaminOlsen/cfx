# Encrypted secret-store component.

add_library(cfx_store STATIC
    cfx_store.c
    cfx_store_crypto.c
    cfx_store_format.c
    cfx_store_grace.c
    cfx_store_helpers.c
)
add_library(cfx::store ALIAS cfx_store)

target_link_libraries(cfx_store PUBLIC cfx::cfx cfx_misc)
target_include_directories(cfx_store PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

if(CFX_BUILD_STORE_TOOL)
    add_executable(cfx_store_tool cfx_store_main.c)
    target_link_libraries(cfx_store_tool PRIVATE cfx_store)
    set_target_properties(cfx_store_tool PROPERTIES OUTPUT_NAME cfx_store)

    install(TARGETS cfx_store_tool
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
