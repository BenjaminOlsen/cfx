# BGE file encryption component.

add_library(cfx_bge STATIC
    cfx_bge.c
    cfx_bge_file.c
)
add_library(cfx::bge ALIAS cfx_bge)

target_link_libraries(cfx_bge PUBLIC cfx::cfx cfx_misc)
target_include_directories(cfx_bge PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

if(CFX_BUILD_BGE_TOOL)
    add_executable(cfx_bge_tool cfx_bge_main.c)
    target_link_libraries(cfx_bge_tool PRIVATE cfx_bge)
    set_target_properties(cfx_bge_tool PROPERTIES OUTPUT_NAME cfx_bge)

    install(TARGETS cfx_bge_tool
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
