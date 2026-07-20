
add_library(cfx_bge STATIC
    cfx_bge_api.c
)
add_library(cfx::bge ALIAS cfx_bge)

target_link_libraries(cfx_bge PUBLIC cfx::cfx)
target_include_directories(cfx_bge
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/..
)

add_library(cfx_bge_command STATIC
    cfx_bge.c
)

target_link_libraries(cfx_bge_command PUBLIC cfx_bge PRIVATE cfx_misc)
target_include_directories(cfx_bge_command PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

install(TARGETS cfx_bge
    EXPORT cfxTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/bge.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

if(CFX_BUILD_BGE_TOOL)
    add_executable(cfx_bge_exe cfx_bge_main.c)
    target_link_libraries(cfx_bge_exe PRIVATE cfx_bge_command)
    set_target_properties(cfx_bge_exe PROPERTIES OUTPUT_NAME cfx_bge)

    install(TARGETS cfx_bge_exe
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
