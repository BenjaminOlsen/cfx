set(CFX_UTIL_TESTS
    test_randart
)

foreach(test_name IN LISTS CFX_UTIL_TESTS)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE cfx_misc cfx::cfx)
    target_include_directories(${test_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../unit)
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    message(STATUS "added util test '${test_name}'")
    target_compile_options(${test_name} PRIVATE
        ${CFX_CFLAGS}
        ${CFX_SILENCED}
    )
    if(NOT MSVC)
        target_compile_options(${test_name} PRIVATE -Wno-missing-prototypes)
    endif()

    # Link bare-metal support for Cortex-M4
    if(CFX_TARGET STREQUAL "arm_cortex_m4")
        target_link_libraries(${test_name} PRIVATE
            -Wl,--whole-archive cfx_m4_support -Wl,--no-whole-archive)
    endif()
endforeach()
