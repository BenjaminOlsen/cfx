# CLI regression tests — call cfx_<name>_run() and verify stdout output
set(CFX_CLI_TESTS
    test_cli_randart
    test_cli_hash
    test_cli_hmac
    test_cli_factor
    test_cli_isprime
    test_cli_choose
    test_cli_pow
    test_cli_xgcd
    test_cli_div
    test_cli_dc
    test_cli_store
)

foreach(test_name IN LISTS CFX_CLI_TESTS)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE cfx_util_lib)
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/../unit
        ${CMAKE_SOURCE_DIR}/utils)
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    message(STATUS "added cli test '${test_name}'")
    target_compile_options(${test_name} PRIVATE
        ${CFX_CFLAGS}
        ${CFX_SILENCED}
    )
    if(NOT MSVC)
        target_compile_options(${test_name} PRIVATE -Wno-missing-prototypes)
    endif()
endforeach()
