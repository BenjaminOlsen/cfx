
set(CFX_UNIT_TESTS
    test_big
    test_big_mont
    test_fac
    test_algo
    test_pollard_rho
    test_mul_rows
    test_poly1305
    test_chacha20
    test_mul_csa_pthreads
)

foreach(test_name IN LISTS CFX_UNIT_TESTS)
    add_executable(${test_name} ${test_name}.c)
    # set_target_properties(${test_name} PROPERTIES
    #     RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    # )
    target_link_libraries(${test_name} PRIVATE cfx::cfx)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
