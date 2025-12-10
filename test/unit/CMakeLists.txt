
set(CFX_UNIT_TESTS
    test_arith
    test_big
    test_big_mont
    test_fac
    test_algo
    test_pollard_rho
    test_mul_rows
    test_poly1305
    test_chacha20
    test_chacha_avx2_transpose
    test_mul_csa_pthreads
    test_rand
    test_aead_chacha20_poly1305
    test_sha256
)

foreach(test_name IN LISTS CFX_UNIT_TESTS)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE cfx::cfx)
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    message(STATUS "added unit test '${test_name}'")
    target_compile_options(${test_name} PRIVATE
        ${CFX_CFLAGS}
        ${CFX_SILENCED}
        -Wno-missing-prototypes
    )
endforeach()
