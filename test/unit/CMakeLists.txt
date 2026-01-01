
set(CFX_UNIT_TESTS
    test_arith
    test_big_init
    test_big_str
    test_big_hex
    test_big_arith
    test_big_mul
    test_big_div
    test_big_shift
    test_big_exp
    test_big_bitwise
    test_big_bytes
    test_big_gcd
    test_big_pollard_rho
    test_big_mont
    test_ecm
    test_fac
    test_algo
    test_pollard_rho
    test_mul_rows
    test_poly1305
    test_chacha20
    test_chacha_avx2_transpose
    test_mul_csa_pthreads
    test_rand
    test_base64
    test_aead_chacha20_poly1305
    test_sha256
    test_vec
    test_compat
    test_ntt
    test_sbig
    test_fe25519
    test_x25519
    test_ge25519
    test_ed25519
)

foreach(test_name IN LISTS CFX_UNIT_TESTS)
    add_executable(${test_name} ${test_name}.c)
    target_link_libraries(${test_name} PRIVATE cfx::cfx)
    add_test(NAME ${test_name} COMMAND $<TARGET_FILE:${test_name}>)
    message(STATUS "added unit test '${test_name}'")
    target_compile_options(${test_name} PRIVATE
        ${CFX_CFLAGS}
        ${CFX_SILENCED}
    )
    if(NOT MSVC)
        target_compile_options(${test_name} PRIVATE -Wno-missing-prototypes)
    endif()
endforeach()
