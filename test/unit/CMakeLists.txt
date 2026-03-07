
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
    test_poly1305
    test_chacha20
    test_base64
    test_aead_chacha20_poly1305
    test_sha256
    test_sha512
    test_sha3
    test_blake2
    test_siphash
    test_vec
    test_compat
    test_ntt
    test_sbig
)

# x86-only tests
if(NOT CMAKE_CROSSCOMPILING OR CFX_ARCH STREQUAL "x86_64")
    list(APPEND CFX_UNIT_TESTS test_chacha_avx2_transpose)
endif()

# Tests that need an OS (pthreads, /dev/urandom, mmap)
if(NOT CFX_TARGET STREQUAL "arm_cortex_m4")
    list(APPEND CFX_UNIT_TESTS
        test_mul_csa_pthreads
        test_mul_rows
        test_rand
    )
endif()

# Tests incompatible with Cortex-M4 bare-metal (RAM limits, no OS, 32-bit perf)
if(CFX_TARGET STREQUAL "arm_cortex_m4")
    list(REMOVE_ITEM CFX_UNIT_TESTS
        test_compat           # requires pthreads
        test_ntt              # large malloc(n * sizeof(uint64_t))
        test_ecm              # many simultaneous big integers
        test_big_exp          # huge intermediate results
        test_big_pollard_rho  # big integer working set
        test_big_arith        # heap exhaustion (large intermediate mul/add)
        test_big_div          # heap exhaustion (division scratch)
        test_big_gcd          # heap exhaustion (uses division internally)
        test_big_mul          # heap exhaustion (large multiplications)
        test_big_str          # heap exhaustion (decimal conversion scratch)
        test_sbig             # heap exhaustion (wraps cfx_big_t)
        test_sha512           # 64-bit ops on 32-bit CPU, too slow under QEMU
        test_algo             # too slow under QEMU emulation
        test_fac              # factorization too slow on 32-bit under QEMU
    )
endif()

# Tests that require 128-bit types (64-bit platforms only)
if(NOT CFX_TARGET MATCHES "arm_cortex_m4|arm_neon" AND NOT CFX_LIMB_BITS STREQUAL "32")
    list(APPEND CFX_UNIT_TESTS
        test_fe25519
        test_x25519
        test_ge25519
        test_ed25519
        test_argon2
        test_poly1271
        test_aead_chacha20_poly1271
    )
endif()

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

    # Link bare-metal support for Cortex-M4
    # --whole-archive ensures syscall stubs are available to libc_nano
    if(CFX_TARGET STREQUAL "arm_cortex_m4")
        target_link_libraries(${test_name} PRIVATE
            -Wl,--whole-archive cfx_m4_support -Wl,--no-whole-archive)
    endif()
endforeach()
