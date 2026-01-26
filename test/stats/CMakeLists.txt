# Constant-time test (no external dependencies)
add_executable(test_ct test_ct.c)
target_link_libraries(test_ct PRIVATE cfx::cfx)
target_compile_options(test_ct PRIVATE ${CFX_CFLAGS})
message(STATUS "added stats test 'test_ct'")

# TestU01 tests (optional)
if(TESTU01_ROOT)
    include_directories(${TESTU01_ROOT}/include)
    message(STATUS "added testu01 root ${TESTU01_ROOT}")

    set(CFX_STATS_TESTS
        test_testu01
    )

    foreach(test_name IN LISTS CFX_STATS_TESTS)
        add_executable(${test_name} ${test_name}.c)

        target_link_libraries(${test_name} PRIVATE
            cfx::cfx
            ${TESTU01_ROOT}/lib/libtestu01.a
            ${TESTU01_ROOT}/lib/libprobdist.a
            ${TESTU01_ROOT}/lib/libmylib.a
            m
        )

        message(STATUS "added stats test '${test_name}'")
    endforeach()

    set(TESTS_USING_OPENSSL
        test_testu01
    )

    foreach(test_name IN LISTS TESTS_USING_OPENSSL)
        if (OpenSSL_FOUND)
            message(STATUS "Enabling OpenSSL for ${test_name}")
            target_link_libraries(${test_name} PRIVATE OpenSSL::Crypto)
        else()
            message(STATUS "Not using OpenSSL for ${test_name}")
        endif()
    endforeach()
else()
    message(STATUS "TestU01 not configured (set TESTU01_ROOT to enable)")
endif()
