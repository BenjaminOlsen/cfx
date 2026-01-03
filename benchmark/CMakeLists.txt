enable_language(CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(benchmark REQUIRED)

set(CFX_BENCHMARKS
    bench_big
    bench_mul_vs_csa
    bench_row_mul
    bench_ntt
    bench_monty
    bench_isqrt
    bench_rand
    bench_rand_fn_ptr
    bench_chacha20
    bench_aead
    bench_hash
    bench_poly_mac
)

find_package(Threads REQUIRED)

foreach(bench_name IN LISTS CFX_BENCHMARKS)
    add_executable(${bench_name} ${bench_name}.cpp)
    if(MSVC)
        # Use generator expressions to avoid /O2 vs /RTC1 conflict in Debug builds
        target_compile_options(${bench_name} PRIVATE
            $<$<CONFIG:Release>:/O2>
            $<$<CONFIG:RelWithDebInfo>:/O2>
        )
    else()
        target_compile_options(${bench_name} PRIVATE -O1 -march=native)
    endif()
    target_link_libraries(${bench_name} PRIVATE
        benchmark::benchmark
        cfx::cfx
        Threads::Threads
    )
endforeach()

set(TESTS_USING_OPENSSL
    bench_rand
)

foreach(bench_name IN LISTS TESTS_USING_OPENSSL)
    if (OpenSSL_FOUND)
        message(STATUS "Enabling OpenSSL for ${bench_name}")
        target_link_libraries(${bench_name} PRIVATE OpenSSL::Crypto)
    else()
        message(STATUS "Disabling OpenSSL for ${bench_name}")
    endif()
endforeach()
