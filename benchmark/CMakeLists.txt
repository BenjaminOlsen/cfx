enable_language(CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(benchmark REQUIRED)

# find_package(OpenSSL QUIET)


set(CFX_BENCHMARKS
    bench_big
    bench_mul_vs_csa
    bench_row_mul
    bench_monty
    bench_isqrt
    bench_rand
    bench_rand_fn_ptr
)

foreach(bench_name IN LISTS CFX_BENCHMARKS)
    add_executable(${bench_name} ${bench_name}.cpp)
    target_compile_options(${bench_name} PRIVATE -O1 -march=native)
    target_link_libraries(${bench_name} PRIVATE
        benchmark::benchmark
        cfx::cfx
        pthread
    )
endforeach()

# Only bench_rand needs OpenSSL + macro
if (OpenSSL_FOUND)
    message(STATUS "OpenSSL found – enabling OpenSSL RNG benchmarks")
    target_link_libraries(bench_rand PRIVATE OpenSSL::Crypto)
    target_compile_definitions(bench_rand PRIVATE CFX_HAVE_OPENSSL=1)
else()
    message(STATUS "OpenSSL not found – OpenSSL RNG benchmarks will be disabled")
endif()
