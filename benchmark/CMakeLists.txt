
enable_language(CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(benchmark REQUIRED)

set(CFX_BENCHMARKS
    bench_big
    bench_mul_vs_csa
    bench_row_mul
    bench_monty
    bench_isqrt
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
