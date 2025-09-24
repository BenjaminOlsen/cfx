cmake_minimum_required(VERSION 3.10)
project(gcd_bench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(benchmark REQUIRED)
find_package(cfx REQUIRED)

set(CFX_BENCHMARKS
    bench_big
    bench_mul_vs_csa
    bench_row_mul
    bench_monty
)

foreach(bench_name IN LISTS CFX_BENCHMARKS)
    add_executable(${bench_name} ${bench_name}.cpp)
    target_compile_options(${bench_name} PRIVATE -O1 -march=native)
    target_link_libraries(${bench_name} PRIVATE benchmark::benchmark cfx::cfx pthread)
endforeach()
