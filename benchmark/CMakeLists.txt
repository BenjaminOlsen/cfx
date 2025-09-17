cmake_minimum_required(VERSION 3.10)
project(gcd_bench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(benchmark REQUIRED)
find_package(cfx REQUIRED)

add_executable(bench_big bench_big.cpp)
target_compile_options(bench_big PRIVATE -O1 -march=native)
target_link_libraries(bench_big PRIVATE benchmark::benchmark cfx::cfx)

add_executable(bench_mul_vs_csa bench_mul_vs_csa.cpp)
target_compile_options(bench_mul_vs_csa PRIVATE -O1 -march=native)
target_link_libraries(bench_mul_vs_csa PRIVATE benchmark::benchmark cfx::cfx)
