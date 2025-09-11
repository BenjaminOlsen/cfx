cmake_minimum_required(VERSION 3.10)
project(gcd_bench LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

find_package(benchmark REQUIRED)
find_package(cfx REQUIRED)

add_executable(big_bench big_bench.cpp)
target_compile_options(big_bench PRIVATE -O1 -march=native)
target_link_libraries(big_bench PRIVATE benchmark::benchmark cfx::cfx)
