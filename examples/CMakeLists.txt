
set(EXAMPLE_NAMES
    cfx_keygen
    prime_test
    factorial
    factor
    choose
    encrypt
    modexp
    primes
    fread
    mac
    div
    pow
    primes_near_pow2
    ulam_spiral
    find_prime
    dcmini
    rand
)

add_library(example_utils utils.h utils.c)

foreach(example_name IN LISTS EXAMPLE_NAMES)
    add_executable(${example_name} ${example_name}.c)
    target_link_libraries(${example_name} PRIVATE cfx::cfx example_utils)
endforeach()