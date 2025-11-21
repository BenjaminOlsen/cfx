
set(EXAMPLE_NAMES
    ulam_spiral
    cfx_keygen
    find_prime
    prime_test
    factorial
    factor
    dcmini
    choose
    encrypt
    primes
    fread
    rand
    mac
    div
    pow
)

add_library(example_utils utils.h utils.c)

foreach(example_name IN LISTS EXAMPLE_NAMES)
    add_executable(${example_name} ${example_name}.c)
    # set_target_properties(${example_name} PROPERTIES
    #     RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
    # )
    target_link_libraries(${example_name} PRIVATE cfx::cfx example_utils)
endforeach()