
set(EXAMPLE_NAMES
    cfx_choose
    cfx_dc
    cfx_div
    cfx_encrypt
    cfx_factor
    cfx_factorial
    cfx_find_prime
    cfx_fread
    cfx_keygen
    cfx_mac
    cfx_modexp
    cfx_pow
    cfx_prime_test
    cfx_primes
    cfx_primes_near_pow2
    cfx_rand
    cfx_ulam_spiral
)

add_library(example_utils utils.h utils.c)

foreach(example_name IN LISTS EXAMPLE_NAMES)
    add_executable(${example_name} ${example_name}.c)
    target_link_libraries(${example_name} PRIVATE cfx::cfx example_utils)
endforeach()

# Install examples as CLI tools
install(TARGETS ${EXAMPLE_NAMES}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)