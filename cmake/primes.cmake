# primes.cmake - Generates primes table at build time

function(cfx_generate_primes PRIMES_TXT OUT_C PRIMES_COUNT)
    # Validate count at configure time (fast check)
    if(PRIMES_COUNT LESS 54)
        message(FATAL_ERROR
            "CFX_PRIMES_COUNT (${PRIMES_COUNT}) < 54 - "
            "this breaks cfx_factor_u64 (requires first 54 primes)")
    endif()

    message(STATUS "CFX_PRIMES_COUNT ${PRIMES_COUNT}")

    # Generate at build time, only when primes.txt changes
    add_custom_command(
        OUTPUT "${OUT_C}"
        COMMAND ${CMAKE_COMMAND}
            -DPRIMES_TXT=${PRIMES_TXT}
            -DOUT_C=${OUT_C}
            -DPRIMES_COUNT=${PRIMES_COUNT}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/generate_primes.cmake
        DEPENDS "${PRIMES_TXT}"
        COMMENT "Generating primes table (${PRIMES_COUNT} primes)"
        VERBATIM
    )

    # Export path back to caller
    set(CFX_PRIMES_GENERATED_SRC "${OUT_C}" PARENT_SCOPE)
endfunction()
