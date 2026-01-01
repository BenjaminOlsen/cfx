# generate_primes.cmake - Run with: cmake -P generate_primes.cmake
# Required variables: PRIMES_TXT, OUT_C, PRIMES_COUNT

if(NOT DEFINED PRIMES_TXT OR NOT DEFINED OUT_C OR NOT DEFINED PRIMES_COUNT)
    message(FATAL_ERROR "Usage: cmake -DPRIMES_TXT=<file> -DOUT_C=<file> -DPRIMES_COUNT=<n> -P generate_primes.cmake")
endif()

file(READ "${PRIMES_TXT}" CFX_PRIMES_TXT)

string(REGEX REPLACE "[, \t\r\n]+" ";" CFX_PRIMES_LIST "${CFX_PRIMES_TXT}")
list(REMOVE_ITEM CFX_PRIMES_LIST "")

list(LENGTH CFX_PRIMES_LIST CFX_PRIMES_TOTAL)
if(PRIMES_COUNT GREATER CFX_PRIMES_TOTAL)
    message(FATAL_ERROR "PRIMES_COUNT (${PRIMES_COUNT}) > total primes (${CFX_PRIMES_TOTAL})")
endif()

if(PRIMES_COUNT LESS 54)
    message(FATAL_ERROR "PRIMES_COUNT (${PRIMES_COUNT}) < 54 - breaks cfx_factor_u64")
endif()

list(SUBLIST CFX_PRIMES_LIST 0 ${PRIMES_COUNT} CFX_PRIMES_USED)
string(JOIN ", " CFX_PRIMES_INIT ${CFX_PRIMES_USED})

file(WRITE "${OUT_C}" "/* Auto-generated. Do not edit. */\n")
file(APPEND "${OUT_C}" "#include \"cfx/primes.h\"\n\n")
file(APPEND "${OUT_C}" "const uint32_t cfx_primes[] = {\n    ${CFX_PRIMES_INIT}\n};\n\n")
file(APPEND "${OUT_C}" "const size_t cfx_primes_len = ${PRIMES_COUNT}u;\n")

message(STATUS "Generated ${OUT_C} with ${PRIMES_COUNT} primes")
