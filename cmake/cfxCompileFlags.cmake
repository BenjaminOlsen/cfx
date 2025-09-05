set(WARNINGS "-Wall -Wextra -Wpedantic")
set(SILENCED "-Wno-gnu-zero-variadic-macro-arguments")

string(APPEND CMAKE_C_FLAGS_DEBUG          " -g -DCFX_DEBUG ${WARNINGS} ${SILENCED}")
string(APPEND CMAKE_C_FLAGS_RELEASE        " -O3 -DNDEBUG ${WARNINGS} ${SILENCED}")
string(APPEND CMAKE_C_FLAGS_RELWITHDEBINFO " -O2 -g -DNDEBUG -DCFX_DEBUG ${WARNINGS} ${SILENCED}")

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")