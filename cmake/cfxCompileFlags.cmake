if (CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    set(CFX_CFLAGS
        -Wall
        -Wextra
        -Werror
        -pedantic
        -Wshadow
        -Wcast-qual
        -Wstrict-prototypes
        -Wmissing-prototypes
        -Wmissing-declarations
        -Wpointer-arith
        -Wold-style-definition
        -Wvla
        -Werror=vla
        -Wredundant-decls
        -Wmissing-field-initializers
    )

    # set(SAN_FLAGS -fsanitize=address,undefined -fno-strict-aliasing -fno-omit-frame-pointer -O1)

    set(CFX_SILENCED
        # -Wno-gnu-zero-variadic-macro-arguments
    )

    target_compile_options(cfx PRIVATE
        ${CFX_CFLAGS}
        ${CFX_SILENCED}

        "$<$<CONFIG:Debug>:-g;-DCFX_DEBUG>"
        "$<$<CONFIG:Release>:-O3;-DNDEBUG;-mbmi2;-madx>"
        "$<$<CONFIG:RelWithDebInfo>:-O2;-g;-DNDEBUG;-DCFX_DEBUG>"
    )
endif()

string(TOUPPER  ${CMAKE_BUILD_TYPE} CFG)

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if (CFG STREQUAL "DEBUG")
    message(STATUS "Flags: ${CMAKE_C_FLAGS_DEBUG} ${CFX_CFLAGS} ${CFX_SILENCED}")
elseif (CFG STREQUAL "RELEASE")
    message(STATUS "Flags: ${CMAKE_C_FLAGS_RELEASE}${CFX_CFLAGS} ${CFX_SILENCED}")
elseif (CFG STREQUAL "RELWITHDEBINFO")
    message(STATUS "Flags: ${CMAKE_C_FLAGS_RELWITHDEBINFO}${CFX_CFLAGS} ${CFX_SILENCED}")
else()
    message(STATUS "Unknown build type" )
endif()