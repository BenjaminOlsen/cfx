
if(NOT TESTU01_ROOT)
  message(FATAL_ERROR "Set TESTU01_ROOT to TestU01 install prefix")
endif()

include_directories(${TESTU01_ROOT}/include)
message(STATUS "added testu01 root ${TESTU01_ROOT}")

set(CFX_STATS_TESTS
    test_testu01
)

foreach(test_name IN LISTS CFX_STATS_TESTS)
    add_executable(${test_name} ${test_name}.c)

    # if(CFX_ENABLE_TESTU01)
    #     target_compile_definitions(${test_name} PRIVATE CFX_ENABLE_TESTU01=1)
    # endif()
    target_link_libraries(test_testu01 PRIVATE
        cfx::cfx
        ${TESTU01_ROOT}/lib/libtestu01.a
        ${TESTU01_ROOT}/lib/libprobdist.a
        ${TESTU01_ROOT}/lib/libmylib.a
        m
    )

    message(STATUS "added stats test '${test_name}'")
endforeach()
