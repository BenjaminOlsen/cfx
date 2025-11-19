
if(NOT TESTU01_ROOT)
  message(FATAL_ERROR "Set TESTU01_ROOT to TestU01 install prefix")
endif()

include_directories(${TESTU01_ROOT}/include)

add_executable(test_chacha20_testu01 test_chacha20_testu01.c)

target_link_libraries(test_chacha20_testu01 PRIVATE
  cfx::cfx
  ${TESTU01_ROOT}/lib/libtestu01.a
  ${TESTU01_ROOT}/lib/libprobdist.a
  ${TESTU01_ROOT}/lib/libmylib.a
  m
)

message(STATUS "added testu01 root ${TESTU01_ROOT}")