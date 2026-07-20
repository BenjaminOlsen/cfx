add_executable(test_bge_api test_bge_api.c)
target_link_libraries(test_bge_api PRIVATE cfx::bge)
add_test(NAME test_bge_api COMMAND test_bge_api)
