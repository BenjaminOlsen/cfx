if(CFX_BUILD_TESTS)
  add_subdirectory(unit)
endif()

if(CFX_ENABLE_TESTU01)
    add_subdirectory(stats)
endif()
