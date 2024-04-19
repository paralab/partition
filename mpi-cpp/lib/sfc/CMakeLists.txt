add_library(sfc sfc.cpp)
target_include_directories(sfc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(sfc PUBLIC util)
