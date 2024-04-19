add_library(graph graph.cpp)
target_include_directories(graph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(graph PUBLIC util)
