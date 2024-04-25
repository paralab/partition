add_library(usort binUtils.cpp parUtils.cpp sort_profiler.cpp)
target_include_directories(usort PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(usort PUBLIC mesh-util)
