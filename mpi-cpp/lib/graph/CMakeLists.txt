# add_library(graph graph.cpp)
# target_include_directories(graph PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# target_link_libraries(graph PUBLIC util)
# list(APPEND MAIN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MAIN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(${MAIN_TARGET_NAME} PRIVATE graph.cpp)

# include_directories(. ../util)