# add_library(bfs-com bfs-com.cpp)
# target_include_directories(bfs-com PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(.)
# list(APPEND MAIN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MAIN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(${MAIN_TARGET_NAME} PRIVATE bfs-com.cpp)