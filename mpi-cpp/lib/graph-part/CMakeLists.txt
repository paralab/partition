# add_library(graph-part graph-part.cpp)
# target_include_directories(graph-part PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# list(APPEND MAIN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MAIN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(${MAIN_TARGET_NAME} PRIVATE graph-part.cpp)

# include_directories(.)