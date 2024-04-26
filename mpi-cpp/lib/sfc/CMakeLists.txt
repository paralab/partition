# add_library(sfc sfc.cpp)
# target_include_directories(sfc PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# target_link_libraries(sfc PUBLIC util)
# list(APPEND MAIN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MAIN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(${MAIN_TARGET_NAME} PRIVATE sfc.cpp)