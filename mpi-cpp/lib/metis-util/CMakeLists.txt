# add_library(metis-util metis-util.cpp)
# target_include_directories(metis-util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


# target_link_libraries(metis-util PUBLIC util ${METIS_LIBRARY} ${GKLIB_LIBRARY})
# list(APPEND MAIN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MAIN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(${MAIN_TARGET_NAME} PRIVATE metis-util.cpp)
target_sources(${OCTREE_TARGET_NAME} PRIVATE metis-util.cpp)


# include_directories(. ../util ${METIS_INCLUDE_DIR})