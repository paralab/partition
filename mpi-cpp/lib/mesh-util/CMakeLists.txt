# add_library(mesh-util mesh-util.cpp)
# target_include_directories(mesh-util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# target_link_libraries(mesh-util PUBLIC util graph usort ${GMSH_LIBRARIES} ${MPI_CXX_LIBRARIES})
# list(APPEND MAIN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MAIN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(${MAIN_TARGET_NAME} PRIVATE mesh-util.cpp)
target_sources(${OCTREE_TARGET_NAME} PRIVATE mesh-util.cpp)


# include_directories(. ../util ../graph ../usort ${MPI_INCLUDE_PATH} ${GMSH_INCLUDE_DIRS})