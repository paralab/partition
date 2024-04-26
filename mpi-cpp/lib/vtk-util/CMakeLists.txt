# add_library(vtk-util vtk-util.cpp)
# target_include_directories(vtk-util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# target_link_libraries(vtk-util PUBLIC util ${VTK_LIBRARIES})
# list(APPEND MAIN_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
# include_directories(${MAIN_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
target_sources(${MAIN_TARGET_NAME} PRIVATE vtk-util.cpp)