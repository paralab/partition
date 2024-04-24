add_library(vtk-util vtk-util.cpp)
target_include_directories(vtk-util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(vtk-util PUBLIC util ${VTK_LIBRARIES})
