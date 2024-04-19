add_library(mesh-util mesh-util.cpp)
target_include_directories(mesh-util PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(mesh-util PUBLIC util graph ${GMSH_LIBRARIES})
