add_subdirectory(util)
add_subdirectory(mesh-util)
# 
add_subdirectory(sfc)
add_subdirectory(metis-util)
add_subdirectory(scotch-util)
add_subdirectory(usort SYSTEM)
add_subdirectory(dist-graph)

add_subdirectory(linalg)

if(ENABLE_VTK_FEATURES)
    add_subdirectory(vtk-util)
endif()




