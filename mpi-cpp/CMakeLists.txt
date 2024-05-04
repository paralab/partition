# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)  # Minimum CMake version required

project(Partitioning)  # Project name

set(CMAKE_COLOR_MAKEFILE ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20) 


# Find MPI package
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})



# Optionally, set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -pedantic -Wall -fpermissive -fdiagnostics-color=always")


# Manually specify the paths to Gmsh headers and libraries
set(GMSH_INCLUDE_DIRS "/home/budvin/bin/gmsh-4.12.2-Linux64-sdk/include")
set(GMSH_LIBRARIES "/home/budvin/bin/gmsh-4.12.2-Linux64-sdk/lib/libgmsh.so")

# Include Gmsh headers
include_directories(${GMSH_INCLUDE_DIRS})

# Set the path to the VTK build directory
set(VTK_DIR "/home/budvin/bin/VTK-9.3.0/build")

# Find VTK package
find_package(VTK COMPONENTS
    CommonColor
    CommonCore
    CommonSystem
    CommonExecutionModel
    IOLegacy)
# include(${VTK_USE_FILE})

# Specify the directory where METIS is built
set(METIS_INCLUDE_DIR "/home/budvin/bin/METIS-5.2.1/build/include")
set(METIS_LIBRARY "/home/budvin/bin/METIS-5.2.1/build/lib/libmetis.so")

include_directories(${METIS_INCLUDE_DIR})


set(GKLIB_INCLUDE_DIR "/home/budvin/bin/GKlib-master/build/include")
set(GKLIB_LIBRARY "/home/budvin/bin/GKlib-master/build/lib/libGKlib.a")
include_directories(${GKLIB_INCLUDE_DIR})

set(MAIN_TARGET_NAME "main-new")

add_executable(${MAIN_TARGET_NAME} main-new.cpp )
add_subdirectory(lib)

include_directories(lib/bfs-com lib/graph lib/graph-part lib/mesh-util lib/metis-util lib/sfc lib/usort lib/util lib/vtk-util lib/dist-graph)





# get_target_property(sources ${MAIN_TARGET_NAME} SOURCES)
# message(sources="${sources}")



target_link_libraries(main-new PRIVATE ${GMSH_LIBRARIES}  ${VTK_LIBRARIES} ${METIS_LIBRARY} ${GKLIB_LIBRARY} ${MPI_CXX_LIBRARIES})






