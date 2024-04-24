# CMakeLists.txt
cmake_minimum_required(VERSION 3.12)  # Minimum CMake version required

project(Partitioning)  # Project name

set(CMAKE_COLOR_MAKEFILE ON)


# Find MPI package
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})



# Optionally, set compiler flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -pedantic -Wall -fdiagnostics-color=always")


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


add_subdirectory(lib)
# Add the executable
# add_executable(main main.cpp)

# add_executable(main-dist main-dist.cpp)
add_executable(main-new main-new.cpp)



# Link against libraries
# target_link_libraries(main PRIVATE graph)
# target_link_libraries(main-dist PRIVATE graph)
# target_link_libraries(main-dist PRIVATE util)
# target_link_libraries(main-new PRIVATE util mesh-util vtk-util graph sfc metis-util ${METIS_LIBRARY} ${GKLIB_LIBRARY} ${MPI_CXX_LIBRARIES})
target_link_libraries(main-new PRIVATE util mesh-util vtk-util graph sfc metis-util ${MPI_CXX_LIBRARIES})






