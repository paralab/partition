# CMakeLists.txt
cmake_minimum_required(VERSION 3.12)  # Minimum CMake version required

project(Partitioning)  # Project name

set(CMAKE_CXX_COMPILER "mpic++")

# Optionally, set compiler flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -Wall")


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
include(${VTK_USE_FILE})

add_subdirectory(lib)
# Add the executable
# add_executable(main main.cpp)

# add_executable(main-dist main-dist.cpp)
add_executable(main-new main-new.cpp)



# Link against libraries
# target_link_libraries(main PRIVATE graph)
# target_link_libraries(main-dist PRIVATE graph)
# target_link_libraries(main-dist PRIVATE util)
target_link_libraries(main-new PRIVATE util mesh-util vtk-util graph sfc)






