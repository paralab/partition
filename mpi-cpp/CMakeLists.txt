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



if(NOT GMSH_SDK_PATH)
    message(FATAL_ERROR "The GMSH_SDK_PATH option is not provided. Please specify the path to the library using -DGMSH_SDK_PATH=/path/to/sdk.")
endif()





# Manually specify the paths to Gmsh headers and libraries
set(GMSH_INCLUDE_DIRS "${GMSH_SDK_PATH}/include")
set(GMSH_LIBRARIES "${GMSH_SDK_PATH}/lib/libgmsh.so")

# Include Gmsh headers
include_directories(${GMSH_INCLUDE_DIRS})


if(NOT VTK_INSTALL_DIR_PATH)
    message(FATAL_ERROR "The VTK_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DVTK_INSTALL_DIR_PATH=/path/to/vtk/build.")
endif()

# Set the path to the VTK build directory
set(VTK_DIR "${VTK_INSTALL_DIR_PATH}")

# Find VTK package
find_package(VTK COMPONENTS
    CommonColor
    CommonCore
    CommonSystem
    CommonExecutionModel
    IOLegacy)

if(NOT METIS_INSTALL_DIR_PATH)
    message(FATAL_ERROR "The METIS_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DMETIS_INSTALL_DIR_PATH=/path/to/metis/build.")
endif()

# Specify the directory where METIS is built
set(METIS_INCLUDE_DIR "${METIS_INSTALL_DIR_PATH}/include")
set(METIS_LIBRARY "${METIS_INSTALL_DIR_PATH}/lib/libmetis.so")

include_directories(${METIS_INCLUDE_DIR})


if(NOT GKLIB_INSTALL_DIR_PATH)
    message(FATAL_ERROR "The GKLIB_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DGKLIB_INSTALL_DIR_PATH=/path/to/gklib/build.")
endif()


set(GKLIB_INCLUDE_DIR "${GKLIB_INSTALL_DIR_PATH}/include")
set(GKLIB_LIBRARY "${GKLIB_INSTALL_DIR_PATH}/lib/libGKlib.a")
include_directories(${GKLIB_INCLUDE_DIR})

if(NOT PARMETIS_INSTALL_DIR_PATH)
    message(FATAL_ERROR "The PARMETIS_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DPARMETIS_INSTALL_DIR_PATH=/path/to/parmetis/build.")
endif()

set(PARMETIS_INCLUDE_DIR "${PARMETIS_INSTALL_DIR_PATH}/include")
set(PARMETIS_LIBRARY "${PARMETIS_INSTALL_DIR_PATH}/lib/libparmetis.so")

include_directories(${PARMETIS_INCLUDE_DIR})

set(MAIN_TARGET_NAME "main-new")

add_executable(${MAIN_TARGET_NAME} main-new.cpp )
add_subdirectory(lib)

include_directories(lib/bfs-com lib/graph lib/graph-part lib/mesh-util lib/metis-util lib/sfc lib/usort lib/util lib/vtk-util lib/dist-graph)



target_link_libraries(main-new PRIVATE ${GMSH_LIBRARIES}  ${VTK_LIBRARIES} ${PARMETIS_LIBRARY} ${METIS_LIBRARY} ${GKLIB_LIBRARY} ${MPI_CXX_LIBRARIES})






