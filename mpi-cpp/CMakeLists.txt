# CMakeLists.txt
cmake_minimum_required(VERSION 3.20)  # Minimum CMake version required

project(Partitioning CXX)  # Project name

set(CMAKE_COLOR_MAKEFILE ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20) 



# indexing data type for graph encoding
set(DEFAULT_GRAPH_INDEXING_TYPE "32")

if(GRAPH_INDEXING_TYPE)
    if(NOT GRAPH_INDEXING_TYPE MATCHES "^(16|32|64)$")
        message(FATAL_ERROR "Invalid GRAPH_INDEXING_TYPE specified.
            Given: ${GRAPH_INDEXING_TYPE}
            Expected one of: 16 32 64")
    endif()
else()
    set(GRAPH_INDEXING_TYPE "${DEFAULT_GRAPH_INDEXING_TYPE}")
endif()
message("GRAPH_INDEXING_TYPE:${GRAPH_INDEXING_TYPE}")
add_compile_definitions(GRAPH_INDEXING_TYPE=${GRAPH_INDEXING_TYPE})



# data type for BFS distance
set(DEFAULT_BFS_DISTANCE_TYPE "32")

if(BFS_DISTANCE_TYPE)
    if(NOT BFS_DISTANCE_TYPE MATCHES "^(16|32|64)$")
        message(FATAL_ERROR "Invalid BFS_DISTANCE_TYPE specified.
            Given: ${BFS_DISTANCE_TYPE}
            Expected one of: 16 32 64")
    endif()
else()
    set(BFS_DISTANCE_TYPE "${DEFAULT_BFS_DISTANCE_TYPE}")
endif()
message("BFS_DISTANCE_TYPE:${BFS_DISTANCE_TYPE}")
add_compile_definitions(BFS_DISTANCE_TYPE=${BFS_DISTANCE_TYPE})


# data type for BFS label
set(DEFAULT_BFS_LABEL_TYPE "32")

if(BFS_LABEL_TYPE)
    if(NOT BFS_LABEL_TYPE MATCHES "^(16|32|64)$")
        message(FATAL_ERROR "Invalid BFS_LABEL_TYPE specified.
            Given: ${BFS_LABEL_TYPE}
            Expected one of: 16 32 64")
    endif()
else()
    set(BFS_LABEL_TYPE "${DEFAULT_BFS_LABEL_TYPE}")
endif()
message("BFS_LABEL_TYPE:${BFS_LABEL_TYPE}")
add_compile_definitions(BFS_LABEL_TYPE=${BFS_LABEL_TYPE})



# Find MPI package
include(FindMPI)

include_directories(${MPI_INCLUDE_PATH})

# Optionally, set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fopenmp -pedantic -Wall -fpermissive -fdiagnostics-color=always")



if(NOT GMSH_SDK_PATH)
    message(FATAL_ERROR "The GMSH_SDK_PATH option is not provided. Please specify the path to the library using -DGMSH_SDK_PATH=/path/to/sdk.")
endif()




# Include Gmsh headers
set(GMSH_INCLUDE_DIRS "${GMSH_SDK_PATH}/include")
include_directories(${GMSH_INCLUDE_DIRS})

# Add Gmsh library
if(EXISTS "${GMSH_SDK_PATH}/lib64/libgmsh.so")
    set(GMSH_LIBRARIES "${GMSH_SDK_PATH}/lib64/libgmsh.so")
elseif(EXISTS "${GMSH_SDK_PATH}/lib/libgmsh.so")
    set(GMSH_LIBRARIES "${GMSH_SDK_PATH}/lib/libgmsh.so")
else()
    message(FATAL_ERROR "Could not find Gmsh library .so ${GMSH_SDK_PATH}/lib64/libgmsh.so or ${GMSH_SDK_PATH}/lib/libgmsh.so")
endif()

option(ENABLE_VTK_FEATURES "Link VTK to output partitioning visualization." OFF)

if(ENABLE_VTK_FEATURES)
    message("VTK features enabled")
    add_compile_definitions(ENABLE_VTK_FEATURES=1)
    if(DEFINED ENV{VTK_ROOT} AND NOT "$ENV{VTK_ROOT}" STREQUAL "")
        set(VTK_DIR ${VTK_ROOT})
    else()
        if(NOT VTK_INSTALL_DIR_PATH)
            message(FATAL_ERROR "The VTK_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DVTK_INSTALL_DIR_PATH=/path/to/vtk/build.")
        endif()
        set(VTK_DIR "${VTK_INSTALL_DIR_PATH}")

    endif()


    # Find VTK package
    find_package(VTK COMPONENTS
        CommonColor
        CommonCore
        CommonSystem
        CommonExecutionModel
        IOLegacy)
endif()






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
set(GKLIB_LIBRARY "${GKLIB_INSTALL_DIR_PATH}/lib/libGKlib.so")
include_directories(${GKLIB_INCLUDE_DIR})

if(NOT PARMETIS_INSTALL_DIR_PATH)
    message(FATAL_ERROR "The PARMETIS_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DPARMETIS_INSTALL_DIR_PATH=/path/to/parmetis/build.")
endif()

set(PARMETIS_INCLUDE_DIR "${PARMETIS_INSTALL_DIR_PATH}/include")
set(PARMETIS_LIBRARY "${PARMETIS_INSTALL_DIR_PATH}/lib/libparmetis.so")

include_directories(${PARMETIS_INCLUDE_DIR})


if(NOT PETSC_INSTALL_DIR_PATH)
    message(FATAL_ERROR "The PETSC_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DPETSC_INSTALL_DIR_PATH=/path/to/petsc/build ")
endif()

set(PETSC_INCLUDE_DIR "${PETSC_INSTALL_DIR_PATH}/include")
set(PETSC_LIBRARY "${PETSC_INSTALL_DIR_PATH}/lib/libpetsc.so")

include_directories(${PETSC_INCLUDE_DIR})


if(ENABLE_PETSC_DRAW_FEATURES)
    message("PETSC draw features enabled")
    add_compile_definitions(ENABLE_PETSC_DRAW_FEATURES=1)
endif()




if(NOT SCOTCH_INSTALL_DIR_PATH)
    message(FATAL_ERROR "The SCOTCH_INSTALL_DIR_PATH option is not provided. Please specify the path to the library using -DSCOTCH_INSTALL_DIR_PATH=/path/to/scotch/build ")
endif()

set(SCOTCH_INCLUDE_DIR "${SCOTCH_INSTALL_DIR_PATH}/include")

set(SCOTCH_LIBRARY "${SCOTCH_INSTALL_DIR_PATH}/lib/libscotch.so")
set(SCOTCHERR_LIBRARY "${SCOTCH_INSTALL_DIR_PATH}/lib/libscotcherr.so")
set(PTSCOTCH_LIBRARY "${SCOTCH_INSTALL_DIR_PATH}/lib/libptscotch.so")
set(PTSCOTCHERR_LIBRARY "${SCOTCH_INSTALL_DIR_PATH}/lib/libptscotcherr.so")


include_directories(${SCOTCH_INCLUDE_DIR})



# a python helper script is used for exporting partition metrics to json
# hence we need python dev
find_package(Python3 COMPONENTS Interpreter Development)

message("Python_FOUND:${Python3_FOUND}")
message("Python_VERSION:${Python3_VERSION}")

include_directories(${Python3_INCLUDE_DIRS})

set(MAIN_TARGET_NAME "main-new")

add_executable(${MAIN_TARGET_NAME} main-new.cpp )
add_subdirectory(lib)

include_directories(lib/util lib/mesh-util lib/metis-util lib/sfc lib/usort lib/dist-graph lib/linalg)

target_link_libraries(main-new PRIVATE ${GMSH_LIBRARIES} ${PARMETIS_LIBRARY} ${METIS_LIBRARY} ${GKLIB_LIBRARY} ${MPI_CXX_LIBRARIES} ${Python3_LIBRARIES} ${PETSC_LIBRARY} ${SCOTCH_LIBRARY} ${SCOTCHERR_LIBRARY} ${PTSCOTCH_LIBRARY} ${PTSCOTCHERR_LIBRARY})

if(ENABLE_VTK_FEATURES)
    include_directories(lib/vtk-util)
    target_link_libraries(main-new PRIVATE ${VTK_LIBRARIES})
endif()

## helper to perform initial gmsh partitioning
add_executable(gmsh_partition gmsh_partition.cpp)
target_link_libraries(gmsh_partition PRIVATE ${GMSH_LIBRARIES})








