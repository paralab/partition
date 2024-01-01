# CMakeLists.txt
cmake_minimum_required(VERSION 3.12)  # Minimum CMake version required

project(Partitioning)  # Project name

set(CMAKE_CXX_COMPILER "mpic++")

# Optionally, set compiler flags
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -Wall")

add_subdirectory(lib)
# Add the executable
add_executable(main main.cpp)

add_executable(main-dist main-dist.cpp)


# Link against libraries
target_link_libraries(main PRIVATE graph)
target_link_libraries(main-dist PRIVATE graph)
target_link_libraries(main-dist PRIVATE util)


