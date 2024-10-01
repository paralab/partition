Cmake dependencies

- Gmsh SDK https://gmsh.info/#Download
- VTK https://vtk.org/download/
- METIS https://github.com/KarypisLab/METIS/tags
- ParMETIS (depends on METIS) https://github.com/KarypisLab/ParMETIS
- GKLIB (required for METIS) https://github.com/KarypisLab/GKlib
- PtScotch https://gitlab.inria.fr/scotch/scotch
    - bison and flex should be installed in the system to build PtScotch


Note:
Build GKLIB before building METIS. On METIS cmake file, You might have to set `set(GKLIB_PATH "/path-to/GKLIB/build-dir")` if GKLIB is not installed globally.



GKlib

in cmake file, set OFF to ON
option(BUILD_SHARED_LIBS "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)" ON)

make config cc=mpiicc prefix=./build openmp=set CFLAGS="-D_POSIX_C_SOURCE=199309L"
make
make install




METIS

make config shared=1 cc=mpiicc prefix=./build i64=1 gklib_path=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build
make install

ParMETIS

make config shared=1 cc=mpiicc prefix=./build gklib_path=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build metis_path=$WORK/partition-project/dependencies/METIS-5.2.1/build/build CFLAGS="-D_POSIX_C_SOURCE=199309L"
make install




gmsh
    add -DENABLE_FLTK=0 on server environments without displays
    opt(UNTANGLE "Enable 2D and 3D UNTANGLER" OFF) on CMakeLists.txt
cd build
cmake -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DENABLE_FLTK=0 -DENABLE_BUILD_DYNAMIC=1 -DCMAKE_INSTALL_PREFIX=/work2/10000/budvin/frontera/partition-project/dependencies/gmsh-repo/build/install  ..
make -j 4
make install


PtScotch
Note: if there is an error related to PRIu64 set(CMAKE_C_STANDARD 99) in the CMake file in scotch root dir, and add -DCMAKE_C_FLAGS="-D__STDC_FORMAT_MACROS" to cmake
in cmake file add
    add_definitions(-DSCOTCH_MPI_ASYNC_COLL)
    turn OFF threads and MPI threads
mkdir build && cd build 

cmake -DCMAKE_C_FLAGS="-D__STDC_FORMAT_MACROS" -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_INSTALL_PREFIX=$WORK/partition-project/dependencies/scotch/build/install/ -DBUILD_SHARED_LIBS=ON -DMPI_HOME=$TACC_IMPI_DIR/ -DSCOTCH_MPI_ASYNC_COLL=ON ..

make -j5 
make install