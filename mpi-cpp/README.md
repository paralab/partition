Cmake dependencies

- Gmsh SDK https://gmsh.info/#Download
- VTK https://vtk.org/download/
- METIS https://github.com/KarypisLab/METIS/tags
- ParMETIS (depends on METIS) https://github.com/KarypisLab/ParMETIS
- GKLIB (required for METIS) https://github.com/KarypisLab/GKlib


Note:
Build GKLIB before building METIS. On METIS cmake file, You might have to set `set(GKLIB_PATH "/path-to/GKLIB/build-dir")` if GKLIB is not installed globally.



GKlib

in cmake file, set OFF to ON
option(BUILD_SHARED_LIBS "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)" ON)

make config cc=mpicc prefix=./build openmp=set CFLAGS="-D_POSIX_C_SOURCE=199309L"
make
make install




METIS

make config shared=1 cc=mpicc prefix=./build gklib_path=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build
make install

ParMETIS

make config shared=1 cc=mpicc prefix=./build gklib_path=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build metis_path=$WORK/partition-project/dependencies/METIS-5.2.1/build/build CFLAGS="-D_POSIX_C_SOURCE=199309L"
make install




gmsh

cd build
cmake -DENABLE_BUILD_DYNAMIC=1 -DCMAKE_INSTALL_PREFIX=/work2/10000/budvin/frontera/partition-project/dependencies/gmsh-4.13.0-source/build/install  ..
make -j 4
make install