set -e

module load cmake/3.26.0
module load ninja
module load ccache
module load gcc/11.2.0
module load openmpi/4.1.4
module load vtk/9.2.6


export OMP_NUM_THREADS=2

export OMP_PROC_BIND=true 

export OMP_PLACES=threads 

GMSH_SDK_PATH=/uufs/chpc.utah.edu/common/home/u1472227/partition-project/dependencies/gmsh-4.13.0-Linux64-sdk
METIS_INSTALL_DIR_PATH=/uufs/chpc.utah.edu/common/home/u1472227/partition-project/dependencies/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=/uufs/chpc.utah.edu/common/home/u1472227/partition-project/dependencies/GKlib/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=/uufs/chpc.utah.edu/common/home/u1472227/partition-project/dependencies/ParMETIS/build/Linux-x86_64/build

mkdir -p build

cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DGMSH_SDK_PATH=${GMSH_SDK_PATH} \
    -DMETIS_INSTALL_DIR_PATH=${METIS_INSTALL_DIR_PATH} \
    -DGKLIB_INSTALL_DIR_PATH=${GKLIB_INSTALL_DIR_PATH} \
    -DPARMETIS_INSTALL_DIR_PATH=${PARMETIS_INSTALL_DIR_PATH}

ninja -C ./build

export LD_LIBRARY_PATH="${VTK_ROOT}/lib:${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${OPENMPI_ROOT}/lib:${LD_LIBRARY_PATH}"

echo -e "compilation done"

mpirun -np 4 ./build/main-new

# export SFC_morton="$PWD/out-sfc.vtk"
# export METIS="$PWD/out-parmetis.vtk"
# export BFS="$PWD/out-bfs.vtk"
# export BFS_grow="$PWD/out-bfs.vtk"

# /home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ../grow/paraview_script.py
