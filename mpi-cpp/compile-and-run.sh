set -e

export OMP_NUM_THREADS=8

mkdir -p build

cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache

ninja -C ./build

mpirun -np 4 --bind-to none --oversubscribe  ./build/main-new

export SFC_morton="$PWD/out-sfc.vtk"
export METIS="$PWD/out-sfc.vtk"
export BFS="$PWD/out-sfc.vtk"
export BFS_grow="$PWD/out-sfc.vtk"

/home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ../grow/paraview_script.py
