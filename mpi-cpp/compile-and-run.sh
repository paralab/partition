set -e

export OMP_NUM_THREADS=1

# export OMP_PROC_BIND=true 

# export OMP_PLACES=threads 

mkdir -p build

cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DGMSH_SDK_PATH=/home/budvin/bin/gmsh-4.12.2-Linux64-sdk \
    -DVTK_INSTALL_DIR_PATH=/home/budvin/bin/VTK-9.3.0/build -DMETIS_INSTALL_DIR_PATH=/home/budvin/bin/METIS-5.2.1/build \
    -DGKLIB_INSTALL_DIR_PATH=/home/budvin/bin/GKlib-master/build \
    -DPARMETIS_INSTALL_DIR_PATH=/home/budvin/bin/ParMETIS/build/Linux-x86_64/build

ninja -C ./build

mpirun -np 25 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh  # octopus
# mpirun -np 4 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh
# mpirun -np 25 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh  # smallest tet
# mpirun -np 8 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/10k_tet/196209_sf_hexa.mesh_73346_289961.obj.mesh #large tet

export SFC_morton="$PWD/out-sfc.vtk"
export METIS="$PWD/out-parmetis.vtk"
export BFS="$PWD/out-bfs.vtk"
export BFS_grow="$PWD/out-grow.vtk"

/home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ../grow/paraview_script.py
