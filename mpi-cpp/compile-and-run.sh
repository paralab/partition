set -e

export OMP_NUM_THREADS=1

# export OMP_PROC_BIND=true 

# export OMP_PLACES=threads 


GMSH_SDK_PATH=/home/budvin/bin/gmsh-4.13.0-source/build/install
METIS_INSTALL_DIR_PATH=/home/budvin/bin/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=/home/budvin/bin/GKlib-master/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=/home/budvin/bin/ParMETIS/build/Linux-x86_64/build
PETSC_INSTALL_DIR_PATH=/home/budvin/bin/petsc/build/install


mkdir -p build

cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DGMSH_SDK_PATH=${GMSH_SDK_PATH} \
    -DENABLE_VTK_FEATURES=ON \
    -DVTK_INSTALL_DIR_PATH=/home/budvin/bin/VTK-9.3.0/build \
    -DMETIS_INSTALL_DIR_PATH=${METIS_INSTALL_DIR_PATH} \
    -DGKLIB_INSTALL_DIR_PATH=${GKLIB_INSTALL_DIR_PATH} \
    -DPARMETIS_INSTALL_DIR_PATH=${PARMETIS_INSTALL_DIR_PATH} \
    -DPETSC_INSTALL_DIR_PATH=${PETSC_INSTALL_DIR_PATH} \
    -DGRAPH_INDEXING_TYPE=32 -DBFS_DISTANCE_TYPE=32 -DBFS_LABEL_TYPE=16

ninja -C ./build

echo "====== complation done ==============="


export LD_LIBRARY_PATH="${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${PETSC_INSTALL_DIR_PATH}/lib:${LD_LIBRARY_PATH}"


# python with pandas needed for metrics export in json format
source ../.venv/bin/activate
dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
export PYTHONPATH=$PYTHONPATH:$dir

metrics_file_path="$PWD/results/$(date +%Y-%m-%d__%H-%M-%S).json"

echo "exporting metrics to file $metrics_file_path"

# File containing list of mesh files
file_list_file="/home/budvin/research/Partitioning/paralab-partition/mpi-cpp/mesh_files_list.txt"

# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")

# mpirun -np 40 --oversubscribe ./build/main-new

# mpirun -np 10 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh  # octopus
# mpirun -np 8 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh
# mpirun -np 40 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh 0 $metrics_file_path  # smallest tet
# mpirun -np 8 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/10k_tet/196209_sf_hexa.mesh_73346_289961.obj.mesh #large tet

# /home/budvin/research/Partitioning/mesh_generator/hex-box-23x23x23.msh

# mpirun -np 40 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/10k_tet/67923_sf_hexa.mesh_2992_10000.obj.mesh 1 $metrics_file_path # small tet

# for file_idx in "${!mesh_file_list[@]}"; do 
#     for np in 4 8
#     do
#         for run_idx in {0..2}; do
#             mpirun -np $np --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/${mesh_file_list[$file_idx]} $file_idx $run_idx $metrics_file_path -no-viz
#         done
#     done
# done

parts_n=4
mesh_file=/home/budvin/research/Partitioning/mesh_generator/hex-box-23x23x23.msh
part_file_prefix="$dir/tmp/part"


time ./build/gmsh_partition $mesh_file $parts_n $part_file_prefix

# exit 0

mpirun -np $parts_n --oversubscribe ./build/main-new $mesh_file $part_file_prefix 0 0 $dir/test.json -viz

rm "${part_file_prefix}"*

# mpirun -np 4 --oversubscribe ./build/main-new /home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh 0 0 $dir/test.json -viz



export SFC_morton="$PWD/out-sfc.vtk"
export METIS="$PWD/out-parmetis.vtk"
export BFS="$PWD/out-bfs.vtk"
export BFS_grow="$PWD/out-grow.vtk"

/home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ../misc/grow/paraview_script.py
