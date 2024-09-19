set -e

export OMP_NUM_THREADS=1

# export OMP_PROC_BIND=true 

# export OMP_PLACES=threads 


GMSH_SDK_PATH=/home/budvin/bin/gmsh-4.13.0-source/build/install
METIS_INSTALL_DIR_PATH=/home/budvin/bin/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=/home/budvin/bin/GKlib-master/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=/home/budvin/bin/ParMETIS/build/Linux-x86_64/build
PETSC_INSTALL_DIR_PATH=/home/budvin/bin/petsc/build/install
SCOTCH_INSTALL_DIR_PATH=/home/budvin/bin/scotch/build/install


mkdir -p build

cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DGMSH_SDK_PATH=${GMSH_SDK_PATH} \
    -DENABLE_VTK_FEATURES=ON \
    -DENABLE_PETSC_DRAW_FEATURES=OFF \
    -DVTK_INSTALL_DIR_PATH=/home/budvin/bin/VTK-9.3.0/build \
    -DMETIS_INSTALL_DIR_PATH=${METIS_INSTALL_DIR_PATH} \
    -DGKLIB_INSTALL_DIR_PATH=${GKLIB_INSTALL_DIR_PATH} \
    -DPARMETIS_INSTALL_DIR_PATH=${PARMETIS_INSTALL_DIR_PATH} \
    -DPETSC_INSTALL_DIR_PATH=${PETSC_INSTALL_DIR_PATH} \
    -DSCOTCH_INSTALL_DIR_PATH=${SCOTCH_INSTALL_DIR_PATH} \
    -DGRAPH_INDEXING_TYPE=32 -DBFS_DISTANCE_TYPE=32 -DBFS_LABEL_TYPE=16

ninja -C ./build

echo "====== complation done ==============="


export LD_LIBRARY_PATH="${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${PETSC_INSTALL_DIR_PATH}/lib:${SCOTCH_INSTALL_DIR_PATH}/lib:${LD_LIBRARY_PATH}"

dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"

# # python with pandas needed for metrics export in json format
# source ../.venv/bin/activate
# 
# export PYTHONPATH=$PYTHONPATH:$dir

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
#     for np in 2 4
#     do
#         for run_idx in {0..2}; do
#             mpirun -np $np --oversubscribe ./build/main-new /home/budvin/research/Partitioning/Meshes/${mesh_file_list[$file_idx]} $file_idx $run_idx $dir/test.json  -no-viz
#         done
#     done
# done



parts_n=20
# mesh_file="/home/budvin/research/Partitioning/mesh_generator/hex-box-60x60x2.msh"
# mesh_file="/home/budvin/research/Partitioning/mesh_generator/generated_tet_50x50x2.mesh"
# mesh_file="/home/budvin/research/Partitioning/Meshes/10k_hex/75651_sf_hexa.mesh"
# mesh_file="/home/budvin/research/Partitioning/Meshes/10k_hex/51508_sf_hexa.mesh"
# mesh_file="/home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh"  #octopus
# mesh_file="/home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh"  #smallest tet


mpirun -np $parts_n --oversubscribe ./build/main-new $mesh_file 0 0 $dir/test.json -viz

# grain_sizes=( 5000 10000 20000 )
# for file_idx in "${!mesh_file_list[@]}"; do 
#     file_path=/home/budvin/research/Partitioning/Meshes/${mesh_file_list[$file_idx]}
    
#     for g in "${grain_sizes[@]}"
#     do
#         np=$(./build/process_count $file_path $g)
#         echo "[$file_idx]  partitioning $file_path for grain size $g with $np processes"
#         for run_idx in {0..2}; do
#             mpirun -np $np --oversubscribe ./build/main-new $file_path $file_idx $run_idx $dir/test_grain_tet_$g.json  -no-viz
#         done
#     done
# done

export SFC_morton="$PWD/out-sfc.vtk"
export parMETIS="$PWD/out-parmetis.vtk"
export BFS="$PWD/out-bfs.vtk"
export ptscotch="$PWD/out-ptscotch.vtk"

/home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ./paraview_script.py
