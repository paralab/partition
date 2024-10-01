#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 320
#SBATCH -N 6
#SBATCH -o /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/output.txt
#SBATCH -e /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/error.txt

#SBATCH -p development



#SBATCH --mail-user=u1472227@utah.edu  
#SBATCH --mail-type=FAIL




set -e

module load cmake/3.24.2


module load intel/23.1.0
module load impi/21.9.0
# module load python3/3.9.2
module load petsc/3.21
module load fftw3/3.3.10


export OMP_NUM_THREADS=1

export OMP_PROC_BIND=true 

export OMP_PLACES=cores 



GMSH_SDK_PATH=$WORK/partition-project/dependencies/gmsh-repo/build/install
METIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/ParMETIS/build/Linux-x86_64/build
SCOTCH_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/scotch/build/install
PETSC_INSTALL_DIR_PATH=${PETSC_DIR}


mkdir -p build

cmake -S . -B build -DCMAKE_C_COMPILER=mpiicc -DCMAKE_CXX_COMPILER=mpiicpc \
    -DGMSH_SDK_PATH=${GMSH_SDK_PATH} \
    -DENABLE_VTK_FEATURES=OFF \
    -DENABLE_PETSC_DRAW_FEATURES=OFF \
    -DMETIS_INSTALL_DIR_PATH=${METIS_INSTALL_DIR_PATH} \
    -DGKLIB_INSTALL_DIR_PATH=${GKLIB_INSTALL_DIR_PATH} \
    -DPARMETIS_INSTALL_DIR_PATH=${PARMETIS_INSTALL_DIR_PATH} \
    -DPETSC_INSTALL_DIR_PATH=${PETSC_INSTALL_DIR_PATH} \
    -DSCOTCH_INSTALL_DIR_PATH=${SCOTCH_INSTALL_DIR_PATH} \
    -DGRAPH_INDEXING_TYPE=32 -DBFS_DISTANCE_TYPE=32 -DBFS_LABEL_TYPE=16

make -C ./build

export LD_LIBRARY_PATH="${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${PETSC_INSTALL_DIR_PATH}/clx/lib:${TACC_FFTW3_LIB}:${MPI_ROOT}/lib:${SCOTCH_INSTALL_DIR_PATH}/lib:${LD_LIBRARY_PATH}"

echo -e "compilation done"

# exit 0


dir=$PWD

# python with pandas needed for metrics export in json format
# source "$dir/.venv/bin/activate"
# export PYTHONPATH=$PYTHONPATH:$dir

# metrics_file_path="$dir/results/frontera-scaling-front_bfs-200spmv-largetet-5rfnrnds___$(date +%Y-%m-%d__%H-%M-%S).json"

# metrics_file_path="$dir/results/test.json"

# echo "exporting metrics to file $metrics_file_path"

# File containing list of mesh files
# file_list_file="$dir/mesh_files_list.txt"
file_list_file="$dir/connected_hex.txt"


# file_list_file="$dir/largest_hex_files.txt"


# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")





# for file_idx in "${!mesh_file_list[@]}"; do 
#     for np in 10 20 40 80 160 320 640 1280 2240
#     do
#         for run_idx in {0..2}; do
#             set +e
#             time ibrun -np $np ./build/main-new $SCRATCH/meshes/${mesh_file_list[$file_idx]} $file_idx $run_idx $metrics_file_path -no-viz < /dev/null
#             set -e
#             sleep 2s
#         done
#     done
# done




# test_mesh_file=$SCRATCH/meshes/10k_tet/75651_sf_hexa.mesh_78608_298692.obj.mesh


# for np in 10 20 40 80 160 320 640 1280 2240
# do
#     for run_idx in {0..2}; do
#         set +e
#         time ibrun -np $np ./build/main-new $test_mesh_file 0 $run_idx $dir/test-new-sort.json -no-viz < /dev/null
#         set -e
#         sleep 2s
#     done 
# done

out_prefix="$(date +%Y-%m-%d__%H-%M-%S)_1rnds_test_grain_hex"
# out_prefix="2024-09-20__01-40-56_5rnds_test_grain_tet"
echo "out_prefix: $out_prefix"
grain_sizes=( 800 1200 1600 )
for file_idx in "${!mesh_file_list[@]}"; do 
# for file_idx in {0..4}; do 
# for ((file_idx=57; file_idx<${#mesh_file_list[@]}; file_idx++)); do 

    file_path=$SCRATCH/meshes/${mesh_file_list[$file_idx]}
    
    for g in "${grain_sizes[@]}"
    do
        np=$(./build/process_count $file_path $g)
        echo "[$file_idx]  partitioning $file_path for grain size $g with $np processes"
        for run_idx in {0..2}; do
            set +e
            time ibrun -np $np ./build/main-new $file_path $file_idx $run_idx $dir/results/${out_prefix}_$g.json  -no-viz < /dev/null
            set -e
            sleep 2s
        done
    done
done


echo "=====done======"
