#!/bin/bash
#SBATCH -t 0:30:00
#SBATCH -n 1280
#SBATCH -N 23
#SBATCH -o /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/output.txt
#SBATCH -e /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/error.txt

#SBATCH -p development



#SBATCH --mail-user=u1472227@utah.edu  
#SBATCH --mail-type=FAIL




set -e

module load cmake/3.24.2


module load intel/19.1.1
module load impi/19.0.9
module load python3/3.9.2
module load petsc/3.19
module load fftw3/3.3.10



export OMP_NUM_THREADS=1

export OMP_PROC_BIND=true 

export OMP_PLACES=cores 



GMSH_SDK_PATH=$WORK/partition-project/dependencies/gmsh-4.13.0-source/build/install
METIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/ParMETIS/build/Linux-x86_64/build
SCOTCH_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/scotch/build/install
PETSC_INSTALL_DIR_PATH=${PETSC_DIR}


mkdir -p build

cmake -S . -B build -DGMSH_SDK_PATH=${GMSH_SDK_PATH} \
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


# dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
dir=$PWD

# python with pandas needed for metrics export in json format
source "$dir/.venv/bin/activate"
export PYTHONPATH=$PYTHONPATH:$dir

metrics_file_path="$dir/results/frontera-scaling-front_bfs-directsolver___$(date +%Y-%m-%d__%H-%M-%S).json"

# metrics_file_path="$dir/results/test.json"

echo "exporting metrics to file $metrics_file_path"

# File containing list of mesh files
file_list_file="$dir/mesh_file_list_scaling.txt"


# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")


part_file_prefix="/scratch1/10000/budvin/partitioning-temp/part"



# for file_idx in "${!mesh_file_list[@]}"; do 
#     for np in 10 20 40 80 160 320 640 1280
#     do
#         time ./build/gmsh_partition $WORK/datasets/Meshes/${mesh_file_list[$file_idx]} $np $part_file_prefix
#         for run_idx in {0..3}; do
#             set +e
#             ibrun -np $np ./build/main-new $WORK/datasets/Meshes/${mesh_file_list[$file_idx]} $part_file_prefix $file_idx $run_idx $metrics_file_path -no-viz < /dev/null
#             set -e
#             sleep 3s
#         done
#         rm "${part_file_prefix}"*
#     done
# done




test_mesh_file=$WORK/datasets/Meshes/10k_tet/75651_sf_hexa.mesh_78608_298692.obj.mesh

for np in 1280
do
    ./build/gmsh_partition $test_mesh_file $np $part_file_prefix
    for run_idx in {0..2}; do
        set +e
        ibrun -np $np ./build/main-new $test_mesh_file $part_file_prefix 0 $run_idx $dir/results/test.json -no-viz < /dev/null
        set -e
        sleep 2s
    done 
    rm "${part_file_prefix}"*
done



echo "=====done======"
