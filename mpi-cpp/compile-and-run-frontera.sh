#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 2240
#SBATCH -N 40
#SBATCH -o /work2/10000/budvin/frontera/partition-project/diffusion-partition-repo/mpi-cpp/bfs_rounds_count_output.txt
#SBATCH -e /work2/10000/budvin/frontera/partition-project/diffusion-partition-repo/mpi-cpp/bfs_rounds_count_error.txt

#SBATCH -p development



#SBATCH --mail-user=budvin.edippuliarachchi@tufts.edu
#SBATCH --mail-type=FAIL




set -e

module load cmake/3.24.2


module load intel/23.1.0
module load impi/21.9.0
module load petsc/3.21
module load fftw3/3.3.10


export OMP_NUM_THREADS=1
export SCOTCH_PTHREAD_NUMBER=1

export OMP_PROC_BIND=true 

export OMP_PLACES=cores 



GMSH_SDK_PATH=$WORK/partition-project/dependencies/gmsh-repo/build/install
METIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/ParMETIS/build/Linux-x86_64/build
SCOTCH_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/scotch/build/install
AMRACUT_INSTALL_DIR_PATH=$WORK/partition-project/amracut/build/install
PETSC_INSTALL_DIR_PATH=${PETSC_DIR}


root_dir=$PWD

cd $WORK/partition-project/amracut

bash compile.sh

cd $root_dir


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
    -DAMRACUT_INSTALL_DIR_PATH=${AMRACUT_INSTALL_DIR_PATH} \
    -DGRAPH_INDEXING_TYPE=32 -DBFS_DISTANCE_TYPE=32 -DBFS_LABEL_TYPE=16

make -C ./build

export LD_LIBRARY_PATH="${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${PETSC_INSTALL_DIR_PATH}/clx/lib:${TACC_FFTW3_LIB}:${MPI_ROOT}/lib:${SCOTCH_INSTALL_DIR_PATH}/lib:${AMRACUT_INSTALL_DIR_PATH}/lib:${LD_LIBRARY_PATH}"

echo -e "===== compilation done ====="

# exit 0


dir=$PWD


# metrics_file_path="$dir/results/diffusion_sc25_cutcells_w-edges_tet_meshes2025-03-26__10-58-34.json"
# metrics_file_path="$dir/results/diffusion_sc25_unweighted_tet_meshes2025-03-26__23-25-06.json"
metrics_file_path="$dir/results/tmp.json"


echo "exporting metrics to file $metrics_file_path"

# File containing list of mesh files
file_list_file="$dir/connected_tet.txt"
# file_list_file="$dir/octree_files.txt"




# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")





# for file_idx in "${!mesh_file_list[@]}"; do 
for ((file_idx=0; file_idx<10; file_idx++)); do
    for np in 10 20 40 80 160 320 640 1280 2240
    do
        for run_idx in {0..1}; do
            set +e
            time ibrun -np $np ./build/main-new $SCRATCH/meshes/${mesh_file_list[$file_idx]} $file_idx $run_idx $metrics_file_path -no-viz < /dev/null
            # time ibrun -np $np ./build/main-octree $SCRATCH/meshes/octree/${mesh_file_list[$file_idx]} $file_idx $run_idx $metrics_file_path -no-viz < /dev/null
            
            set -e
            sleep 2s
        done
    done
done



echo "=====done======"
