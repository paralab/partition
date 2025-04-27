#!/bin/bash
#SBATCH -t 48:00:00
#SBATCH -n 672
#SBATCH -N 12
#SBATCH -o /work2/10000/budvin/frontera/partition-project/partition-fixed-g-size-repo/mpi-cpp/2400g-output.txt
#SBATCH -e /work2/10000/budvin/frontera/partition-project/partition-fixed-g-size-repo/mpi-cpp/2400g-error.txt

#SBATCH -p normal



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

echo -e "===== compilation done ====="

export LD_LIBRARY_PATH="${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${PETSC_INSTALL_DIR_PATH}/clx/lib:${TACC_FFTW3_LIB}:${MPI_ROOT}/lib:${SCOTCH_INSTALL_DIR_PATH}/lib:${AMRACUT_INSTALL_DIR_PATH}/lib:${LD_LIBRARY_PATH}"


# exit 0


dir=$PWD


file_list_file="$dir/connected_tet.txt"


# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")





out_prefix="tet_unweighted_fixed_g_size"
echo "out_prefix: $out_prefix"
grain_size="2400"
for ((file_idx=0; file_idx<${#mesh_file_list[@]}; file_idx++)); do 

    file_path=$SCRATCH/meshes/${mesh_file_list[$file_idx]}
    

    np=$(./build/process_count $file_path $grain_size)
    echo "[$file_idx]  partitioning $file_path for grain size $grain_size with $np processes"
    for run_idx in {0..3}; do
        set +e
        time ibrun -np $np ./build/main-new $file_path $file_idx $run_idx $dir/results/${out_prefix}_$grain_size.json  -no-viz < /dev/null
        set -e
        sleep 2s
    done
    
done


echo "=====done======"