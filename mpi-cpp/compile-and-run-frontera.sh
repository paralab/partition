#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH -n 1280
#SBATCH -N 23
#SBATCH --cpus-per-task 1
#SBATCH -o /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/output.txt
#SBATCH -e /work2/10000/budvin/frontera/partition-project/paralab-partition-repo/mpi-cpp/error.txt

#SBATCH -p development



##SBATCH --mail-user=u1472227@utah.edu  
##SBATCH --mail-type=FAIL




set -e

module load cmake/3.24.2

# module load gcc/9.1.0
# module load mvapich2-x/2.3
# module load python3/3.9.2


module load intel/19.1.1
module load impi/19.0.9
module load python3/3.9.2




export OMP_NUM_THREADS=1

export OMP_PROC_BIND=true 

export OMP_PLACES=cores 

GMSH_SDK_PATH=$WORK/partition-project/dependencies/gmsh-4.13.0-source/build/install
METIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/GKlib/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=$WORK/partition-project/dependencies/ParMETIS/build/Linux-x86_64/build



mkdir -p build

cmake -S . -B build -DGMSH_SDK_PATH=${GMSH_SDK_PATH} \
    -DMETIS_INSTALL_DIR_PATH=${METIS_INSTALL_DIR_PATH} \
    -DGKLIB_INSTALL_DIR_PATH=${GKLIB_INSTALL_DIR_PATH} \
    -DPARMETIS_INSTALL_DIR_PATH=${PARMETIS_INSTALL_DIR_PATH}

make -C ./build

# export LD_LIBRARY_PATH="${VTK_ROOT}/lib:${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${OPENMPI_ROOT}/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${MPI_ROOT}/lib:${LD_LIBRARY_PATH}"

echo -e "compilation done"


# dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
dir=$PWD

# python with pandas needed for metrics export in json format
source "$dir/.venv/bin/activate"
export PYTHONPATH=$PYTHONPATH:$dir

metrics_file_path="$dir/results/frontera-intelmpi-scaling-$(date +%Y-%m-%d__%H-%M-%S).json"

echo "exporting metrics to file $metrics_file_path"

# File containing list of mesh files
file_list_file="$dir/mesh_files_list_scaling.txt"

# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")




for file_idx in "${!mesh_file_list[@]}"; do 
    for np in 10 20 40 80 160 320 640 1280
    do
        for run_idx in {0..4}; do
            ibrun -np $np ./build/main-new $WORK/datasets/Meshes/${mesh_file_list[$file_idx]} $file_idx $run_idx $metrics_file_path -no-viz
            sleep 5s
        done
    done
done

# ibrun -np 4 ./build/main-new $WORK/datasets/Meshes/hex-box-5x5x2.msh 0 0 $metrics_file_path -no-viz


# export SFC_morton="$PWD/out-sfc.vtk"
# export METIS="$PWD/out-parmetis.vtk"
# export BFS="$PWD/out-bfs.vtk"
# export BFS_grow="$PWD/out-bfs.vtk"

# /home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ../grow/paraview_script.py

echo "=====done======"
