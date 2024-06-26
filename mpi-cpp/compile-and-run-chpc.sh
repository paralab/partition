#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --ntasks=320
#SBATCH --cpus-per-task=1
#SBATCH -o /uufs/chpc.utah.edu/common/home/u1472227/partition-project/paralab-partition-repo/mpi-cpp/output.txt
#SBATCH -e /uufs/chpc.utah.edu/common/home/u1472227/partition-project/paralab-partition-repo/mpi-cpp/error.txt
#SBATCH --account=soc-np
#SBATCH --partition=soc-np
#SBATCH --mem-per-cpu=4GB


#---------- SBATCH --exclude=notch[385-389,282-292]




set -e

module load cmake/3.26.0
module load ninja
module load ccache
module load gcc/11.2.0
module load openmpi/4.1.4
module load vtk/9.2.6
module load python/3.10


export OMP_NUM_THREADS=1

export OMP_PROC_BIND=true 

export OMP_PLACES=cores 

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

# dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"
dir=$PWD

# python with pandas needed for metrics export in json format
source "$dir/../.venv/bin/activate"
export PYTHONPATH=$PYTHONPATH:$dir

metrics_file_path="$dir/results/chpc-scaling-$(date +%Y-%m-%d__%H-%M-%S).json"

echo "exporting metrics to file $metrics_file_path"

# File containing list of mesh files
file_list_file="$dir/mesh_files_list_scaling.txt"

# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")

# mpirun -np 320 ./build/main-new /scratch/general/vast/u1472227/meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh  #smallest
# mpirun -np 25 ./build/main-new /scratch/general/vast/u1472227/meshes/10k_tet/196209_sf_hexa.mesh_73346_289961.obj.mesh   #large tet
# mpirun -np 320 ./build/main-new /scratch/general/vast/u1472227/meshes/10k_hex/75651_sf_hexa.mesh   #largest hex
# mpirun -np 320 ./build/main-new /scratch/general/vast/u1472227/meshes/10k_tet/75651_sf_hexa.mesh_78608_298692.obj.mesh   #largest tet

# mpirun -np 80 ./build/main-new /scratch/general/vast/u1472227/meshes/10k_tet/63461_sf_hexa.mesh_44484_213801.obj.mesh 0 $metrics_file_path   #1.1M tet

# /scratch/general/vast/u1472227/meshes/10k_tet/472002_sf_hexa.mesh_44020_172566.obj.mesh     # 843,280 tetrahedra

# for i in 10 20 40 80 160 320 640 1280
# do
#    mpirun -np $i ./build/main-new /scratch/general/vast/u1472227/meshes/10k_tet/75651_sf_hexa.mesh_78608_298692.obj.mesh  # largest tet 1,598,423
# done


for file_idx in "${!mesh_file_list[@]}"; do 
    for np in 10 20 40 80 160 320
    do
        mpirun -np $np ./build/main-new /scratch/general/vast/u1472227/meshes/${mesh_file_list[$file_idx]} $file_idx $metrics_file_path -no-viz
    done
done


# export SFC_morton="$PWD/out-sfc.vtk"
# export METIS="$PWD/out-parmetis.vtk"
# export BFS="$PWD/out-bfs.vtk"
# export BFS_grow="$PWD/out-bfs.vtk"

# /home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ../grow/paraview_script.py

echo "=====done======"
