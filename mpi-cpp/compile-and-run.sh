set -e

root_dir=$PWD

# first compile AMRaCut

cd /home/budvin/research/Partitioning/AMRaCut

bash compile.sh

cd $root_dir



export OMP_NUM_THREADS=1


GMSH_SDK_PATH=/home/budvin/bin/gmsh-4.13.0-source/build/install
METIS_INSTALL_DIR_PATH=/home/budvin/bin/METIS-5.2.1/build/build
GKLIB_INSTALL_DIR_PATH=/home/budvin/bin/GKlib-master/build/Linux-x86_64/build
PARMETIS_INSTALL_DIR_PATH=/home/budvin/bin/ParMETIS/build/Linux-x86_64/build
PETSC_INSTALL_DIR_PATH=/home/budvin/bin/petsc/build/install
SCOTCH_INSTALL_DIR_PATH=/home/budvin/bin/scotch/build/install
AMRACUT_INSTALL_DIR_PATH=/home/budvin/research/Partitioning/AMRaCut/build/install
VTK_INSTALL_DIR_PATH=/home/budvin/bin/VTK-9.3.0/build

mkdir -p build

cmake -G Ninja -S . -B build -DCMAKE_CXX_COMPILER_LAUNCHER=ccache -DGMSH_SDK_PATH=${GMSH_SDK_PATH} \
    -DENABLE_VTK_FEATURES=ON \
    -DENABLE_PETSC_DRAW_FEATURES=OFF \
    -DVTK_INSTALL_DIR_PATH=${VTK_INSTALL_DIR_PATH} \
    -DMETIS_INSTALL_DIR_PATH=${METIS_INSTALL_DIR_PATH} \
    -DGKLIB_INSTALL_DIR_PATH=${GKLIB_INSTALL_DIR_PATH} \
    -DPARMETIS_INSTALL_DIR_PATH=${PARMETIS_INSTALL_DIR_PATH} \
    -DPETSC_INSTALL_DIR_PATH=${PETSC_INSTALL_DIR_PATH} \
    -DSCOTCH_INSTALL_DIR_PATH=${SCOTCH_INSTALL_DIR_PATH} \
    -DAMRACUT_INSTALL_DIR_PATH=${AMRACUT_INSTALL_DIR_PATH} \
    -DGRAPH_INDEXING_TYPE=32 -DBFS_DISTANCE_TYPE=32 -DBFS_LABEL_TYPE=16

ninja -C ./build

echo "====== complation done ==============="

# exit 0

export LD_LIBRARY_PATH="${GMSH_SDK_PATH}/lib:${METIS_INSTALL_DIR_PATH}/lib:${GKLIB_INSTALL_DIR_PATH}/lib:${PARMETIS_INSTALL_DIR_PATH}/lib:${PETSC_INSTALL_DIR_PATH}/lib:${SCOTCH_INSTALL_DIR_PATH}/lib:${AMRACUT_INSTALL_DIR_PATH}/lib:${LD_LIBRARY_PATH}"

dir="$( dirname -- "$( readlink -f -- "$0"; )"; )"



metrics_file_path="$PWD/results/$(date +%Y-%m-%d__%H-%M-%S).json"

echo "exporting metrics to file $metrics_file_path"

# File containing list of mesh files
file_list_file="$PWD/connected_tet.txt"

# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")




parts_n=12


mesh_file="/home/budvin/research/Partitioning/mesh_generator/generated_tet_100x100x2.mesh"

mpirun -np $parts_n --oversubscribe ./build/main-new $mesh_file 0 0 $dir/tmp.json -viz
# mpirun -np $parts_n --oversubscribe ./build/main-octree $mesh_file 0 0 $dir/tmp.json -viz




export SFC_morton="$PWD/out-sfc.vtk"
export parMETIS="$PWD/out-parmetis.vtk"
export amracut="$PWD/out-amracut.vtk"
export ptscotch="$PWD/out-ptscotch.vtk"

/home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview ./paraview_script.py
