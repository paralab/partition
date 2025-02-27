set -e

output_dir=cut-cell-meshes

mkdir -p $output_dir

MAX_DIVISIONS_PER_DIM=40

# File containing list of unstructured mesh files
file_list_file="/home/budvin/research/Partitioning/siam-cut-cell-presentation/mpi-cpp/meshes-for-cut-cell.txt"

# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")

mesh_path_prefix="/home/budvin/research/Partitioning/Meshes/"

set +e

# for file_idx in {0..1}; do
for file_idx in "${!mesh_file_list[@]}"; do 
    file_path=${mesh_path_prefix}${mesh_file_list[$file_idx]}
    echo "starting embedding for ${file_path}"
    ./build/embed_mesh $file_path ${output_dir}/embedded_mesh_${file_idx}.msh ${MAX_DIVISIONS_PER_DIM}
    echo "done embedding"
done

echo "====== EMBEDDING TASK FINISHED =========="

