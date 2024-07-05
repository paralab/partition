# File containing list of mesh files
file_list_file="/home/budvin/research/Partitioning/paralab-partition/misc/connected_meshes_tet.txt"

# Read the file list into an array, skipping empty lines
mapfile -t mesh_file_list < <(grep -v '^$' "$file_list_file")


# for file_idx in "${!mesh_file_list[@]}"; do 
for file_idx in {0..20}; do 

    echo "$file_idx copying file ${mesh_file_list[$file_idx]}"
    scp ${mesh_file_list[$file_idx]} budvin@frontera.tacc.utexas.edu:/work2/10000/budvin/frontera/datasets/Meshes/10k_tet/
done