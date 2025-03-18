import functools
import numpy as np
import os
os.environ["METIS_DLL"] = "/home/budvin/research/Partitioning/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.so"
import metis

import subprocess

from gmsh_utils import read_mesh_file
from sfc_utils import to_integer_coords, morton_compare
from vtk_utils import export_points_with_parts_to_vtk, export_points_scalar_to_vtk
from misc_utils import part_counts, print_metrics
from BFS_Partition import get_BFS_partitions
from diffusion_partition import diffuse_from_seeds, diffuse_from_partitions, refine_by_diffusion



def run_single_file(file_name: str, parts_n: int):
    data = read_mesh_file(file_name)
    assert(data["valid"])
    idx_to_element = data["idx_to_element"]
    element_to_idx = data["element_to_idx"]
    e2e_graph = data["e2e_graph"]
    element_coords = data["element_coords"]

    element_count = len(element_coords)

    integer_element_coords = to_integer_coords(element_coords)
    integer_element_coords_with_index = [ [coord,i] for (i,coord) in enumerate(integer_element_coords)]
    integer_element_coords_with_index_morton_order = sorted(integer_element_coords_with_index, key=functools.cmp_to_key(morton_compare))
    morton_order = [d[1] for d in integer_element_coords_with_index_morton_order]   ## morton order in indices

    morton_sfc_part_sizes = part_counts(element_count, parts_n)
    morton_sfc_partition_labels = [None for _ in range(element_count)]
    morton_i = 0
    for p_i, p_c in enumerate(morton_sfc_part_sizes):
        for _ in range(p_c):
            morton_sfc_partition_labels[morton_order[morton_i]] = p_i
            morton_i+=1
    print("SFC partition done")
    
    morton_sfc_part_sizes_scanned = [0] + np.cumsum(morton_sfc_part_sizes).tolist()
    morton_seeds_indices = [morton_sfc_part_sizes_scanned[p_i] + morton_sfc_part_sizes[p_i]//2   for p_i in range(parts_n)]
    morton_seeds = [idx_to_element[morton_order[s_i]] for s_i in morton_seeds_indices]

    element_to_BFS_partition = get_BFS_partitions(e2e_graph,morton_seeds ,parts_n)

    BFS_partition_labels = [None for _ in range(element_count)]

    for elem in element_to_BFS_partition:
        BFS_partition_labels[element_to_idx[elem]] = element_to_BFS_partition[elem]

    assert(None not in BFS_partition_labels)
    print("BFS partition done")

    # element_to_diffusion_value = diffuse_from_partitions(e2e_graph, element_to_BFS_partition, parts_n)


    # export_points_with_parts_to_vtk(element_coords,BFS_partition_labels, "BFS.vtk")
    # exit()    

    element_to_fastpart_partition, element_to_diffusion_value  = refine_by_diffusion(e2e_graph,element_to_BFS_partition, parts_n)

    fastpart_partition_labels = [None for _ in range(element_count)]

    for elem in element_to_fastpart_partition:
        fastpart_partition_labels[element_to_idx[elem]] = element_to_fastpart_partition[elem]

    assert(None not in element_to_fastpart_partition)
    print("fastpart partition done")

    diffusion_values = [None for _ in range(element_count)]

    for elem in element_to_diffusion_value:
        diffusion_values[element_to_idx[elem]] = element_to_diffusion_value[elem]
    export_points_scalar_to_vtk(element_coords, diffusion_values, "diffusion.vtk")


    (_, METIS_partition_labels) = metis.part_graph(e2e_graph, parts_n)
    assert(None not in METIS_partition_labels)
    print("METIS partition done")

    print()
    print("SFC metrics")
    print_metrics(e2e_graph, morton_sfc_partition_labels, element_to_idx, parts_n)
    print()
    print("BFS metrics")
    print_metrics(e2e_graph, BFS_partition_labels, element_to_idx, parts_n)
    print()
    print("fastpart metrics")
    print_metrics(e2e_graph, fastpart_partition_labels, element_to_idx, parts_n)
    print()
    print("METIS metrics")
    print_metrics(e2e_graph, METIS_partition_labels, element_to_idx, parts_n)



    export_points_with_parts_to_vtk(element_coords,morton_sfc_partition_labels, "SFC_morton.vtk")
    export_points_with_parts_to_vtk(element_coords,BFS_partition_labels, "BFS.vtk")
    export_points_with_parts_to_vtk(element_coords,fastpart_partition_labels, "fastpart.vtk")
    export_points_with_parts_to_vtk(element_coords,METIS_partition_labels, "METIS.vtk")

    my_env = os.environ.copy()
    # method_names = ['SFC_morton','METIS','BFS','fastpart']


    base_dir = os.getcwd()
    my_env["SFC_morton"] = base_dir + "/SFC_morton.vtk"
    my_env["METIS"] = base_dir + "/METIS.vtk"
    my_env["BFS"] = base_dir + "/BFS.vtk"
    my_env["fastpart"] = base_dir + "/fastpart.vtk"


    subprocess.run(['/home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview','paraview_script.py'],env=my_env)


    return


def run_multiple(file_names: list[str]):
    return

# input_file_name = "/home/budvin/research/Partitioning/mesh_generator/hex-box-50x50x2.msh"
input_file_name = "/home/budvin/research/Partitioning/mesh_generator/generated_tet_50x50x2.mesh"
input_file_name = "/home/budvin/research/Partitioning/Meshes/10k_tet/57181_sf_hexa.mesh_5006_17194.obj.mesh"
input_file_name = "/home/budvin/research/Partitioning/Meshes/10k_tet/130968_sf_hexa.mesh_9940_38180.obj.mesh"       # disc
input_file_name = "/home/budvin/research/Partitioning/Meshes/10k_tet/90280_sf_hexa.mesh_6608_22267.obj.mesh"
# input_file_name = "/home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh"        # octopus
run_single_file(input_file_name, 11)


