import pandas as pd
# from ctypes import c_uint32, c_uint64

import typing

def export_metrics(mesh_file: str, file_idx: int, partition_count: int, vertex_count: int,
                    sfc_partition_sizes: list[int], sfc_partition_boundaries: list[int],
                    bfs_partition_sizes: list[int], bfs_partition_boundaries: list[int], bfs_time: int,
                    grow_partition_sizes: list[int], grow_partition_boundaries: list[int], grow_time: int,
                    parmetis_partition_sizes: list[int], parmetis_partition_boundaries: list[int], parmetis_time: int,
                    metrics_out_file_path: str):
    method_names = ['SFC_morton','BFS','BFS_grow','METIS']
    times = [0, bfs_time, grow_time, parmetis_time]
    result_row = pd.Series()
    result_row['mesh_idx'] = file_idx
    result_row['mesh_file'] = mesh_file
    result_row['np'] = partition_count
    result_row['n'] = vertex_count

    for (sizes,boundaries), time_,  method in zip([
        [sfc_partition_sizes,sfc_partition_boundaries],
        [bfs_partition_sizes,bfs_partition_boundaries],
        [grow_partition_sizes,grow_partition_boundaries],
        [parmetis_partition_sizes,parmetis_partition_boundaries]],times, method_names):
        total_boundary_size = sum(boundaries)
        max_part_size = max(sizes)
        min_part_size = min(sizes)
        ideal_size = vertex_count//partition_count

        result_row[f"{method}_boundary_ratio"] = int(total_boundary_size)/int(vertex_count)
        result_row[f"{method}_boundary_ratio_expr"] = f"{total_boundary_size}/{vertex_count}"
        result_row[f"{method}_rho_max"] = int(max_part_size)/int(ideal_size)
        result_row[f"{method}_rho_max_expr"] = f"{max_part_size}/{ideal_size}"
        result_row[f"{method}_rho_min"] = int(min_part_size)/int(ideal_size)
        result_row[f"{method}_rho_min_expr"] = f"{min_part_size}/{ideal_size}"
        result_row[f"{method}_partition_sizes"] = sizes
        result_row[f"{method}_partition_boundaries"] = boundaries
        result_row[f"{method}_time"] = time_


    pd.DataFrame([result_row]).to_json(metrics_out_file_path,index=False,mode='a',lines=True,orient='records')
