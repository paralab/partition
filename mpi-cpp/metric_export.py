import pandas as pd
# from ctypes import c_uint32, c_uint64

import typing

def export_metrics(mesh_file: str, file_idx: int, run_idx: int, partition_count: int, vertex_count: int,
                    graph_setup_time : int,
                    sfc_partition_sizes: list, sfc_partition_boundaries: list, sfc_partition_time: int, sfc_mat_assembly_time: int, sfc_matvec_time: int,
                    bfs_partition_sizes: list, bfs_partition_boundaries: list, bfs_labeling_time: int, bfs_redistribution_time: int, bfs_mat_assembly_time: int, bfs_matvec_time: int,
                    parmetis_partition_sizes: list, parmetis_partition_boundaries: list, parmetis_labeling_time: int, parmetis_redistribution_time: int, parmetis_mat_assembly_time: int, parmetis_matvec_time: int,
                    ptscotch_partition_sizes: list, ptscotch_partition_boundaries: list, ptscotch_labeling_time: int, ptscotch_redistribution_time: int, ptscotch_mat_assembly_time: int, ptscotch_matvec_time: int,
                    metrics_out_file_path: str):

    result_row = pd.Series()
    result_row['mesh_idx'] = file_idx
    result_row['run_idx'] = run_idx
    result_row['mesh_file'] = mesh_file
    result_row['np'] = partition_count
    result_row['n'] = vertex_count


    method_names = ['SFC_morton','BFS','parMETIS', 'ptscotch']
    mat_assembly_times = [sfc_mat_assembly_time, bfs_mat_assembly_time, parmetis_mat_assembly_time, ptscotch_mat_assembly_time]
    matvec_times = [sfc_matvec_time, bfs_matvec_time, parmetis_matvec_time, ptscotch_matvec_time]

    for (sizes,boundaries), mat_assembly_time, matvec_time, method in zip([
        [sfc_partition_sizes,sfc_partition_boundaries],
        [bfs_partition_sizes,bfs_partition_boundaries],
        [parmetis_partition_sizes,parmetis_partition_boundaries],
        [ptscotch_partition_sizes,ptscotch_partition_boundaries]], mat_assembly_times, matvec_times, method_names):
        
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

        result_row[f"{method}_mat_assembly_time"] = mat_assembly_time
        result_row[f"{method}_matvec_time"] = matvec_time



    result_row['graph_setup_time'] = graph_setup_time

    result_row['SFC_morton_partition_time'] = sfc_partition_time

    result_row['BFS_labeling_time'] = bfs_labeling_time
    result_row['parMETIS_labeling_time'] = parmetis_labeling_time
    result_row['ptscotch_labeling_time'] = ptscotch_labeling_time


    result_row['BFS_redistribution_time'] = bfs_redistribution_time
    result_row['parMETIS_redistribution_time'] = parmetis_redistribution_time
    result_row['ptscotch_redistribution_time'] = ptscotch_redistribution_time




    pd.DataFrame([result_row]).to_json(metrics_out_file_path,index=False,mode='a',lines=True,orient='records')
