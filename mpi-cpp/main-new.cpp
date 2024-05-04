#include "util.hpp"
#include "mesh-util.hpp"
#include "vtk-util.hpp"
#include "graph.hpp"
#include "sfc.hpp"
#include "metis-util.hpp"
#include "dist-graph.hpp"
#include <string>
#include <stdexcept>

#include <mpi.h>
#include <algorithm>
#include <chrono>
#include <omp.h>

#include <parUtils.h>
#include <dtypes.h>

#include "ompUtils.h"

int main(int argc, char *argv[])
{
    int numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // uint64_t partition_count=3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(hostname, &len);
    // if (true || taskid == 0)
    // {
    // omp_set_num_threads(8);
    const std::string file_path("/home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh");
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh");
    // const std::string file_path("/home/budvin/research/Partitioning/mesh_generator/hex-box-3x3x3.msh");

    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh");   //octopus
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/196209_sf_hexa.mesh_73346_289961.obj.mesh");     //large tet
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_t/home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.mshet/75651_sf_hexa.mesh_78608_298692.obj.mesh");  //largest tet
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_hex/75651_sf_hexa.mesh");  //largest hex

    ElementType elementType = GetElementType(file_path, MPI_COMM_WORLD);
    print_log("element type: ", elementType);
    uint64_t local_element_count;
    uint64_t global_element_count;

    std::vector<double> global_element_coords;
    std::vector<uint64_t> global_element_partition_labels_morton;

    std::vector<uint64_t> proc_element_counts(numtasks);
    std::vector<uint64_t> proc_element_counts_scanned(numtasks);
    std::vector<ElementWithCoord> local_elements;       // contains {element_tag, global_idx, [x,y,z]}
    std::vector<std::pair<ElementWithTag, ElementWithTag>> local_connected_element_pairs;
    std::vector<std::pair<ElementWithTag, ElementWithTag>> boundary_connected_element_pairs;

    std::vector<ElementWithFace> local_unconnected_elements_faces;



    switch (elementType)
    {
    case ElementType::TET:
    {
        // std::vector<TetElementWithFaces> localElementsAllData;
        // GetElementsWithFacesCentroids<TetElementWithFaces>(file_path, localElementsAllData, ElementType::TET, MPI_COMM_WORLD);
        // SetMortonEncoding(localElementsAllData,ElementType::TET,MPI_COMM_WORLD);
        // std::vector<TetElementWithFaces> localElementsAllDataSorted(localElementsAllData.size());
        // MPI_Barrier(MPI_COMM_WORLD);

        // par::sampleSort<TetElementWithFaces>(localElementsAllData,localElementsAllDataSorted,MPI_COMM_WORLD);
        // MPI_Barrier(MPI_COMM_WORLD);
        // print_log("[", taskid, "]: ", "global SFC sort done");

        // local_element_count = localElementsAllDataSorted.size();
        // local_elements.resize(local_element_count);
        // for (size_t local_elem_i = 0; local_elem_i < local_element_count; local_elem_i++)
        // {
        //     local_elements[local_elem_i].element_tag = localElementsAllDataSorted[local_elem_i].element_tag;
        //     local_elements[local_elem_i].x = localElementsAllDataSorted[local_elem_i].x;
        //     local_elements[local_elem_i].y = localElementsAllDataSorted[local_elem_i].y;
        //     local_elements[local_elem_i].z = localElementsAllDataSorted[local_elem_i].z;


        // }


        // MPI_Gather(&local_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,0,MPI_COMM_WORLD);
        // MPI_Reduce(&local_element_count,&global_element_count,1,MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);

        // ResolveLocalElementConnectivity<TetElementWithFaces>(localElementsAllDataSorted,ElementType::TET,local_connected_element_pairs,local_unconnected_elements_faces);

        // print_log("[", taskid, "]:", "local_connected_element_pairs = ", VectorToString(local_connected_element_pairs));
        // print_log("[", taskid, "]:", "local_unconnected_element_face = ", VectorToString(local_unconnected_elements_faces));





        // auto displacements = GetDisplacementsFromCounts(proc_element_counts);
        // std::vector<int> proc_element_counts_(proc_element_counts.begin(), proc_element_counts.end());
        
        // std::vector<TetElementWithFaces> globalElementsAllData;
        // if (!taskid)
        // {
        //     globalElementsAllData.resize(global_element_count);
        //     global_element_coords.resize(global_element_count*3);
        // }
        // MPI_Gatherv(localElementsAllDataSorted.data(), local_element_count, par::Mpi_datatype<TetElementWithFaces>::value(), 
        //     globalElementsAllData.data(), proc_element_counts_.data(), displacements.data(), par::Mpi_datatype<TetElementWithFaces>::value(), 0, MPI_COMM_WORLD);
        
        // if (!taskid)
        // {
        //     for (size_t global_element_i = 0; global_element_i < global_element_count; global_element_i++)
        //     {
        //         global_element_coords[3*global_element_i] = globalElementsAllData[global_element_i].x;
        //         global_element_coords[3*global_element_i+1] = globalElementsAllData[global_element_i].y;
        //         global_element_coords[3*global_element_i+2] = globalElementsAllData[global_element_i].z;

        //     }
            
        // }
        break;
    }
    case ElementType::HEX:
    {
        std::vector<HexElementWithFaces> localElementsAllData;
        GetElementsWithFacesCentroids<HexElementWithFaces>(file_path, localElementsAllData, ElementType::HEX, MPI_COMM_WORLD);
        SetMortonEncoding(localElementsAllData,ElementType::HEX,MPI_COMM_WORLD);
        std::vector<HexElementWithFaces> localElementsAllDataSorted(localElementsAllData.size());
        MPI_Barrier(MPI_COMM_WORLD);

        par::sampleSort<HexElementWithFaces>(localElementsAllData,localElementsAllDataSorted,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        print_log("[", taskid, "]: ", "global SFC sort done");


        local_element_count = localElementsAllDataSorted.size();
        local_elements.resize(local_element_count);
        MPI_Allgather(&local_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,MPI_COMM_WORLD);

        MPI_Allreduce(&local_element_count,&global_element_count,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);

        omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],numtasks);
        uint64_t global_idx_start = proc_element_counts_scanned[taskid];
        // print_log("[", taskid, "]: proc_element_counts_scanned ", VectorToString(proc_element_counts_scanned));
        for (size_t local_elem_i = 0; local_elem_i < local_element_count; local_elem_i++)
        {
            localElementsAllDataSorted[local_elem_i].global_idx = global_idx_start + local_elem_i;
            local_elements[local_elem_i].global_idx = global_idx_start + local_elem_i;

            local_elements[local_elem_i].element_tag = localElementsAllDataSorted[local_elem_i].element_tag;
            local_elements[local_elem_i].x = localElementsAllDataSorted[local_elem_i].x;
            local_elements[local_elem_i].y = localElementsAllDataSorted[local_elem_i].y;
            local_elements[local_elem_i].z = localElementsAllDataSorted[local_elem_i].z;

        }
        // print_log("[", taskid, "]:", "localElementsAllDataSorted = ", VectorToString(localElementsAllDataSorted));
        
    



        ResolveLocalElementConnectivity<HexElementWithFaces>(localElementsAllDataSorted,ElementType::HEX,local_connected_element_pairs,local_unconnected_elements_faces);
        print_log("[", taskid, "]:", "local_elements = ", VectorToString(local_elements));

        // print_log("[", taskid, "]:", "local_connected_element_pairs = ", VectorToString(local_connected_element_pairs));
        // print_log("[", taskid, "]:", "local_unconnected_element_face = ", VectorToString(local_unconnected_elements_faces));    



        ResolveBoundaryElementConnectivity(local_unconnected_elements_faces,proc_element_counts, boundary_connected_element_pairs,MPI_COMM_WORLD);
        // std::vector<int> testvec = {456,34547,56,67,8967,956,85,6867,93,6,3452,3524};
        // omp_par::merge_sort(&testvec[0],&testvec[testvec.size()]);
        // print_log(VectorToString(testvec));
        
        
        // auto displacements = GetDisplacementsFromCounts(proc_element_counts);

        // std::vector<int> proc_element_counts_(proc_element_counts.begin(), proc_element_counts.end());
        // std::vector<HexElementWithFaces> globalElementsAllData;;
        // if (!taskid)
        // {
        //     globalElementsAllData.resize(global_element_count);
        //     global_element_coords.resize(global_element_count*3);
        // }
        
        // MPI_Gatherv(localElementsAllDataSorted.data(), local_element_count, par::Mpi_datatype<HexElementWithFaces>::value(), 
        //     globalElementsAllData.data(), proc_element_counts_.data(), displacements.data(), par::Mpi_datatype<HexElementWithFaces>::value(), 0, MPI_COMM_WORLD);

        // if (!taskid)
        // {
        //     for (size_t global_element_i = 0; global_element_i < global_element_count; global_element_i++)
        //     {
        //         global_element_coords[3*global_element_i] = globalElementsAllData[global_element_i].x;
        //         global_element_coords[3*global_element_i+1] = globalElementsAllData[global_element_i].y;
        //         global_element_coords[3*global_element_i+2] = globalElementsAllData[global_element_i].z;

        //     }
            
        // }
        

        break;
    }

    default:
        break;
    }

    std::vector<ElementWithTag> ghost_elements;
    std::vector<int> ghost_element_counts(numtasks);
    ExtractGhostElements(boundary_connected_element_pairs,proc_element_counts,proc_element_counts_scanned,ghost_elements,ghost_element_counts,MPI_COMM_WORLD);

    DistGraph dist_graph(local_elements,ghost_elements,local_connected_element_pairs,boundary_connected_element_pairs,proc_element_counts,
                            proc_element_counts_scanned,ghost_element_counts,MPI_COMM_WORLD);

    print_log("[", taskid, "]:\n", dist_graph.GraphToString());

    // std::vector<ElementWithTag> ghost_elements_dups(boundary_connected_element_pairs.size());

    // //TODO: can be parallelized
    // for (size_t boudary_connect_i = 0; boudary_connect_i < boundary_connected_element_pairs.size(); boudary_connect_i++)
    // {
    //     ghost_elements_dups[boudary_connect_i].element_tag = boundary_connected_element_pairs[boudary_connect_i].second.element_tag;
    //     ghost_elements_dups[boudary_connect_i].global_idx = boundary_connected_element_pairs[boudary_connect_i].second.global_idx;

    // }

    // omp_par::merge_sort(&ghost_elements_dups[0], &ghost_elements_dups[ghost_elements_dups.size()]);

    // // print_log("[", taskid, "]:", "ghost_elements_dups sorted = ", VectorToString(ghost_elements_dups));    
    
    // std::vector<ElementWithTag> ghost_elements;
    // ghost_elements.push_back(ghost_elements_dups[0]);
    // for (size_t ghost_element_dup_i = 1; ghost_element_dup_i < ghost_elements_dups.size(); ghost_element_dup_i++)
    // {
    //     if (ghost_elements_dups[ghost_element_dup_i].global_idx != ghost_elements_dups[ghost_element_dup_i-1].global_idx)
    //     {
    //         ghost_elements.push_back(ghost_elements_dups[ghost_element_dup_i]);
    //     }
        
    // }

    // print_log("[", taskid, "]:", "ghost_elements = ", VectorToString(ghost_elements));    
    // std::vector<int> ghost_element_counts(numtasks);
    // {
    //     uint64_t current_proc = 0;
    //     for (size_t ghost_element_i = 0; ghost_element_i < ghost_elements.size(); ghost_element_i++)
    //     {
    //         if (ghost_elements[ghost_element_i].global_idx >= proc_element_counts_scanned[current_proc] &&
    //             ghost_elements[ghost_element_i].global_idx < (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc]))
    //         {
    //             ghost_element_counts[current_proc]++;
    //         } else
    //         {
    //             while (1)
    //             {
    //                 current_proc++;
    //                 if (ghost_elements[ghost_element_i].global_idx >= proc_element_counts_scanned[current_proc] &&
    //                     ghost_elements[ghost_element_i].global_idx < (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc]))
    //                 {
    //                     ghost_element_counts[current_proc]++;
    //                     break;
    //                 }                    
    //             }
                
    //         }   
    //     }       

    // }
    // print_log("[", taskid, "]:", "ghost_element_counts = ", VectorToString(ghost_element_counts));    



    

    // if (!taskid)
    // {
    //     global_element_partition_labels_morton.resize(global_element_count);
    //     uint64_t sum=0;
    //     for (size_t part_i = 0; part_i < numtasks; part_i++)
    //     {
    //         uint64_t start = sum;
    //         uint64_t end = sum+proc_element_counts[part_i];
    //         for (size_t element_i = start; element_i < end; element_i++)
    //         {
    //             global_element_partition_labels_morton[element_i] = part_i;
    //         }
    //         sum+=proc_element_counts[part_i];
            
    //     }
    //     // print_log(VectorToString(global_element_coords));

    //     // print_log(VectorToString(global_element_partition_labels_morton));
    //     PointsWithPartitionsToVtk(global_element_coords,global_element_partition_labels_morton,global_element_count,"out-sfc.vtk");

        
    // }
    
    // std::vector<TetElementWithFaces> localElementsAllData(8);
    // print_log(VectorToString(localElementsAllData));
    // std::vector<double> elem_coords;
    // std::vector<size_t> elem_tags;

    // Graph element_connectivity_graph =  GmshGetElementGraph(file_path, elem_coords, elem_tags);
    // uint64_t element_count = element_connectivity_graph.GetSize();

    // // print_log(VectorToString(elem_coords));
    // std::vector<uint64_t> morton_order_indices = SortMorton(elem_coords, element_count);
    // std::vector<uint64_t> SFC_partition_labels(element_count);
    // AssignPartitionLabelsInOrder(morton_order_indices, element_count,partition_count,SFC_partition_labels);
    // PointsWithPartitionsToVtk(elem_coords,SFC_partition_labels,element_count,"out-sfc.vtk");

    // std::vector<uint64_t> BFS_seeds(partition_count);

    // GetSamplesFromOrdered(morton_order_indices, elem_tags,partition_count, BFS_seeds);
    // element_connectivity_graph.InitMultiBFS(BFS_seeds, partition_count);
    // auto start = std::chrono::high_resolution_clock::now();
    // element_connectivity_graph.RunMultiBFSToStable();
    // auto end = std::chrono::high_resolution_clock::now();
    // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // print_log("Multi BFS time:\t", duration.count(), " ms");

    // auto BFS_partition_labels = element_connectivity_graph.GetMultiBFSLabels();


    // // std::vector<double> BFS_partition_labels_(BFS_partition_labels.begin(), BFS_partition_labels.end());
    // PointsWithPartitionsToVtk(elem_coords,BFS_partition_labels,element_count,"out-bfs.vtk");
    // // PointsToVtk(elem_integer_coords_,element_count, "out-integer.vtk");
    // auto xadj = element_connectivity_graph.GetCSR_xadj();
    // auto adjncy = element_connectivity_graph.GetCSR_adjncy();
    // print_log("csr done");
    // auto metis_partition_labels = GetMETISPartitions(xadj, adjncy, (int32_t)element_count, (int32_t)partition_count);
    // PointsWithPartitionsToVtk(elem_coords,metis_partition_labels,element_count,"out-metis.vtk");
    // }
    
    MPI_Finalize();
    return 0;
}
