#include "util.hpp"
#include "mesh-util.hpp"
#include "vtk-util.hpp"
#include "graph.hpp"
#include "sfc.hpp"
#include "metis-util.hpp"
#include <string>
#include <stdexcept>

#include <mpi.h>
#include <algorithm>
#include <chrono>
#include <omp.h>

#include <parUtils.h>
#include <dtypes.h>

int main(int argc, char *argv[])
{
    int numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    uint64_t partition_count=3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(hostname, &len);
    // if (true || taskid == 0)
    // {
    // omp_set_num_threads(8);
    // const std::string file_path("/home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh");
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh");
    const std::string file_path("/home/budvin/research/Partitioning/mesh_generator/hex-box-25x25x3.msh");

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


    switch (elementType)
    {
    case ElementType::TET:
    {
        std::vector<TetElementWithFaces> localElementsAllData;
        GetElementsWithFacesCentroids<TetElementWithFaces>(file_path, localElementsAllData, ElementType::TET, MPI_COMM_WORLD);
        SetMortonEncoding(localElementsAllData,ElementType::TET,MPI_COMM_WORLD);
        std::vector<TetElementWithFaces> localElementsAllDataSorted(localElementsAllData.size());
        MPI_Barrier(MPI_COMM_WORLD);

        par::sampleSort<TetElementWithFaces>(localElementsAllData,localElementsAllDataSorted,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        local_element_count = localElementsAllDataSorted.size();
        MPI_Gather(&local_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,0,MPI_COMM_WORLD);
        MPI_Reduce(&local_element_count,&global_element_count,1,MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);
        print_log("[",taskid,"] ", VectorToString(proc_element_counts));
        print_log("[",taskid,"] ", global_element_count);
        auto displacements = GetDisplacementsFromCounts(proc_element_counts);
        std::vector<int> proc_element_counts_(proc_element_counts.begin(), proc_element_counts.end());
        
        std::vector<TetElementWithFaces> globalElementsAllData;
        if (!taskid)
        {
            globalElementsAllData.resize(global_element_count);
            global_element_coords.resize(global_element_count*3);
        }
        MPI_Gatherv(localElementsAllDataSorted.data(), local_element_count, par::Mpi_datatype<TetElementWithFaces>::value(), 
            globalElementsAllData.data(), proc_element_counts_.data(), displacements.data(), par::Mpi_datatype<TetElementWithFaces>::value(), 0, MPI_COMM_WORLD);
        
        if (!taskid)
        {
            for (size_t global_element_i = 0; global_element_i < global_element_count; global_element_i++)
            {
                global_element_coords[3*global_element_i] = globalElementsAllData[global_element_i].x;
                global_element_coords[3*global_element_i+1] = globalElementsAllData[global_element_i].y;
                global_element_coords[3*global_element_i+2] = globalElementsAllData[global_element_i].z;

            }
            
        }
        break;
    }
    case ElementType::HEX:
    {
        std::vector<HexElementWithFaces> localElementsAllData;
        GetElementsWithFacesCentroids<HexElementWithFaces>(file_path, localElementsAllData, ElementType::HEX, MPI_COMM_WORLD);
        SetMortonEncoding(localElementsAllData,ElementType::TET,MPI_COMM_WORLD);
        std::vector<HexElementWithFaces> localElementsAllDataSorted(localElementsAllData.size());
        MPI_Barrier(MPI_COMM_WORLD);

        par::sampleSort<HexElementWithFaces>(localElementsAllData,localElementsAllDataSorted,MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        local_element_count = localElementsAllDataSorted.size();
        MPI_Gather(&local_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,0,MPI_COMM_WORLD);

        MPI_Reduce(&local_element_count,&global_element_count,1,MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);
        print_log("[",taskid,"] ", VectorToString(proc_element_counts));
        print_log("[",taskid,"] ", global_element_count);
        auto displacements = GetDisplacementsFromCounts(proc_element_counts);
        print_log("here...");

        std::vector<int> proc_element_counts_(proc_element_counts.begin(), proc_element_counts.end());
        std::vector<HexElementWithFaces> globalElementsAllData;;
        if (!taskid)
        {
            globalElementsAllData.resize(global_element_count);
            global_element_coords.resize(global_element_count*3);
        }
        
        MPI_Gatherv(localElementsAllDataSorted.data(), local_element_count, par::Mpi_datatype<HexElementWithFaces>::value(), 
            globalElementsAllData.data(), proc_element_counts_.data(), displacements.data(), par::Mpi_datatype<HexElementWithFaces>::value(), 0, MPI_COMM_WORLD);

        if (!taskid)
        {
            for (size_t global_element_i = 0; global_element_i < global_element_count; global_element_i++)
            {
                global_element_coords[3*global_element_i] = globalElementsAllData[global_element_i].x;
                global_element_coords[3*global_element_i+1] = globalElementsAllData[global_element_i].y;
                global_element_coords[3*global_element_i+2] = globalElementsAllData[global_element_i].z;

            }
            
        }
        

        break;
    }

    default:
        break;
    }
    if (!taskid)
    {
        global_element_partition_labels_morton.resize(global_element_count);
        uint64_t sum=0;
        for (size_t part_i = 0; part_i < numtasks; part_i++)
        {
            uint64_t start = sum;
            uint64_t end = sum+proc_element_counts[part_i];
            for (size_t element_i = start; element_i < end; element_i++)
            {
                global_element_partition_labels_morton[element_i] = part_i;
            }
            sum+=proc_element_counts[part_i];
            
        }
        // print_log(VectorToString(global_element_coords));

        // print_log(VectorToString(global_element_partition_labels_morton));
        PointsWithPartitionsToVtk(global_element_coords,global_element_partition_labels_morton,global_element_count,"out-sfc.vtk");

        
    }
    
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
