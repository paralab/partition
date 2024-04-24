#include "util.hpp"
#include "mesh-util.hpp"
#include "vtk-util.hpp"
#include "graph.hpp"
#include "sfc.hpp"
#include "metis-util.hpp"
#include <string>

#include <mpi.h>
#include <algorithm>
#include <chrono>
#include <omp.h>

int main(int argc, char *argv[])
{
    int numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    uint64_t partition_count=10;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(hostname, &len);
    if (taskid == 0)
    {
        // omp_set_num_threads(8);
        // const std::string file_path("/home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh");
        // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh");
        const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh");   //octopus
        // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/196209_sf_hexa.mesh_73346_289961.obj.mesh");     //large tet
        // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/75651_sf_hexa.mesh_78608_298692.obj.mesh");  //largest tet
        // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_hex/75651_sf_hexa.mesh");  //largest hex

        std::vector<double> elem_coords;
        std::vector<size_t> elem_tags;

        Graph element_connectivity_graph =  GmshGetElementGraph(file_path, elem_coords, elem_tags);
        uint64_t element_count = element_connectivity_graph.GetSize();

        // print_log(VectorToString(elem_coords));
        std::vector<uint64_t> morton_order_indices = SortMorton(elem_coords, element_count);
        std::vector<uint64_t> SFC_partition_labels(element_count);
        AssignPartitionLabelsInOrder(morton_order_indices, element_count,partition_count,SFC_partition_labels);
        PointsWithPartitionsToVtk(elem_coords,SFC_partition_labels,element_count,"out-sfc.vtk");

        std::vector<uint64_t> BFS_seeds(partition_count);

        GetSamplesFromOrdered(morton_order_indices, elem_tags,partition_count, BFS_seeds);
        element_connectivity_graph.InitMultiBFS(BFS_seeds, partition_count);
        auto start = std::chrono::high_resolution_clock::now();
        element_connectivity_graph.RunMultiBFSToStable();
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        print_log("Multi BFS time:\t", duration.count(), " ms");

        auto BFS_partition_labels = element_connectivity_graph.GetMultiBFSLabels();


        // std::vector<double> BFS_partition_labels_(BFS_partition_labels.begin(), BFS_partition_labels.end());
        PointsWithPartitionsToVtk(elem_coords,BFS_partition_labels,element_count,"out-bfs.vtk");
        // PointsToVtk(elem_integer_coords_,element_count, "out-integer.vtk");
        auto xadj = element_connectivity_graph.GetCSR_xadj();
        auto adjncy = element_connectivity_graph.GetCSR_adjncy();
        print_log("csr done");
        auto metis_partition_labels = GetMETISPartitions(xadj, adjncy, (int32_t)element_count, (int32_t)partition_count);
        PointsWithPartitionsToVtk(elem_coords,metis_partition_labels,element_count,"out-metis.vtk");
    }
    
    MPI_Finalize();
    return 0;
}
