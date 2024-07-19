#include "util.hpp"
#include "mesh-util.hpp"
// #include "graph.hpp"
#include "sfc.hpp"
#include "metis-util.hpp"
#include "dist-graph.hpp"
#include "linalg.hpp"
#include <string>
#include <stdexcept>

#include <mpi.h>
#include <algorithm>
#include <chrono>
#include <omp.h>

#include <parUtils.h>
#include <dtypes.h>

#include "ompUtils.h"

#ifdef ENABLE_VTK_FEATURES
#include "vtk-util.hpp"
#endif

struct SFCStatus
{
    int return_code;
    int time_us;
};

template <class T>
SFCStatus ReadAndDistributeSFC(std::string mesh_file_path, ElementType element_type, std::vector<T>& elements_out, MPI_Comm comm);

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
    // const std::string file_path("/home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh");
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh");
    // const std::string file_path("/home/budvin/research/Partitioning/mesh_generator/hex-box-3x3x3.msh");

    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh");   //octopus
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/196209_sf_hexa.mesh_73346_289961.obj.mesh");     //large tet
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_tet/75651_sf_hexa.mesh_78608_298692.obj.mesh");  //largest tet
    // const std::string file_path("/home/budvin/research/Partitioning/Meshes/10k_hex/75651_sf_hexa.mesh");  //largest hex

    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <mesh file path> <file index> <run index> <metrics out file path> <-viz or -no-viz>" << std::endl;
        return 1; // indicating an error
    }
    const std::string mesh_file_path = argv[1];


    int file_idx = std::stoi(argv[2]);
    int run_idx = std::stoi(argv[3]);


    const std::string metrics_out_file_path = argv[4];
    const std::string viz_flag_str = argv[5];
    bool viz_flag;

    if (viz_flag_str == "-viz")
    {
        viz_flag = true;
    }else if (viz_flag_str == "-no-viz")
    {
        viz_flag = false;
    }else
    {
        std::cerr << "Invalid flag: " << viz_flag_str << std::endl;
        std::cerr << "Usage: " << argv[0] << " <mesh file path> <file index> <run index> <metrics out file path> <-viz or -no-viz>" << std::endl;
        return 1;
    }
    
    
    
    


    if(!taskid){
        print_log("running on", numtasks, "MPI processs");
        print_log("partitioning: ", mesh_file_path);
    }

    ElementType elementType = GetElementType(mesh_file_path, MPI_COMM_WORLD);
    // print_log("element type: ", elementType);
    uint64_t local_element_count;
    uint64_t global_element_count;

    std::vector<double> global_element_coords;
    std::vector<uint64_t> global_element_partition_labels_morton;

    std::vector<uint64_t> proc_element_counts(numtasks);
    std::vector<uint64_t> proc_element_counts_scanned(numtasks);


    std::vector<TetElementWithFacesNodes> localElementsAllData_Tet;     // contains all information about face tags, node tags
    std::vector<HexElementWithFacesNodes> localElementsAllData_Hex;

    std::vector<ElementWithCoord> local_elements;       // contains {element_tag, global_idx, [x,y,z]}
    std::vector<std::pair<ElementWithTag, ElementWithTag>> local_connected_element_pairs;
    std::vector<std::pair<ElementWithTag, ElementWithTag>> boundary_connected_element_pairs;

    std::vector<ElementWithFace> local_unconnected_elements_faces;


    SFCStatus sfc_status;
    
    switch (elementType)
    {
    case ElementType::TET:
    {
        sfc_status = ReadAndDistributeSFC(mesh_file_path, ElementType::TET, localElementsAllData_Tet, MPI_COMM_WORLD);
        local_element_count = localElementsAllData_Tet.size();
        break;
    }
    case ElementType::HEX:
    {
        sfc_status = ReadAndDistributeSFC(mesh_file_path, ElementType::HEX, localElementsAllData_Hex, MPI_COMM_WORLD);
        local_element_count = localElementsAllData_Hex.size();
        break;
    }
    
    default:
    {
        throw std::runtime_error("Unknown element type");
        break;
    }
        
    }
    
   
    
    MPI_Allgather(&local_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,MPI_COMM_WORLD);

    MPI_Allreduce(&local_element_count,&global_element_count,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);

    omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],numtasks);

    local_elements.resize(local_element_count);
    uint64_t global_idx_start = proc_element_counts_scanned[taskid];


    MPI_Barrier(MPI_COMM_WORLD);
    auto graph_setup_start = std::chrono::high_resolution_clock::now();

    switch (elementType)
    {
    case ElementType::TET:
    {
        for (size_t local_elem_i = 0; local_elem_i < local_element_count; local_elem_i++)
        {
            localElementsAllData_Tet[local_elem_i].global_idx = global_idx_start + local_elem_i;
            local_elements[local_elem_i].global_idx = global_idx_start + local_elem_i;

            local_elements[local_elem_i].element_tag = localElementsAllData_Tet[local_elem_i].element_tag;
            local_elements[local_elem_i].x = localElementsAllData_Tet[local_elem_i].x;
            local_elements[local_elem_i].y = localElementsAllData_Tet[local_elem_i].y;
            local_elements[local_elem_i].z = localElementsAllData_Tet[local_elem_i].z;

        }
        ResolveLocalElementConnectivity(localElementsAllData_Tet,ElementType::TET,local_connected_element_pairs,local_unconnected_elements_faces);
        break;
    }
    case ElementType::HEX:
    {

        for (size_t local_elem_i = 0; local_elem_i < local_element_count; local_elem_i++)
        {
            localElementsAllData_Hex[local_elem_i].global_idx = global_idx_start + local_elem_i;
            local_elements[local_elem_i].global_idx = global_idx_start + local_elem_i;

            local_elements[local_elem_i].element_tag = localElementsAllData_Hex[local_elem_i].element_tag;
            local_elements[local_elem_i].x = localElementsAllData_Hex[local_elem_i].x;
            local_elements[local_elem_i].y = localElementsAllData_Hex[local_elem_i].y;
            local_elements[local_elem_i].z = localElementsAllData_Hex[local_elem_i].z;

        }
        ResolveLocalElementConnectivity(localElementsAllData_Hex,ElementType::HEX,local_connected_element_pairs,local_unconnected_elements_faces);
        break;
    }
    
    default:
    {
        throw std::runtime_error("Unknown element type");
        break;
    }
        
    }


    // print_log("[", taskid, "]:", "local_elements = ", VectorToString(local_elements));


    ResolveBoundaryElementConnectivity(local_unconnected_elements_faces,proc_element_counts, boundary_connected_element_pairs,MPI_COMM_WORLD);



    std::vector<ElementWithTag> ghost_elements;
    std::vector<int> ghost_element_counts(numtasks);
    ExtractGhostElements(boundary_connected_element_pairs,proc_element_counts,proc_element_counts_scanned,ghost_elements,ghost_element_counts,MPI_COMM_WORLD);

    DistGraph dist_graph(local_elements,ghost_elements,local_connected_element_pairs,boundary_connected_element_pairs,proc_element_counts,
                            proc_element_counts_scanned,ghost_element_counts,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    auto graph_setup_end = std::chrono::high_resolution_clock::now();

    if(!taskid) print_log("graph setup done");
    auto graph_setup_duration = std::chrono::duration_cast<std::chrono::microseconds>(graph_setup_end - graph_setup_start);
    if(!taskid) print_log("graph setup time:", graph_setup_duration.count(), "us");



    if(!taskid) print_log("starting BFS partitioning");
    std::vector<uint16_t> local_bfs_partition_labels(local_element_count);      // TODO: make the label type consistent with bfs_label_t or int32
    auto bfs_status = dist_graph.PartitionBFS(local_bfs_partition_labels);
    if(!taskid) print_log("BFS partitioning done");

    if(!taskid) print_log("starting parmetis");
    std::vector<uint16_t> local_parmetis_partition_labels(local_element_count);
    auto parmetis_status = dist_graph.PartitionParmetis(local_parmetis_partition_labels);
    if(!taskid) print_log("parmetis done");

    if(!taskid) print_log("starting ptscotch");
    std::vector<uint16_t> local_ptscotch_partition_labels(local_element_count);
    auto ptscotch_status = dist_graph.PartitionPtScotch(local_ptscotch_partition_labels);
    if(!taskid) print_log("ptscotch done");


    std::vector<uint16_t> local_sfc_partition_labels(local_element_count,taskid);



    // collecting partitioning metrics

    std::vector<uint32_t> global_bfs_partition_sizes;
    std::vector<uint32_t> global_bfs_partition_boundaries;
    dist_graph.GetPartitionMetrics(local_bfs_partition_labels,global_bfs_partition_sizes, global_bfs_partition_boundaries);


    std::vector<uint32_t> global_sfc_partition_sizes;
    std::vector<uint32_t> global_sfc_partition_boundaries;
    dist_graph.GetPartitionMetrics(local_sfc_partition_labels,global_sfc_partition_sizes, global_sfc_partition_boundaries);

    std::vector<uint32_t> global_parmetis_partition_sizes;
    std::vector<uint32_t> global_parmetis_partition_boundaries;
    dist_graph.GetPartitionMetrics(local_parmetis_partition_labels,global_parmetis_partition_sizes, global_parmetis_partition_boundaries);

    std::vector<uint32_t> global_ptscotch_partition_sizes;
    std::vector<uint32_t> global_ptscotch_partition_boundaries;
    dist_graph.GetPartitionMetrics(local_ptscotch_partition_labels,global_ptscotch_partition_sizes, global_ptscotch_partition_boundaries);



    #ifdef ENABLE_VTK_FEATURES
    if (viz_flag)
    {

        // collecting to root process for visualization
        std::vector<ElementWithCoord> global_all_elements;
        std::vector<uint16_t> global_all_elements_bfs_partition_labels;
        std::vector<uint16_t> global_all_elements_sfc_partition_labels;
        std::vector<uint16_t> global_all_elements_parmetis_partition_labels;
        std::vector<uint16_t> global_all_elements_ptscotch_partition_labels;


        if (! taskid)
        {
            global_all_elements.resize(global_element_count);
            global_all_elements_bfs_partition_labels.resize(global_element_count);
            global_all_elements_sfc_partition_labels.resize(global_element_count);
            global_all_elements_parmetis_partition_labels.resize(global_element_count);
            global_all_elements_ptscotch_partition_labels.resize(global_element_count);


        }
        
        

        std::vector<int> proc_element_counts_(proc_element_counts.begin(), proc_element_counts.end());
        std::vector<int> proc_element_counts_scanned_(proc_element_counts_scanned.begin(), proc_element_counts_scanned.end());


        MPI_Gatherv(local_elements.data(),local_element_count,par::Mpi_datatype<ElementWithCoord>::value(),global_all_elements.data(),
                    proc_element_counts_.data(),proc_element_counts_scanned_.data(),par::Mpi_datatype<ElementWithCoord>::value(), 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Gatherv(local_bfs_partition_labels.data(),local_element_count,MPI_UINT16_T,global_all_elements_bfs_partition_labels.data(),
                    proc_element_counts_.data(),proc_element_counts_scanned_.data(),MPI_UINT16_T, 0, MPI_COMM_WORLD);

        
        
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gatherv(local_sfc_partition_labels.data(),local_element_count,MPI_UINT16_T,global_all_elements_sfc_partition_labels.data(),
                    proc_element_counts_.data(),proc_element_counts_scanned_.data(),MPI_UINT16_T, 0, MPI_COMM_WORLD);


        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Gatherv(local_parmetis_partition_labels.data(),local_element_count,MPI_UINT16_T,global_all_elements_parmetis_partition_labels.data(),
                    proc_element_counts_.data(),proc_element_counts_scanned_.data(),MPI_UINT16_T, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Gatherv(local_ptscotch_partition_labels.data(),local_element_count,MPI_UINT16_T,global_all_elements_ptscotch_partition_labels.data(),
                    proc_element_counts_.data(),proc_element_counts_scanned_.data(),MPI_UINT16_T, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        if (!taskid)
        {
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_bfs_partition_labels, global_element_count, "out-bfs.vtk");
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_sfc_partition_labels, global_element_count, "out-sfc.vtk");
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_parmetis_partition_labels, global_element_count, "out-parmetis.vtk");
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_ptscotch_partition_labels, global_element_count, "out-ptscotch.vtk");
        }
        
    }

    #endif
    
    DistributionStatus bfs_distribution_status;
    DistributionStatus parmetis_distribution_status;
    DistributionStatus ptscotch_distribution_status;


    SpMVStatus sfc_spmv_status;
    SpMVStatus bfs_spmv_status;
    SpMVStatus parmetis_spmv_status;
    SpMVStatus ptscotch_spmv_status;




    switch (elementType)
    {
    case ElementType::TET:
    {
        if(!taskid) print_log("testing SFC partitioning");
        sfc_spmv_status = TestSpMV(localElementsAllData_Tet, ElementType::TET, viz_flag ,MPI_COMM_WORLD);

        if(!taskid) print_log("testing BFS partitioning");
        std::vector<TetElementWithFacesNodes> localElementsAllData_Tet_2;
        bfs_distribution_status = Redistribute<TetElementWithFacesNodes>(localElementsAllData_Tet,local_bfs_partition_labels,localElementsAllData_Tet_2, MPI_COMM_WORLD);
        bfs_spmv_status = TestSpMV(localElementsAllData_Tet_2, ElementType::TET,viz_flag,MPI_COMM_WORLD);


        if(!taskid) print_log("testing parMETIS partitioning");
        localElementsAllData_Tet_2.clear();
        parmetis_distribution_status = Redistribute<TetElementWithFacesNodes>(localElementsAllData_Tet,local_parmetis_partition_labels,localElementsAllData_Tet_2, MPI_COMM_WORLD);
        parmetis_spmv_status = TestSpMV(localElementsAllData_Tet_2, ElementType::TET,viz_flag,MPI_COMM_WORLD);   

        if(!taskid) print_log("testing ptscotch partitioning");
        localElementsAllData_Tet_2.clear();
        ptscotch_distribution_status = Redistribute<TetElementWithFacesNodes>(localElementsAllData_Tet,local_ptscotch_partition_labels,localElementsAllData_Tet_2, MPI_COMM_WORLD);
        ptscotch_spmv_status = TestSpMV(localElementsAllData_Tet_2, ElementType::TET,viz_flag,MPI_COMM_WORLD);   


        break;
    }
    case ElementType::HEX:
    {
        if(!taskid) print_log("testing SFC partitioning");
        sfc_spmv_status = TestSpMV(localElementsAllData_Hex, ElementType::HEX,viz_flag,MPI_COMM_WORLD);

        if(!taskid) print_log("testing BFS partitioning");
        std::vector<HexElementWithFacesNodes> localElementsAllData_Hex_2;
        bfs_distribution_status = Redistribute<HexElementWithFacesNodes>(localElementsAllData_Hex,local_bfs_partition_labels,localElementsAllData_Hex_2, MPI_COMM_WORLD);
        bfs_spmv_status = TestSpMV(localElementsAllData_Hex_2, ElementType::HEX,viz_flag,MPI_COMM_WORLD);


        if(!taskid) print_log("testing parMETIS partitioning");
        localElementsAllData_Hex_2.clear();
        parmetis_distribution_status = Redistribute<HexElementWithFacesNodes>(localElementsAllData_Hex,local_parmetis_partition_labels,localElementsAllData_Hex_2, MPI_COMM_WORLD);
        parmetis_spmv_status = TestSpMV(localElementsAllData_Hex_2, ElementType::HEX,viz_flag,MPI_COMM_WORLD);

        if(!taskid) print_log("testing ptscotch partitioning");
        localElementsAllData_Hex_2.clear();
        ptscotch_distribution_status = Redistribute<HexElementWithFacesNodes>(localElementsAllData_Hex,local_ptscotch_partition_labels,localElementsAllData_Hex_2, MPI_COMM_WORLD);
        ptscotch_spmv_status = TestSpMV(localElementsAllData_Hex_2, ElementType::HEX,viz_flag,MPI_COMM_WORLD);
        break;
    }
    
    default:
    {
        throw std::runtime_error("Unknown element type");
        break;
    }
        
    }
    if (! taskid)
    {   

        ExportMetricsToJson(mesh_file_path, file_idx, run_idx, numtasks, global_element_count,
                                graph_setup_duration.count(),
                                global_sfc_partition_sizes, global_sfc_partition_boundaries, sfc_status.time_us, sfc_spmv_status.mat_assembly_time_us, sfc_spmv_status.matvec_time_us,
                                global_bfs_partition_sizes,global_bfs_partition_boundaries, bfs_status.time_us, bfs_distribution_status.time_us, bfs_spmv_status.mat_assembly_time_us, bfs_spmv_status.matvec_time_us,
                                global_parmetis_partition_sizes,global_parmetis_partition_boundaries, parmetis_status.time_us, parmetis_distribution_status.time_us, parmetis_spmv_status.mat_assembly_time_us, parmetis_spmv_status.matvec_time_us,
                                global_ptscotch_partition_sizes,global_ptscotch_partition_boundaries, ptscotch_status.time_us, ptscotch_distribution_status.time_us, ptscotch_spmv_status.mat_assembly_time_us, ptscotch_spmv_status.matvec_time_us,
                                metrics_out_file_path);

    }
   
    MPI_Finalize();
    return 0;
}

template <class T>
SFCStatus ReadAndDistributeSFC(std::string mesh_file_path, ElementType element_type, std::vector<T>& elements_out, MPI_Comm comm){

    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    std::vector<T> initial_elements;        // before SFC
    GetInitialElementsDistribution<T>(mesh_file_path, initial_elements, element_type, MPI_COMM_WORLD);
    SetMortonEncoding(initial_elements,element_type,MPI_COMM_WORLD);
    elements_out.resize(initial_elements.size());
    if (! my_rank)
    {
        print_log("starting sfc sort");
    }
    MPI_Barrier(MPI_COMM_WORLD);
    auto start = std::chrono::high_resolution_clock::now();
    par::sampleSort<T>(initial_elements,elements_out,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (! my_rank)
    {
        print_log("global sfc sort done");
        print_log("SFC sort time:\t", duration.count(), " us");
    }

    SFCStatus status;
    status.return_code = 0;
    status.time_us = duration.count();

    return status;
}
