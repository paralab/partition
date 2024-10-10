#include "util.hpp"
#include "mesh-util.hpp"
// #include "graph.hpp"
#include "sfc.hpp"
#include "metis-util.hpp"
#include "dist-graph.hpp"
#include "linalg.hpp"
#include <string>
#include <stdexcept>
#include <fstream>

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

#define OCTREE_BOUNDARY UINT_MAX

typedef long int D_INT_L;
struct SFCStatus
{
    int return_code;
    int time_us;
};

// template <class T>
SFCStatus ReadAndDistributeSFC(std::string mesh_file_path, std::vector<OctreeElementWithNeigh>& elements_out, MPI_Comm comm);

void ResolveOctreeElementConnectivity(std::vector<OctreeElementWithNeigh> &local_elements,
                                std::vector<uint64_t> &proc_element_counts,
                                std::vector<uint64_t> &proc_element_counts_scanned,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &connected_element_pairs_out,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &boundary_connected_element_pairs_out,
                                MPI_Comm comm);

int main(int argc, char *argv[])
{
    int numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // uint64_t partition_count=3;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(hostname, &len);

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

    // ElementType elementType = GetElementType(mesh_file_path, MPI_COMM_WORLD);
    // print_log("element type: ", elementType);
    uint64_t local_element_count = 0;
    uint64_t global_element_count = 0;

    // std::vector<double> global_element_coords;
    // std::vector<uint64_t> global_element_partition_labels_morton;

    std::vector<uint64_t> proc_element_counts(numtasks);
    std::vector<uint64_t> proc_element_counts_scanned(numtasks);


    std::vector<OctreeElementWithNeigh> localElementsAllData;     // contains all information with neighbors

    std::vector<ElementWithCoord> local_elements;       // contains {element_tag, global_idx, [x,y,z]}
    std::vector<std::pair<ElementWithTag, ElementWithTag>> local_connected_element_pairs;
    std::vector<std::pair<ElementWithTag, ElementWithTag>> boundary_connected_element_pairs;

    std::vector<ElementWithFace> local_unconnected_elements_faces;


    SFCStatus sfc_status = ReadAndDistributeSFC(mesh_file_path, localElementsAllData, MPI_COMM_WORLD);

   
    local_element_count = localElementsAllData.size();
    
    MPI_Allgather(&local_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,MPI_COMM_WORLD);

    MPI_Allreduce(&local_element_count,&global_element_count,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);

    omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],numtasks);

    local_elements.resize(local_element_count);
    uint64_t global_idx_start = proc_element_counts_scanned[taskid];


    MPI_Barrier(MPI_COMM_WORLD);
    auto graph_setup_start = std::chrono::high_resolution_clock::now();


    for (size_t local_elem_i = 0; local_elem_i < local_element_count; local_elem_i++)
    {
        localElementsAllData[local_elem_i].global_idx = global_idx_start + local_elem_i;
        local_elements[local_elem_i].global_idx = global_idx_start + local_elem_i;

        local_elements[local_elem_i].element_tag = localElementsAllData[local_elem_i].element_tag;
        local_elements[local_elem_i].x = localElementsAllData[local_elem_i].x;
        local_elements[local_elem_i].y = localElementsAllData[local_elem_i].y;
        local_elements[local_elem_i].z = localElementsAllData[local_elem_i].z;

    }
    // if (taskid == 2)
    // {
    //     for (size_t i = 0; i < 10; i++)
    //     {
    //         print_log(localElementsAllData[i]);
    //     }
        
    // }
    
    

    ResolveOctreeElementConnectivity(localElementsAllData,proc_element_counts,proc_element_counts_scanned,local_connected_element_pairs,boundary_connected_element_pairs,MPI_COMM_WORLD);


    // // print_log("[", taskid, "]:", "local_elements = ", VectorToString(local_elements));

    // print_log("[", taskid, "]:","resolve done");



    std::vector<ElementWithTag> ghost_elements;
    std::vector<int> ghost_element_counts(numtasks);
    ExtractGhostElements(boundary_connected_element_pairs,proc_element_counts,proc_element_counts_scanned,ghost_elements,ghost_element_counts,MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    // print_log("[", taskid, "]:","ghost extract done");
    // MPI_Finalize();
    // return 0;

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

    //     // collecting to root process for visualization
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
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_bfs_partition_labels, global_element_count, "out-fastPart.vtk");
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_sfc_partition_labels, global_element_count, "out-sfc.vtk");
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_parmetis_partition_labels, global_element_count, "out-parmetis.vtk");
            ElementsWithPartitionsToVtk(global_all_elements, global_all_elements_ptscotch_partition_labels, global_element_count, "out-ptscotch.vtk");
        }
        
    }

    #endif

    
    DistributionStatus bfs_distribution_status;
    DistributionStatus parmetis_distribution_status;
    DistributionStatus ptscotch_distribution_status;



    if(!taskid) print_log("testing BFS partitioning");
    std::vector<OctreeElementWithNeigh> localElementsAllData_2;
    bfs_distribution_status = Redistribute<OctreeElementWithNeigh>(localElementsAllData,local_bfs_partition_labels,localElementsAllData_2, MPI_COMM_WORLD);


    if(!taskid) print_log("testing parMETIS partitioning");
    localElementsAllData_2.clear();
    parmetis_distribution_status = Redistribute<OctreeElementWithNeigh>(localElementsAllData,local_parmetis_partition_labels,localElementsAllData_2, MPI_COMM_WORLD);

    if(!taskid) print_log("testing ptscotch partitioning");
    localElementsAllData_2.clear();
    ptscotch_distribution_status = Redistribute<OctreeElementWithNeigh>(localElementsAllData,local_ptscotch_partition_labels,localElementsAllData_2, MPI_COMM_WORLD);




    if (! taskid)
    {   

        ExportMetricsToJson(mesh_file_path, file_idx, run_idx, numtasks, global_element_count,
                                graph_setup_duration.count(),
                                global_sfc_partition_sizes, global_sfc_partition_boundaries, sfc_status.time_us, 0, 0,
                                global_bfs_partition_sizes,global_bfs_partition_boundaries, bfs_status.time_us, bfs_distribution_status.time_us, 0, 0,
                                global_parmetis_partition_sizes,global_parmetis_partition_boundaries, parmetis_status.time_us, parmetis_distribution_status.time_us, 0, 0,
                                global_ptscotch_partition_sizes,global_ptscotch_partition_boundaries, ptscotch_status.time_us, ptscotch_distribution_status.time_us, 0, 0,
                                metrics_out_file_path);

    }
   
    MPI_Finalize();
    return 0;
}


void GetInitialElementsDistributionOct(const std::string &mesh_file_path, std::vector<OctreeElementWithNeigh> &local_elements_out, 
                                    MPI_Comm comm)
{
    struct oct_data
    {
        D_INT_L eid;
        D_INT_L coord[3];
        D_INT_L e2e[6];

        oct_data(){};

        oct_data(D_INT_L eid, D_INT_L coord[3], D_INT_L e2e[6])
        {
            this->eid=eid;

            this->coord[0] = coord[0];
            this->coord[1] = coord[1];
            this->coord[2] = coord[2];

            this->e2e[0]   = e2e[0];
            this->e2e[1]   = e2e[1];
            this->e2e[2]   = e2e[2];

            this->e2e[3]   = e2e[3];
            this->e2e[4]   = e2e[4];
            this->e2e[5]   = e2e[5];

        }

    };
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    int local_element_count;

    std::vector<int> proc_element_counts(procs_n);              // populated only at root   
    std::vector<int> proc_element_counts_scanned(procs_n);      // populated only at root
    std::vector<OctreeElementWithNeigh> all_elements;           // populated only at root



    if (!my_rank)
    {
        uint64_t raw_element_count = 0;
        {
            // Find the position of the underscore and the dot
            std::size_t underscorePos = mesh_file_path.rfind('_'); // Find last underscore
            std::size_t dotPos = mesh_file_path.rfind('.');        // Find last dot

            // Ensure that both underscore and dot exist in the filename and are in the correct order
            if (underscorePos != std::string::npos && dotPos != std::string::npos && underscorePos < dotPos) {
                // Extract the substring between the underscore and the dot
                std::string numberStr = mesh_file_path.substr(underscorePos + 1, dotPos - underscorePos - 1);

                // Convert the substring to uint64_t
                std::stringstream ss(numberStr);
                ss >> raw_element_count;

                // Output the extracted number
            } else {
                throw std::runtime_error("Incorrect file name format. Use *_n.octree format.");;
            }
        }
        print_log(raw_element_count, "\toctree elements");
        std::vector<oct_data> all_elements_raw(raw_element_count);
        all_elements.resize(raw_element_count);
        // all_elements.reserve(raw_element_count);

        std::ifstream inputFile(mesh_file_path, std::ios::binary);
        if (!inputFile) {
            throw std::runtime_error("Error opening mesh file.");
        }

        inputFile.read(reinterpret_cast<char*>(all_elements_raw.data()), raw_element_count * sizeof(oct_data));


        if (!inputFile) {
            throw std::runtime_error("Error reading from mesh file.");
        }

        inputFile.close();

        for (size_t i = 0; i < raw_element_count; i++)
        {

            all_elements[i].element_tag = all_elements_raw[i].eid;
            all_elements[i].x = (double)all_elements_raw[i].coord[0];
            all_elements[i].y = (double)all_elements_raw[i].coord[1];
            all_elements[i].z = (double)all_elements_raw[i].coord[2];

            all_elements[i].neigh[0] = all_elements_raw[i].e2e[0];
            all_elements[i].neigh[1] = all_elements_raw[i].e2e[1];
            all_elements[i].neigh[2] = all_elements_raw[i].e2e[2];
            all_elements[i].neigh[3] = all_elements_raw[i].e2e[3];
            all_elements[i].neigh[4] = all_elements_raw[i].e2e[4];
            all_elements[i].neigh[5] = all_elements_raw[i].e2e[5];

        }
        uint64_t total_element_count = raw_element_count;

        std::fill(proc_element_counts.begin(), proc_element_counts.end(), total_element_count/procs_n);

        // Distribute the remainder element counts
        int rem = total_element_count % procs_n;
        for (int proc_i = 0; proc_i < rem; proc_i++) {
            proc_element_counts[proc_i]++;
        }
        omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],procs_n);
        print_log("mesh reading done");
        
    }

    MPI_Scatter(proc_element_counts.data(), 1, MPI_INT, &local_element_count, 1, MPI_INT, 0, comm);
    local_elements_out.resize(local_element_count);
    MPI_Scatterv(all_elements.data(),
        proc_element_counts.data(), proc_element_counts_scanned.data(), par::Mpi_datatype<OctreeElementWithNeigh>::value(),
        local_elements_out.data(),local_element_count,par::Mpi_datatype<OctreeElementWithNeigh>::value(), 
        0, comm);
    

}                          

SFCStatus ReadAndDistributeSFC(std::string mesh_file_path, std::vector<OctreeElementWithNeigh>& elements_out, MPI_Comm comm){

    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    std::vector<OctreeElementWithNeigh> initial_elements;        // before SFC
    GetInitialElementsDistributionOct(mesh_file_path, initial_elements, MPI_COMM_WORLD);
    SetMortonEncoding(initial_elements,ElementType::HEX,MPI_COMM_WORLD);
    elements_out.resize(initial_elements.size());
    if (! my_rank)
    {
        print_log("starting sfc sort");
    }

    MPI_Barrier(MPI_COMM_WORLD);
    std::vector<OctreeElementWithNeigh> initial_elements_cpy(initial_elements);
    // par::sampleSort<T>(initial_elements_cpy,elements_out,MPI_COMM_WORLD);       //warmup?
    SampleSortMorton(initial_elements_cpy,elements_out,MPI_COMM_WORLD);     // warmup

    MPI_Barrier(MPI_COMM_WORLD);
    auto start = std::chrono::high_resolution_clock::now();
    // par::sampleSort<T>(initial_elements,elements_out,MPI_COMM_WORLD);
    SampleSortMorton(initial_elements,elements_out,MPI_COMM_WORLD);
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

void ResolveOctreeElementConnectivity(std::vector<OctreeElementWithNeigh> &local_elements,
                                std::vector<uint64_t> &proc_element_counts,
                                std::vector<uint64_t> &proc_element_counts_scanned,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &connected_element_pairs_out,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &boundary_connected_element_pairs_out,
                                MPI_Comm comm)
{
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    uint64_t global_element_count = proc_element_counts[procs_n-1]+ proc_element_counts_scanned[procs_n-1];
    uint64_t local_element_count = local_elements.size();

    std::vector<int> proc_element_counts_(proc_element_counts.begin(), proc_element_counts.end());
    std::vector<int> proc_element_counts_scanned_(proc_element_counts_scanned.begin(),proc_element_counts_scanned.end());


    std::vector<OctreeElementWithNeigh> global_all_elements(proc_element_counts[procs_n-1]+ proc_element_counts_scanned[procs_n-1]);


    //collecting all elements
    par::Mpi_Allgatherv(local_elements.data(),static_cast<int>(local_element_count),global_all_elements.data(),proc_element_counts_.data(),proc_element_counts_scanned_.data(),comm);

    std::unordered_map<uint64_t,uint64_t> elemtag_globidx_map;
    elemtag_globidx_map.reserve(global_element_count);

    for (auto elem : global_all_elements)
    {
        elemtag_globidx_map[elem.element_tag] = elem.global_idx;
        
    }
    




    
    std::vector<std::pair<uint64_t, uint64_t>> local_connected;
    std::vector<std::pair<uint64_t, uint64_t>> boundary_connected;


    

    for (size_t i = 0; i < global_element_count; i++)
    {
        bool is_my_i = global_all_elements[i].global_idx >= proc_element_counts_scanned[my_rank] && (global_all_elements[i].global_idx < (proc_element_counts_scanned[my_rank] + proc_element_counts[my_rank]));
        for (int n_i = 0; n_i < 6; n_i++)
        {
            if (global_all_elements[i].neigh[n_i] == OCTREE_BOUNDARY)
            {
                continue;
            }
            
            
            uint64_t j = elemtag_globidx_map[global_all_elements[i].neigh[n_i]];
            assert(i!=j);

            bool is_my_j = j>= proc_element_counts_scanned[my_rank] && (j < (proc_element_counts_scanned[my_rank] + proc_element_counts[my_rank]));

            if(!is_my_i && !is_my_j)        //edge outside my range (i.e. not local, not ghost)
            {
                continue;
            }

            if (is_my_i && is_my_j)  //local edge
            {

                local_connected.push_back({std::min(i, j),std::max(i, j)});

                
            }else if(is_my_i)      // boundary edge 
            {   
                boundary_connected.push_back({i,j});

            }else if(is_my_j)      // boundary edge 
            {   
                boundary_connected.push_back({j,i});

            }
            
            

        }
        

    }
    omp_par::merge_sort(&local_connected[0], &local_connected[local_connected.size()],
                        [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) { 
                            return a.second < b.second; 
                                
                        });

    stable_sort(&local_connected[0], &local_connected[local_connected.size()],
                        [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) { 
                            return a.first < b.first; 
                                
                        });


    omp_par::merge_sort(&boundary_connected[0], &boundary_connected[boundary_connected.size()],
                        [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) { 
                            return a.second < b.second; 
                                
                        });

    stable_sort(&boundary_connected[0], &boundary_connected[boundary_connected.size()],
                        [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) { 
                            return a.first < b.first; 
                                
                        });

    
    connected_element_pairs_out.clear();
    {
        std::pair<ElementWithTag, ElementWithTag> pair;
        pair.first.element_tag = global_all_elements[local_connected[0].first].element_tag;
        pair.first.global_idx = global_all_elements[local_connected[0].first].global_idx;

        pair.second.element_tag = global_all_elements[local_connected[0].second].element_tag;
        pair.second.global_idx = global_all_elements[local_connected[0].second].global_idx;
        connected_element_pairs_out.push_back(pair);

        for (size_t pair_i = 1; pair_i < local_connected.size(); pair_i++)
        {
            if ((local_connected[pair_i].first == local_connected[pair_i-1].first) && (local_connected[pair_i].second == local_connected[pair_i-1].second))
            {
                // duplicate local edge found
                continue;
            }
            pair.first.element_tag = global_all_elements[local_connected[pair_i].first].element_tag;
            pair.first.global_idx = global_all_elements[local_connected[pair_i].first].global_idx;

            pair.second.element_tag = global_all_elements[local_connected[pair_i].second].element_tag;
            pair.second.global_idx = global_all_elements[local_connected[pair_i].second].global_idx;
            connected_element_pairs_out.push_back(pair);
            
        }

    }

    boundary_connected_element_pairs_out.clear();
    {
        std::pair<ElementWithTag, ElementWithTag> pair;
        pair.first.element_tag = global_all_elements[boundary_connected[0].first].element_tag;
        pair.first.global_idx = global_all_elements[boundary_connected[0].first].global_idx;

        pair.second.element_tag = global_all_elements[boundary_connected[0].second].element_tag;
        pair.second.global_idx = global_all_elements[boundary_connected[0].second].global_idx;
        boundary_connected_element_pairs_out.push_back(pair);

        for (size_t pair_i = 1; pair_i < boundary_connected.size(); pair_i++)
        {
            if ((boundary_connected[pair_i].first == boundary_connected[pair_i-1].first) && (boundary_connected[pair_i].second == boundary_connected[pair_i-1].second))
            {
                // duplicate bdry edge found
                continue;
            }
            pair.first.element_tag = global_all_elements[boundary_connected[pair_i].first].element_tag;
            pair.first.global_idx = global_all_elements[boundary_connected[pair_i].first].global_idx;

            pair.second.element_tag = global_all_elements[boundary_connected[pair_i].second].element_tag;
            pair.second.global_idx = global_all_elements[boundary_connected[pair_i].second].global_idx;
            boundary_connected_element_pairs_out.push_back(pair);
            
        }

    }

}

