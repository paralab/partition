
#include "dist-graph.hpp"
#include "../mesh-util/mesh-util.hpp"
#include "../usort/ompUtils.h"
#include "../usort/parUtils.h"

#include "mpi.h"
#include "../metis-util/metis-util.hpp"

#include <chrono>
#include <unordered_map>


// Overloading the << operator for BFSValue
std::ostream& operator<<(std::ostream& os, const BFSValue& obj) {
    os << "(" << obj.distance << "," << obj.label  << ")";
    // os << obj.global_idx;
    
    return os;
}

// Overloading the << operator for PageRankValue
std::ostream& operator<<(std::ostream& os, const PageRankValue& obj) {
    os << "(" << obj.value << "," << obj.label  << ")";
    // os << obj.global_idx;
    
    return os;
}

/**
 * Assumes the global_idx of elements should be in a contiguous range
 * For rank i, global_idx of its elements should should be in interval
 * [proc_element_counts_scanned[i], proc_element_counts_scanned[i] +
 * proc_element_counts[i] -1 )
 * 
 * both own_elements and ghost_elements should be given sorted by global_idx. ascending order
 * both local_connectivity and boundary_connectivity_cpy should not contain duplicate edges. 
 * (i.e.) if u and v are connected, either (u,v) or (v,u) should be in connectivity vectors
 * in boudary connectivity, for each (u,v) u should belong to this process and v should belong to another process
 */
DistGraph::DistGraph(const std::vector<ElementWithCoord>& own_elements,const std::vector<ElementWithTag>& ghost_elements,
                     const std::vector<std::pair<ElementWithTag, ElementWithTag>>& local_connectivity,
                     const std::vector<std::pair<ElementWithTag, ElementWithTag>>& boundary_connectivity,
                     const std::vector<uint64_t>& proc_element_counts,
                     const std::vector<uint64_t>& proc_element_counts_scanned, 
                     const std::vector<int>& ghost_element_counts,
                     MPI_Comm comm) {
    this->comm = comm;
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    this->local_count = own_elements.size();
    this->ghost_count = ghost_elements.size();

    this->vtx_dist.assign(proc_element_counts_scanned.begin(), proc_element_counts_scanned.end());

    this->vtx_dist.push_back(proc_element_counts_scanned[procs_n-1] + proc_element_counts[procs_n-1]);

    // print_log("[", my_rank, "]: vtx_dist ", VectorToString(vtx_dist));




    this->local_degrees.resize(own_elements.size() + ghost_elements.size());
    std::fill(this->local_degrees.begin(), this->local_degrees.end(),0);
    this->local_xdj.resize(own_elements.size() + ghost_elements.size() + 1);
    std::fill(this->local_xdj.begin(), this->local_xdj.end(),0);

    this->local_adjncy.resize(2 * (local_connectivity.size() + boundary_connectivity.size()));
    this->dist_adjncy.resize(2*local_connectivity.size()  + boundary_connectivity.size());
    this->ghost_counts.assign(ghost_element_counts.begin(), ghost_element_counts.end());

    this->ghost_counts_scanned.resize(procs_n);
    omp_par::scan(&this->ghost_counts[0], &this->ghost_counts_scanned[0], procs_n);



    // std::vector<uint64_t> vertex_degrees(own_elements.size() + ghost_elements.size());
    for (auto& edge : local_connectivity) {
        this->local_degrees[edge.first.global_idx - proc_element_counts_scanned[my_rank]]++;         // this assumes global_idx is in correct range
        this->local_degrees[edge.second.global_idx - proc_element_counts_scanned[my_rank]]++;         // this assumes global_idx is in correct range
    
    }
    std::vector<std::pair<ElementWithTag, ElementWithTag>> boundary_connectivity_cpy(boundary_connectivity);
    // print_log("[", my_rank, "]: ghost_elements ", VectorToString(ghost_elements));
    // print_log("[", my_rank, "]: local_elements ", VectorToString(own_elements));

    // print_log("[", my_rank, "]: local_connectivity ", VectorToString(local_connectivity));
    // print_log("[", my_rank, "]: boundary_connectivity_cpy ", VectorToString(boundary_connectivity_cpy));


    //sort boundary edges by other process elements (i.e. the second element)
    //this helps in adding edges for local ordering
    
    omp_par::merge_sort(&boundary_connectivity_cpy[0], &boundary_connectivity_cpy[boundary_connectivity_cpy.size()],
                        [](const std::pair<ElementWithTag, ElementWithTag>& a, const std::pair<ElementWithTag, ElementWithTag>& b) { 
                            return a.second.global_idx < b.second.global_idx; 
                                
                        });
    // print_log("[", my_rank, "]: boundary_connectivity_cpy sorted ", VectorToString(boundary_connectivity_cpy));
    
    if(!boundary_connectivity_cpy.empty()){
        uint64_t current_ghost_element_local_index = own_elements.size();
        //assumes first element belongs ro this process and second element belongs to another process
        this->local_degrees[boundary_connectivity_cpy[0].first.global_idx - proc_element_counts_scanned[my_rank]]++;
        this->local_degrees[current_ghost_element_local_index]++;
        

        for (int boundary_edge_i =1; boundary_edge_i <boundary_connectivity_cpy.size(); boundary_edge_i++)
        {
            // since we sorted boundary edges by second element
            if (boundary_connectivity_cpy[boundary_edge_i].second.global_idx != boundary_connectivity_cpy[boundary_edge_i-1].second.global_idx)
            {
                current_ghost_element_local_index++;
            }
            
            //assumes first element belongs ro this process and second element belongs to another process
            this->local_degrees[boundary_connectivity_cpy[boundary_edge_i].first.global_idx - proc_element_counts_scanned[my_rank]]++;
            this->local_degrees[current_ghost_element_local_index]++;
        }        
    }

    // print_log("[", my_rank, "]: local_degrees ", VectorToString(this->local_degrees));

    omp_par::scan(&this->local_degrees[0], &this->local_xdj[0], own_elements.size() + ghost_elements.size());
    //scan operation does not populate the last extra entry for CSR encoding, hence we have to manually populate the last element
    this->local_xdj[own_elements.size() + ghost_elements.size()] = 
            this->local_xdj[own_elements.size() + ghost_elements.size() -1 ] + this->local_degrees[own_elements.size() + ghost_elements.size()-1];
    // print_log("[", my_rank, "]: local_xdj ", VectorToString(local_xdj));

    /**
     * populating adjacency structure
    */
    std::vector<uint64_t> next_index;
    next_index.assign(this->local_xdj.begin(), this->local_xdj.end()-1);

    for (auto& edge : local_connectivity) {
        auto local_index_1 = edge.first.global_idx - proc_element_counts_scanned[my_rank];
        auto local_index_2 = edge.second.global_idx - proc_element_counts_scanned[my_rank];

        this->local_adjncy[next_index[local_index_1]] = local_index_2;
        this->local_adjncy[next_index[local_index_2]] = local_index_1;

        this->dist_adjncy[next_index[local_index_1]] = edge.second.global_idx;
        this->dist_adjncy[next_index[local_index_2]] = edge.first.global_idx;


        next_index[local_index_1]++;
        next_index[local_index_2]++;
    }
    // print_log("[", my_rank, "]: local_xdj done for local elements");


    if(!boundary_connectivity_cpy.empty()){
        uint64_t current_ghost_element_local_index = own_elements.size();
        //assumes first element belongs to this process and second element belongs to another process
        {
            auto local_index_1 = boundary_connectivity_cpy[0].first.global_idx - proc_element_counts_scanned[my_rank];
            auto local_index_2 = current_ghost_element_local_index;
            this->local_adjncy[next_index[local_index_1]] = local_index_2;
            this->local_adjncy[next_index[local_index_2]] = local_index_1;

            this->dist_adjncy[next_index[local_index_1]] = boundary_connectivity_cpy[0].second.global_idx;

            next_index[local_index_1]++;
            next_index[local_index_2]++;
        }

        for (int boundary_edge_i =1; boundary_edge_i <boundary_connectivity_cpy.size(); boundary_edge_i++)
        {
            //since we sorted boundary edges by second element
            if (boundary_connectivity_cpy[boundary_edge_i].second.global_idx != boundary_connectivity_cpy[boundary_edge_i-1].second.global_idx)
            {
                current_ghost_element_local_index++;
            }
            
            //assumes first element belongs ro this process and second element belongs to another process
            auto local_index_1 = boundary_connectivity_cpy[boundary_edge_i].first.global_idx - proc_element_counts_scanned[my_rank];
            auto local_index_2 = current_ghost_element_local_index;
            this->local_adjncy[next_index[local_index_1]] = local_index_2;
            this->local_adjncy[next_index[local_index_2]] = local_index_1;

            this->dist_adjncy[next_index[local_index_1]] = boundary_connectivity_cpy[boundary_edge_i].second.global_idx;



            next_index[local_index_1]++;
            next_index[local_index_2]++;
        }





    }
    /**
     * building the scatter map
     * */     
    this->sending_scatter_map.clear();  
    this->send_counts.resize(procs_n);
    std::fill(this->send_counts.begin(), this->send_counts.end(),0);
    this->send_counts_scanned.resize(procs_n);
    std::fill(this->send_counts_scanned.begin(), this->send_counts_scanned.end(),0);
    if(!boundary_connectivity_cpy.empty()){
        uint64_t current_other_proc = 0;

        std::vector<ElementWithTag> send_elements_with_dups;
        std::vector<uint64_t> send_elements_with_dups_counts(procs_n, 0);
        std::vector<uint64_t> send_elements_with_dups_counts_scanned(procs_n, 0);


        for (auto& edge : boundary_connectivity_cpy){
            if (edge.second.global_idx >= proc_element_counts_scanned[current_other_proc] &&
                edge.second.global_idx < proc_element_counts_scanned[current_other_proc] + proc_element_counts[current_other_proc])
            {
                send_elements_with_dups.push_back(edge.first);
                send_elements_with_dups_counts[current_other_proc]++;
            }else
            {
                while (1)
                {
                    current_other_proc++;
                    if (edge.second.global_idx >= proc_element_counts_scanned[current_other_proc] &&
                        edge.second.global_idx < proc_element_counts_scanned[current_other_proc] + proc_element_counts[current_other_proc])
                    {
                        send_elements_with_dups.push_back(edge.first);
                        send_elements_with_dups_counts[current_other_proc]++;
                        break;
                    }
                }
                
            }          
            
        }
        omp_par::scan(&send_elements_with_dups_counts[0], &send_elements_with_dups_counts_scanned[0], procs_n);

        std::vector<ElementWithTag> send_elements;

        //sort each send elements per each process
        for (size_t proc_i = 0; proc_i < procs_n; proc_i++)
        {
            if (send_elements_with_dups_counts[proc_i])
            {
                omp_par::merge_sort(&send_elements_with_dups[send_elements_with_dups_counts_scanned[proc_i]], 
                                    &send_elements_with_dups[send_elements_with_dups_counts_scanned[proc_i]+send_elements_with_dups_counts[proc_i]]);
            
                send_elements.push_back(send_elements_with_dups[send_elements_with_dups_counts_scanned[proc_i]]);
                send_counts[proc_i]++;
                for (size_t elem_i = send_elements_with_dups_counts_scanned[proc_i]+1; elem_i < (send_elements_with_dups_counts_scanned[proc_i] + send_elements_with_dups_counts[proc_i]); elem_i++)
                {
                    if (send_elements_with_dups[elem_i].global_idx != send_elements_with_dups[elem_i-1].global_idx)
                    {
                        send_elements.push_back(send_elements_with_dups[elem_i]);
                        send_counts[proc_i]++;
                    }
                    
                }
                
            }    


        }
        

        // print_log("[", my_rank, "]: send_elements ", VectorToString(send_elements));
        omp_par::scan(&this->send_counts[0], &this->send_counts_scanned[0], procs_n);
        this->send_count = this->send_counts_scanned[procs_n-1]+ this->send_counts[procs_n-1];
        this->sending_scatter_map.resize(this->send_count);

        #pragma omp parallel for
        for (size_t send_elem_i = 0; send_elem_i < send_elements.size(); send_elem_i++)
        {
            this->sending_scatter_map[send_elem_i] = send_elements[send_elem_i].global_idx - proc_element_counts_scanned[my_rank];
        }
        

    }


    

}

// DistGraph::~DistGraph(){
//     // this->local_xdj.clear();
//     // this->local_adjncy.clear();

// }


std::string DistGraph::PrintLocal(){
    std::ostringstream output;
    for (size_t vertex_i = 0; vertex_i < this->local_xdj.size()-1; vertex_i++)
    {
        output << vertex_i << "\t->";
        for (size_t neigh_i = this->local_xdj[vertex_i]; neigh_i < this->local_xdj[vertex_i+1]; neigh_i++)
        {
            output << this->local_adjncy[neigh_i] << ",";
        }
        output << "\n";
        
    }

    return output.str();
    
}

std::string DistGraph::PrintDist(){
    std::ostringstream output;
    for (size_t vertex_i = 0; vertex_i < this->local_count; vertex_i++)
    {
        output << vertex_i << "\t->";
        for (size_t neigh_i = this->local_xdj[vertex_i]; neigh_i < this->local_xdj[vertex_i+1]; neigh_i++)
        {
            output << this->dist_adjncy[neigh_i] << ",";
        }
        output << "\n";
        
    }

    return output.str();
    
}


PartitionStatus DistGraph::PartitionBFS(std::vector<uint16_t>& partition_labels_out){
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);
    std::vector<BFSValue> bfs_vector(this->local_xdj.size()-1, {.label = DIST_GRAPH_BFS_NO_LABEL, .distance =  DIST_GRAPH_BFS_INFINITY});

    //temp buffers to be used in each local iteration in BFS. (to minimize re-allocation overhead in each iteration)
    std::vector<BFSValue> bfs_vector_tmp(bfs_vector.size());
    std::vector<bool> vector_diff(bfs_vector.size());
    
    std::vector<BFSValue> ghost_send_buffer(this->send_count);
    std::vector<BFSValue> ghost_recv_buffer(this->ghost_count);

    std::vector<bool> ghost_is_not_stable(this->ghost_count);


    /**
     * using sfc seeds
     * elements are already ordered to morton SFC
     * get the 'middle' local element as seed
     * TODO: oversampling can be implemented here
    */

    bfs_vector[this->local_count/2].distance = 0;
    bfs_vector[this->local_count/2].label = my_rank;
    bool is_not_stable_global = true;      // global BFS stability
    int round_counter = 0;
    MPI_Barrier(this->comm);
    auto com_duration = std::chrono::microseconds(0);
    auto reduce_duration = std::chrono::microseconds(0);

    auto start = std::chrono::high_resolution_clock::now();
    while (is_not_stable_global)
    {
        // if (!my_rank)
        // {
        //     print_log("BFS round: ", ++round_counter);
        // }
        round_counter++;
        
        is_not_stable_global = false;
        bool is_not_stable_local = this->RunLocalMultiBFSToStable(bfs_vector, bfs_vector_tmp, vector_diff);
        // print_log("[", my_rank, "]: BFS iteration done");

        #pragma omp parallel for
        for (size_t send_i = 0; send_i < this->send_count; send_i++)
        {
            ghost_send_buffer[send_i] = bfs_vector[this->sending_scatter_map[send_i]];
        }
        auto com_start = std::chrono::high_resolution_clock::now();
        //ghost exchange


        // MPI_Alltoallv(ghost_send_buffer.data(),
        //             this->send_counts.data(), this->send_counts_scanned.data(), par::Mpi_datatype<BFSValue>::value(), 
        //             ghost_recv_buffer.data(), this->ghost_counts.data(), this->ghost_counts_scanned.data(),
        //             par::Mpi_datatype<BFSValue>::value(), comm);

        MPI_Barrier(this->comm);
        par::Mpi_Alltoallv_sparse(ghost_send_buffer.data(), this->send_counts.data(), this->send_counts_scanned.data(), 
                ghost_recv_buffer.data(), this->ghost_counts.data(), this->ghost_counts_scanned.data(), comm);

        MPI_Barrier(this->comm);

        auto com_end = std::chrono::high_resolution_clock::now();
        com_duration += std::chrono::duration_cast<std::chrono::microseconds>(com_end - com_start);

        //ghost update for received values
        #pragma omp parallel for
        for (size_t recv_i = 0; recv_i < this->ghost_count; recv_i++)
        {
            ghost_is_not_stable[recv_i] = false;
            auto offset = this->local_count;        // ghost elements are in the last section of the vector, in sorted order
            if (bfs_vector[offset+recv_i].distance > ghost_recv_buffer[recv_i].distance)
            {
                bfs_vector[offset+recv_i].distance = ghost_recv_buffer[recv_i].distance;
                bfs_vector[offset+recv_i].label = ghost_recv_buffer[recv_i].label;
                // is_not_stable_local = true;
                ghost_is_not_stable[recv_i] = true;
            }
            
        }

        bool ghost_any_is_not_stable = false;
        #pragma omp parallel for reduction(||:ghost_any_is_not_stable)
        for (size_t recv_i = 0; recv_i < this->ghost_count; recv_i++) {
            ghost_any_is_not_stable = ghost_any_is_not_stable|| ghost_is_not_stable[recv_i]; 
        }

        is_not_stable_local = is_not_stable_local || ghost_any_is_not_stable;

        auto reduce_com_start = std::chrono::high_resolution_clock::now();

        MPI_Allreduce(&is_not_stable_local,&is_not_stable_global,1,MPI_CXX_BOOL,MPI_LOR,this->comm);
        auto reduce_com_end = std::chrono::high_resolution_clock::now();
        reduce_duration += std::chrono::duration_cast<std::chrono::microseconds>(reduce_com_end - reduce_com_start);
        // print_log("[", my_rank, "]: BFS vector", VectorToString(bfs_vector));

        
    } 

    // print_log("[", my_rank, "]: BFS done");
    // print_log("[", my_rank, "]: BFS vector", VectorToString(bfs_vector));

    MPI_Barrier(this->comm);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (!my_rank)
    {
        print_log("BFS sync rounds: ", round_counter);
        print_log("BFS comm (ghost exchange) time:\t", com_duration.count(), " us");
        print_log("BFS comm (reduce) time:\t\t", reduce_duration.count(), " us");

        print_log("BFS total time:\t\t\t", duration.count(), " us");
    }

    partition_labels_out.resize(this->local_count);


    #pragma omp parallel for
    for (size_t local_i = 0; local_i < this->local_count; local_i++)
    {
        partition_labels_out[local_i] = bfs_vector[local_i].label;
    }

    return {.return_code = 0, .time_ms = duration.count()};

}

bool DistGraph::RunLocalMultiBFSToStable(std::vector<BFSValue>& bfs_vector, std::vector<BFSValue>& bfs_vector_tmp, std::vector<bool> vector_diff){
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);
    bool changed = false;       // to detect if the BFS incremented

    bool is_not_stable = true; // to detect if the BFS incremented in each increment
    // std::vector<BFSValue> bfs_vector_tmp(bfs_vector.size());
    // std::vector<bool> vector_diff(bfs_vector.size());
    while (is_not_stable)
    {
        // print_log(VectorToString(multi_bfs_distances));

        is_not_stable = false;
        #pragma omp parallel
        {
            // print_log("thread count ", omp_get_num_threads());
            #pragma omp for
            for (graph_indexing_t v_i = 0; v_i < bfs_vector.size(); v_i++)
            {
                // bfs_status_new_temp[v_i] = NULL;
                auto best_distance = bfs_vector[v_i].distance;
                auto best_label = bfs_vector[v_i].label;
                vector_diff[v_i] = false;
                // print_log(best_distance, best_label);
                for (graph_indexing_t neighbor_i = this->local_xdj[v_i]; neighbor_i < this->local_xdj[v_i+1];neighbor_i++)
                {
                    auto neighbor = local_adjncy[neighbor_i];
                    if (bfs_vector[neighbor].label == DIST_GRAPH_BFS_NO_LABEL)
                    {
                        continue;
                    }
                    // print_log(multi_bfs_labels[neighbor_i]);

                    if (best_distance > (bfs_vector[neighbor].distance + 1))
                    {
                        best_distance = bfs_vector[neighbor].distance + 1;
                        best_label = bfs_vector[neighbor].label;
                        vector_diff[v_i] = true;


                    }
                }
                bfs_vector_tmp[v_i].distance = best_distance;
                bfs_vector_tmp[v_i].label = best_label;

            }
            #pragma omp for
            for (size_t v_i = 0; v_i < bfs_vector.size(); v_i++){
                bfs_vector[v_i].distance = bfs_vector_tmp[v_i].distance;
                bfs_vector[v_i].label = bfs_vector_tmp[v_i].label;
            }
            #pragma omp for reduction(||:is_not_stable)
            for (size_t v_i = 0; v_i < bfs_vector.size(); v_i++) {
                is_not_stable = is_not_stable|| vector_diff[v_i]; // Perform logical OR operation
            }

            changed = changed || is_not_stable;



        }
    }

    return changed;
}

/**
 * routine to collect the vextex degrees, including ghost degrees of ghost vertices
 * ghost vertex degrees are required to correctly run the pagerank
*/
void DistGraph::GetVertexDegrees(std::vector<graph_indexing_t>& degrees_out){
    int my_rank;
    MPI_Comm_rank(this->comm, &my_rank);
    degrees_out.resize(this->local_count + this->ghost_count);
    // std::fill(degrees_out.begin(), degrees_out.end(),0);

    #pragma omp parallel for
    for (size_t local_i = 0; local_i < this->local_count ; local_i++)
    {
        degrees_out[local_i] = this->local_xdj[local_i+1] - this->local_xdj[local_i];
    }
    std::vector<graph_indexing_t> send_buffer(this->send_count);
    #pragma omp parallel for
    for (size_t send_i = 0; send_i < this->send_count; send_i++)
    {
        send_buffer[send_i] = degrees_out[this->sending_scatter_map[send_i]];
    }
    MPI_Barrier(this->comm);
    par::Mpi_Alltoallv_sparse(send_buffer.data(), this->send_counts.data(), this->send_counts_scanned.data(), 
            &degrees_out[this->local_count], this->ghost_counts.data(), this->ghost_counts_scanned.data(), comm);
    MPI_Barrier(this->comm);
    // print_log("[", my_rank, "]: degrees", VectorToString(degrees_out));

}



void DistGraph::PartitionPageRank(std::vector<uint16_t>& partition_labels_out){
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);
    std::vector<PageRankValue> pagerank_vector(this->local_xdj.size()-1, {.label = DIST_GRAPH_PAGERANK_NO_LABEL, .value =  0});
    std::vector<PageRankValue> ghost_send_buffer(this->send_count);
    std::vector<PageRankValue> ghost_recv_buffer(this->ghost_count);

    std::vector<graph_indexing_t> vertex_degrees(this->local_count + this->ghost_count);
    this->GetVertexDegrees(vertex_degrees);

    std::vector<bool> ghost_is_not_stable(this->ghost_count);


    /**
     * using sfc seeds
     * elements are already ordered to morton SFC
     * get the 'middle' local element as seed
     * TODO: oversampling can be implemented here
    */

   const float pagerank_initial_value = 10000;
   const float min_relative_change = 1.05;

    pagerank_vector[this->local_count/2].value = pagerank_initial_value;
    pagerank_vector[this->local_count/2].label = my_rank;
    bool is_not_stable_global = true;      // global PageRank stability
    int round_counter = 0;
    MPI_Barrier(this->comm);
    auto com_duration = std::chrono::microseconds(0);
    auto reduce_duration = std::chrono::microseconds(0);

    auto start = std::chrono::high_resolution_clock::now();
    while (is_not_stable_global)
    {
        if (!my_rank)
        {
            print_log("PageRank round: ", ++round_counter);
        }
        
        is_not_stable_global = false;
        bool is_not_stable_local = this->RunLocalMultiPageRankToStable(pagerank_vector,vertex_degrees, min_relative_change);
        // print_log("[", my_rank, "]: PageRank iteration done");

        #pragma omp parallel for
        for (size_t send_i = 0; send_i < this->send_count; send_i++)
        {
            ghost_send_buffer[send_i] = pagerank_vector[this->sending_scatter_map[send_i]];
        }
        auto com_start = std::chrono::high_resolution_clock::now();
        //ghost exchange


        // MPI_Alltoallv(ghost_send_buffer.data(),
        //             this->send_counts.data(), this->send_counts_scanned.data(), par::Mpi_datatype<PageRankValue>::value(), 
        //             ghost_recv_buffer.data(), this->ghost_counts.data(), this->ghost_counts_scanned.data(),
        //             par::Mpi_datatype<PageRankValue>::value(), comm);

        MPI_Barrier(this->comm);
        par::Mpi_Alltoallv_sparse(ghost_send_buffer.data(), this->send_counts.data(), this->send_counts_scanned.data(), 
                ghost_recv_buffer.data(), this->ghost_counts.data(), this->ghost_counts_scanned.data(), comm);

        MPI_Barrier(this->comm);

        auto com_end = std::chrono::high_resolution_clock::now();
        com_duration += std::chrono::duration_cast<std::chrono::microseconds>(com_end - com_start);

        //ghost update for received values
        #pragma omp parallel for
        for (size_t recv_i = 0; recv_i < this->ghost_count; recv_i++)
        {
            ghost_is_not_stable[recv_i] = false;
            auto offset = this->local_count;        // ghost elements are in the last section of the vector, in sorted order
            if ((ghost_recv_buffer[recv_i].value/pagerank_vector[offset+recv_i].value) > min_relative_change )
            {
                pagerank_vector[offset+recv_i].value = ghost_recv_buffer[recv_i].value;
                pagerank_vector[offset+recv_i].label = ghost_recv_buffer[recv_i].label;
                // is_not_stable_local = true;
                ghost_is_not_stable[recv_i] = true;
            }
            
        }

        bool ghost_any_is_not_stable = false;
        #pragma omp parallel for reduction(||:ghost_any_is_not_stable)
        for (size_t recv_i = 0; recv_i < this->ghost_count; recv_i++) {
            ghost_any_is_not_stable = ghost_any_is_not_stable|| ghost_is_not_stable[recv_i]; 
        }

        is_not_stable_local = is_not_stable_local || ghost_any_is_not_stable;

        auto reduce_com_start = std::chrono::high_resolution_clock::now();

        MPI_Allreduce(&is_not_stable_local,&is_not_stable_global,1,MPI_CXX_BOOL,MPI_LOR,this->comm);
        auto reduce_com_end = std::chrono::high_resolution_clock::now();
        reduce_duration += std::chrono::duration_cast<std::chrono::microseconds>(reduce_com_end - reduce_com_start);
        // print_log("[", my_rank, "]: PageRank vector", VectorToString(pagerank_vector));

        
    } 

    // print_log("[", my_rank, "]: PageRank done");
    // print_log("[", my_rank, "]: PageRank vector", VectorToString(pagerank_vector));

    MPI_Barrier(this->comm);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (!my_rank)
    {
        print_log("PageRank comm time:\t", com_duration.count(), " us");
        print_log("PageRank comm (reduce) time:\t", reduce_duration.count(), " us");

        print_log("PageRank time:\t", duration.count(), " us");
    }

    partition_labels_out.resize(this->local_count);


    #pragma omp parallel for
    for (size_t local_i = 0; local_i < this->local_count; local_i++)
    {
        partition_labels_out[local_i] = pagerank_vector[local_i].label;
    }

}


bool DistGraph::RunLocalMultiPageRankToStable(std::vector<PageRankValue>& pagerank_vector,
                std::vector<graph_indexing_t> vertex_degrees, const float min_relative_change){
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);
    bool changed = false;       // to detect if the pagerank continoued at least once

    bool is_not_stable = true; // to detect if the pagerank continoued incremented in each increment
    std::vector<PageRankValue> pagerank_vector_tmp(pagerank_vector.size());
    std::vector<bool> vector_diff(pagerank_vector.size());
    while (is_not_stable)
    {
        if (!my_rank)
        {
            // print_log("[", my_rank, "]: PageRank vector", VectorToString(pagerank_vector));
            // print_log(pagerank_vector[0]);
        }
        

        is_not_stable = false;
        #pragma omp parallel
        {
            // print_log("thread count ", omp_get_num_threads());
            #pragma omp for
            for (graph_indexing_t v_i = 0; v_i < pagerank_vector.size(); v_i++)
            {
                // bfs_status_new_temp[v_i] = NULL;
                auto curr_value = pagerank_vector[v_i].value;
                auto curr_label = pagerank_vector[v_i].label;
                vector_diff[v_i] = false;
                // print_log(best_distance, best_label);
                std::unordered_map<decltype(PageRankValue::label),decltype(PageRankValue::value)> incoming;
                incoming.reserve(this->local_xdj[v_i+1] - this->local_xdj[v_i]);
                
                for (graph_indexing_t neighbor_i = this->local_xdj[v_i]; neighbor_i < this->local_xdj[v_i+1];neighbor_i++)
                {
                    auto neighbor = local_adjncy[neighbor_i];
                    if (pagerank_vector[neighbor].label == DIST_GRAPH_PAGERANK_NO_LABEL)
                    {
                        continue;
                    }
                    // print_log(multi_bfs_labels[neighbor_i]);
                    // float neighbor_contribution_factor = 1.0/(this->local_xdj[neighbor+1] - this->local_xdj[neighbor]);
                    float neighbor_contribution_factor = 1.0/vertex_degrees[neighbor];
                    // incoming[pagerank_vector[neighbor].label]+= (pagerank_vector[neighbor].value)*neighbor_contribution_factor;

                    if (incoming.find(pagerank_vector[neighbor].label) == incoming.end())
                    {
                        incoming[pagerank_vector[neighbor].label]= (pagerank_vector[neighbor].value)*neighbor_contribution_factor;
                    }else
                    {
                        incoming[pagerank_vector[neighbor].label]+= (pagerank_vector[neighbor].value)*neighbor_contribution_factor;
                    }
                    
                    
                    
                }
                pagerank_vector_tmp[v_i].label = curr_label;
                pagerank_vector_tmp[v_i].value = curr_value;

                if (incoming.size() == 0)
                {

                    continue;
                }
                
                auto best_incoming = get_max(incoming);
                // if(!my_rank) print_log(best_incoming.first, best_incoming.second);
                if (curr_label == DIST_GRAPH_PAGERANK_NO_LABEL)
                {
                    pagerank_vector_tmp[v_i].label = best_incoming.first;
                    pagerank_vector_tmp[v_i].value = best_incoming.second;
                    vector_diff[v_i] = true;

                } else if ((best_incoming.second/curr_value) > min_relative_change)
                {
                    pagerank_vector_tmp[v_i].label = best_incoming.first;
                    pagerank_vector_tmp[v_i].value = best_incoming.second;
                    vector_diff[v_i] = true;
                }

            }
            #pragma omp for
            for (size_t v_i = 0; v_i < pagerank_vector.size(); v_i++){
                pagerank_vector[v_i].value = pagerank_vector_tmp[v_i].value;
                pagerank_vector[v_i].label = pagerank_vector_tmp[v_i].label;
            }
            #pragma omp for reduction(||:is_not_stable)
            for (size_t v_i = 0; v_i < pagerank_vector.size(); v_i++) {
                is_not_stable = is_not_stable|| vector_diff[v_i]; // Perform logical OR operation
            }

            changed = changed || is_not_stable;



        }
    }

    return changed;
}

PartitionStatus DistGraph::PartitionParmetis(std::vector<uint16_t>& partition_labels_out){
    int procs_n;
    MPI_Comm_size(this->comm, &procs_n);
    std::vector<uint64_t> dist_xadj(this->local_xdj.begin(), this->local_xdj.begin()+ (this->local_count + 1));
    return GetParMETISPartitions(this->vtx_dist,dist_xadj,this->dist_adjncy,this->local_count,procs_n,partition_labels_out,this->comm);
}

/**
 * given a partition labelling, get global partition sizes and boundary sizes
 * result will be present in partition_sizes_out and partition_boundaries_out, only in MPI rank 0 process.
*/
void DistGraph::GetPartitionMetrics(std::vector<uint16_t>& local_partition_labels,
                                    std::vector<uint32_t>& partition_sizes_out,
                                    std::vector<uint32_t>& partition_boundaries_out) {
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);
    assert(local_partition_labels.size() == this->local_count);

    std::vector<uint32_t> local_partition_sizes(procs_n, 0);
    for (size_t local_i = 0; local_i < this->local_count; local_i++) {
        local_partition_sizes[local_partition_labels[local_i]]++;
    }

    /**
     * exchanging partition labels of ghost vertices to calculate partition boundary counts
     */
    std::vector<uint16_t> local_and_ghost_partition_labels(this->local_count + this->ghost_count);
    std::copy(local_partition_labels.begin(), local_partition_labels.end(), local_and_ghost_partition_labels.begin());

    std::vector<uint16_t> ghost_send_buffer(this->send_count);
#pragma omp parallel for
    for (size_t send_i = 0; send_i < this->send_count; send_i++) {
        ghost_send_buffer[send_i] = local_and_ghost_partition_labels[this->sending_scatter_map[send_i]];
    }

    MPI_Barrier(this->comm);
    par::Mpi_Alltoallv_sparse(ghost_send_buffer.data(), this->send_counts.data(), this->send_counts_scanned.data(),
                              &local_and_ghost_partition_labels[this->local_count], this->ghost_counts.data(),
                              this->ghost_counts_scanned.data(), comm);

    MPI_Barrier(this->comm);

    // now calculating partition boundaries
    std::vector<uint32_t> local_partition_boundaries(procs_n, 0);

    for (graph_indexing_t local_vertex = 0; local_vertex < this->local_count; local_vertex++) {

        for (graph_indexing_t neighbor_i = this->local_xdj[local_vertex]; neighbor_i < this->local_xdj[local_vertex + 1];
             neighbor_i++) {
            auto neighbor = local_adjncy[neighbor_i];
            if (local_and_ghost_partition_labels[local_vertex] != local_and_ghost_partition_labels[neighbor]) {
                local_partition_boundaries[local_and_ghost_partition_labels[local_vertex]]++;
                break;
            }
        }
    }

    // collect results to 0 MPI proc

    if (!my_rank) {
        partition_sizes_out.resize(procs_n);
        std::fill(partition_sizes_out.begin(), partition_sizes_out.end(), 0);
        partition_boundaries_out.resize(procs_n);
        std::fill(partition_boundaries_out.begin(), partition_boundaries_out.end(), 0);
    }

    MPI_Reduce(local_partition_sizes.data(), partition_sizes_out.data(), procs_n, MPI_UINT32_T, MPI_SUM, 0, this->comm);
    MPI_Barrier(this->comm);
    MPI_Reduce(local_partition_boundaries.data(), partition_boundaries_out.data(), procs_n, MPI_UINT32_T, MPI_SUM, 0,
               this->comm);
}