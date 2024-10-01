
#include "dist-graph.hpp"
#include "../mesh-util/mesh-util.hpp"
#include "../usort/ompUtils.h"
#include "../usort/parUtils.h"

#include "mpi.h"
#include "../metis-util/metis-util.hpp"
#include "../scotch-util/scotch-util.hpp"


#include <chrono>
#include <numeric>


// Overloading the << operator for BFSValue
std::ostream& operator<<(std::ostream& os, const BFSValue& obj) {
    os << "(" << obj.distance << "," << obj.label  << ")";
    // os << obj.global_idx;
    
    return os;
}

// Overloading the << operator for GhostBFSValue
std::ostream& operator<<(std::ostream& os, const GhostBFSValue& obj) {
    os << "(" << obj.distance << "," << obj.label << ", [" <<obj.offset  << "])";
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
    this->global_count = proc_element_counts_scanned[procs_n-1] + proc_element_counts[procs_n-1];
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

    for (int proc_i = 0; proc_i < procs_n; proc_i++)
    {
        if (this->ghost_counts[proc_i] > 0)
        {
            this->ghost_procs.push_back(proc_i);
        }
        
    }
    



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
        

        for (size_t boundary_edge_i =1; boundary_edge_i <boundary_connectivity_cpy.size(); boundary_edge_i++)
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

        for (size_t boundary_edge_i =1; boundary_edge_i <boundary_connectivity_cpy.size(); boundary_edge_i++)
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
        for (int proc_i = 0; proc_i < procs_n; proc_i++)
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

        for (int proc_i = 0; proc_i < procs_n; proc_i++)
        {
            if (this->send_counts[proc_i] > 0)
            {
                this->send_procs.push_back(proc_i);
            }
            
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
    std::vector<BFSValue> bfs_vector(this->local_xdj.size()-1);

    BFSValue bfs_init_value = {.label = DIST_GRAPH_BFS_NO_LABEL, .distance =  DIST_GRAPH_BFS_INFINITY};

    int BFS_stop_guess = std::min(procs_n, 7);      // TODO: this bound could be improved
    const float partition_size_imbalance_additive_factor = 0.2;       // BFS will try to keep the paritition sizes in range [ideal_size-factor, ideal_size+factor]
    assert( 0<= partition_size_imbalance_additive_factor && 1 > partition_size_imbalance_additive_factor);
    // // TODO: openmp can be used to populate
    // std::fill(bfs_vector.begin(),bfs_vector.end(), init_value);

    /**
     * using sfc seeds
     * elements are already ordered to morton SFC
     * get the 'middle' local element as seed
     * TODO: oversampling can be implemented here
    */
    graph_indexing_t seed = static_cast<graph_indexing_t>(this->local_count/2);

    // bfs_vector[seed].distance = 0;
    // bfs_vector[seed].label = my_rank;



    bool is_not_stable_global = true;      // global BFS stability
    int round_counter = 0;
    MPI_Barrier(this->comm);
    auto com_duration = std::chrono::microseconds(0);

    auto start = std::chrono::high_resolution_clock::now();



    //special first iteration using standard BFS frontier method with a queue
    this->RunFirstBFSIteration(bfs_vector, seed, my_rank);

    // when we receive ghost updates, keep track of the update with minimum value
    // then, in the next inner BFS, vertices can be filtered based on this value
    bfs_distance_t ghost_min_update = 0; 

    int refinement_rounds = 0;   
    int refinement_rounds_stop = 1;    

    if(procs_n > 1)
    {
        
        std::vector<BFSValue> ghost_send_buffer(this->send_count);
        std::vector<BFSValue> ghost_send_buffer_prev(this->send_count, bfs_init_value);

        std::vector<BFSValue> ghost_recv_buffer(this->ghost_count, bfs_init_value);

        std::vector<bool> ghost_updated(this->ghost_count,false);

        this->ghost_count_requests = new MPI_Request[this->ghost_procs.size() + this->send_procs.size()];
        this->ghost_count_statuses = new MPI_Status[this->ghost_procs.size() + this->send_procs.size()];
        this->updated_only_recv_counts.resize(procs_n);

        
        
        while (is_not_stable_global)
        {
            // if (!my_rank)
            // {
            //     print_log("BFS round: ", ++round_counter);
            // }
            round_counter++;
            
            is_not_stable_global = false;
            bool is_not_stable_local = true;

            if (round_counter > 1)
            {
                auto start_ = std::chrono::high_resolution_clock::now();
                this->StartReceivingUpdatedOnlyGhostCounts();
                auto end_ = std::chrono::high_resolution_clock::now();
                com_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_);
            }
            
            if (round_counter > 1)      // first round is already handled by the "special first iteration"
            {   
                is_not_stable_local = this->RunLocalMultiBFSToStable2(bfs_vector, ghost_updated, ghost_min_update);

            }
            
            // print_log("[", my_rank, "]: BFS iteration done");

            #pragma omp parallel for
            for (size_t send_i = 0; send_i < this->send_count; send_i++)
            {
                ghost_send_buffer[send_i] = bfs_vector[this->sending_scatter_map[send_i]];
            }
            auto com_start = std::chrono::high_resolution_clock::now();
            
            
            /**
             * ghost exchange
            */
            if (round_counter > 1)          // after the first round, exchanging updated only ghosts is better for communication
            {
                this->ExchangeUpdatedOnlyBFSGhost(ghost_send_buffer,ghost_send_buffer_prev,ghost_recv_buffer);
                
            }else
            {
                // MPI_Barrier(this->comm);
                par::Mpi_Alltoallv_sparse(ghost_send_buffer.data(), this->send_counts.data(), this->send_counts_scanned.data(), 
                        ghost_recv_buffer.data(), this->ghost_counts.data(), this->ghost_counts_scanned.data(), comm);

                // MPI_Barrier(this->comm);
            }    

            auto com_end = std::chrono::high_resolution_clock::now();
            com_duration += std::chrono::duration_cast<std::chrono::microseconds>(com_end - com_start);

            std::copy(ghost_send_buffer.begin(), ghost_send_buffer.end(), ghost_send_buffer_prev.begin());      // for the next iteration

            /**
             * ghost update using received values
            */
            bool ghost_any_is_not_stable = false;

            ghost_min_update = DIST_GRAPH_BFS_INFINITY;
            std::fill(ghost_updated.begin(), ghost_updated.end(), false);
            for (size_t recv_i = 0; recv_i < this->ghost_count; recv_i++)
            {
                auto offset = this->local_count;        // ghost elements are in the last section of the vector, in sorted order
                if (bfs_vector[offset+recv_i].distance > ghost_recv_buffer[recv_i].distance)
                {
                    bfs_vector[offset+recv_i].distance = ghost_recv_buffer[recv_i].distance;
                    bfs_vector[offset+recv_i].label = ghost_recv_buffer[recv_i].label;

                    ghost_any_is_not_stable = true;
                    ghost_min_update = std::min(ghost_min_update, ghost_recv_buffer[recv_i].distance);
                    ghost_updated[recv_i] = true;
                }            
            }
            is_not_stable_local = is_not_stable_local || ghost_any_is_not_stable;

            if(round_counter >= BFS_stop_guess){
                auto start_ = std::chrono::high_resolution_clock::now();
                MPI_Allreduce(&is_not_stable_local,&is_not_stable_global,1,MPI_CXX_BOOL,MPI_LOR,this->comm);
                auto end_ = std::chrono::high_resolution_clock::now();
                com_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_);
            } else {
                is_not_stable_global = true;
            }
            // print_log("[", my_rank, "]: BFS vector", VectorToString(bfs_vector));


            // partition refinement after reaching global stable status
            if (!is_not_stable_global && refinement_rounds < refinement_rounds_stop)
            {
                refinement_rounds++;

               
                std::vector<uint32_t> local_partition_sizes(procs_n, 0);
                for (size_t local_i = 0; local_i < this->local_count; local_i++) {
                    local_partition_sizes[bfs_vector[local_i].label]++;
                }
                std::vector<uint32_t> global_partition_sizes(procs_n, 0);
                
                {
                    auto start_ = std::chrono::high_resolution_clock::now();
                    MPI_Allreduce(local_partition_sizes.data(),global_partition_sizes.data(),procs_n,MPI_UINT32_T,MPI_SUM,this->comm);
                    auto end_ = std::chrono::high_resolution_clock::now();
                    com_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_);
                }
                

                uint32_t ideal_partition_size = this->global_count / procs_n;
                
                // auto partition_size_cutoff_max = static_cast<uint32_t>(ideal_partition_size*(1.0 + partition_size_imbalance_additive_factor));   
                // auto partition_size_cutoff_min = static_cast<uint32_t>(ideal_partition_size*(1.0 - partition_size_imbalance_additive_factor)); 

                auto cutoff_max_1 =   static_cast<uint32_t>(ideal_partition_size*1.2);
                auto cutoff_max_2 =   static_cast<uint32_t>(ideal_partition_size*1.4);
                auto cutoff_max_3 =   static_cast<uint32_t>(ideal_partition_size*1.6);
                auto cutoff_max_4 =   static_cast<uint32_t>(ideal_partition_size*1.8);
                auto cutoff_max_5 =   static_cast<uint32_t>(ideal_partition_size*2.0);
                auto cutoff_max_6 =   static_cast<uint32_t>(ideal_partition_size*2.2);


                auto cutoff_min_1 =   static_cast<uint32_t>(ideal_partition_size*0.8);
                auto cutoff_min_2 =   static_cast<uint32_t>(ideal_partition_size*0.6);
                auto cutoff_min_3 =   static_cast<uint32_t>(ideal_partition_size*0.4);
                auto cutoff_min_4 =   static_cast<uint32_t>(ideal_partition_size*0.2);
                

                std::vector<int> part_adjust(procs_n,0);

                for (int proc_i = 0; proc_i < procs_n; proc_i++)
                {
                    if (global_partition_sizes[proc_i] >= cutoff_max_6)
                    {
                        part_adjust[proc_i] = 6;
                    }else if (global_partition_sizes[proc_i] >= cutoff_max_5)
                    {
                        part_adjust[proc_i] = 5;
                    }else if (global_partition_sizes[proc_i] >= cutoff_max_4)
                    {
                        part_adjust[proc_i] = 4;
                    }else if (global_partition_sizes[proc_i] >= cutoff_max_3)
                    {
                        part_adjust[proc_i] = 3;
                    }else if (global_partition_sizes[proc_i] >= cutoff_max_2)
                    {
                        part_adjust[proc_i] = 2;
                    }else if (global_partition_sizes[proc_i] >= cutoff_max_1)
                    {
                        part_adjust[proc_i] = 1;
                    }else if (global_partition_sizes[proc_i] <= cutoff_min_4)
                    {
                        part_adjust[proc_i] = -4;
                    }else if (global_partition_sizes[proc_i] <= cutoff_min_3)
                    {
                        part_adjust[proc_i] = -3;
                    }else if (global_partition_sizes[proc_i] <= cutoff_min_2)
                    {
                        part_adjust[proc_i] = -2;
                    }else if (global_partition_sizes[proc_i] <= cutoff_min_1)
                    {
                        part_adjust[proc_i] = -1;
                    }
                    
                    
                }

                // if(! my_rank) print_log(VectorToString(global_partition_sizes));
                // if(! my_rank) print_log(VectorToString(part_adjust));

                for (size_t vec_i = 0; vec_i < bfs_vector.size(); vec_i++)
                {
                    auto adjust_val = part_adjust[bfs_vector[vec_i].label];
                    if (adjust_val > 0)
                    {
                        bfs_vector[vec_i].distance += adjust_val;
                    }else if (adjust_val < 0)
                    {
                        bfs_vector[vec_i].distance = bfs_vector[vec_i].distance >= std::abs(adjust_val) ? bfs_vector[vec_i].distance - std::abs(adjust_val) : 0;   // unsigned underflow prevention
                    }                    
                    
                }
                round_counter++;

                this->RunLocalMultiBFSToStable(bfs_vector);
                
                // one ghost exchange
                for (size_t send_i = 0; send_i < this->send_count; send_i++)
                {
                    ghost_send_buffer[send_i] = bfs_vector[this->sending_scatter_map[send_i]];
                }

                {
                    auto start_ = std::chrono::high_resolution_clock::now();
                    // force all ghosts to be exchanged because we want to reset send buffers
                    par::Mpi_Alltoallv_sparse(ghost_send_buffer.data(), this->send_counts.data(), this->send_counts_scanned.data(), 
                            ghost_recv_buffer.data(), this->ghost_counts.data(), this->ghost_counts_scanned.data(), comm);
                    auto end_ = std::chrono::high_resolution_clock::now();
                    com_duration += std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_);
                }


                std::copy(ghost_send_buffer.begin(), ghost_send_buffer.end(), ghost_send_buffer_prev.begin());      // for the next iteration             

                ghost_min_update = DIST_GRAPH_BFS_INFINITY;
                std::fill(ghost_updated.begin(), ghost_updated.end(), false);
                for (size_t recv_i = 0; recv_i < this->ghost_count; recv_i++)
                {
                    auto offset = this->local_count;        // ghost elements are in the last section of the vector, in sorted order
                    if (bfs_vector[offset+recv_i].distance > ghost_recv_buffer[recv_i].distance)
                    {
                        bfs_vector[offset+recv_i].distance = ghost_recv_buffer[recv_i].distance;
                        bfs_vector[offset+recv_i].label = ghost_recv_buffer[recv_i].label;

                        ghost_any_is_not_stable = true;
                        ghost_min_update = std::min(ghost_min_update, ghost_recv_buffer[recv_i].distance);
                        ghost_updated[recv_i] = true;
                    }            
                }


                is_not_stable_global = true;
            }
            

            
        } 

        delete[] this->ghost_count_requests;
        delete[] this->ghost_count_statuses;

    }

    // print_log("[", my_rank, "]: BFS done");
    // print_log("[", my_rank, "]: BFS vector", VectorToString(bfs_vector));

    MPI_Barrier(this->comm);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (!my_rank)
    {
        print_log("BFS sync rounds: ", round_counter);
        print_log("BFS comm  time:\t\t\t", com_duration.count(), " us");
        print_log("BFS total time:\t\t\t", duration.count(), " us");
    }

    partition_labels_out.resize(this->local_count);


    #pragma omp parallel for
    for (size_t local_i = 0; local_i < this->local_count; local_i++)
    {
        partition_labels_out[local_i] = bfs_vector[local_i].label;
    }

    return {.return_code = 0, .time_us = static_cast<int>(duration.count())};

}

/**
 * 
 * The first BFS iteration in the local partition only has one frontier.
 * This function implements the BFS in the standard way using frontier in a queue
*/
void DistGraph::RunFirstBFSIteration(std::vector<BFSValue>& bfs_vector, graph_indexing_t seed, bfs_label_t label) {
    
    BFSValue init_value = {.label = DIST_GRAPH_BFS_NO_LABEL, .distance =  DIST_GRAPH_BFS_INFINITY};

    // TODO: openmp can be used to populate
    std::fill(bfs_vector.begin(),bfs_vector.end(), init_value);

    //mark the seed
    bfs_vector[seed].label = label;
    bfs_vector[seed].distance = 0;

    std::vector<graph_indexing_t> frontier_buffer(bfs_vector.size());
    graph_indexing_t curr_frontier_start = 0;
    graph_indexing_t curr_frontier_size = 1;

    frontier_buffer[curr_frontier_start] = seed; // first frontier has only the seed

    bfs_distance_t curr_distance = 1;

    while (curr_frontier_size != 0) {

        graph_indexing_t next_frontier_size = 0;
        graph_indexing_t next_frontier_slot = curr_frontier_start + curr_frontier_size;
        for (graph_indexing_t frontier_i = curr_frontier_start; frontier_i < curr_frontier_start + curr_frontier_size;
             frontier_i++) {
            graph_indexing_t frontier_vertex = frontier_buffer[frontier_i];

            // looping neighbors
            for (graph_indexing_t neighbor_i = this->local_xdj[frontier_vertex];
                 neighbor_i < this->local_xdj[frontier_vertex + 1]; neighbor_i++) 
            {
                auto neighbor = local_adjncy[neighbor_i];

                if(neighbor >= this->local_count) continue;      // we dont update the ghost vertices

                if (bfs_vector[neighbor].label == DIST_GRAPH_BFS_NO_LABEL) {
                    // neighbor is not visited, visit it now
                    bfs_vector[neighbor].label = label;
                    bfs_vector[neighbor].distance = curr_distance;

                    // add to next frontier
                    frontier_buffer[next_frontier_slot++] = neighbor;
                    next_frontier_size++;
                }
            }
        }
        curr_distance++;
        curr_frontier_start+= curr_frontier_size;
        curr_frontier_size = next_frontier_size;
    }
}

/**
 * to be used as an optimization step
 * this function calculates the minimum distance from last updated ghost vertices
 * afterwards this distances can be used for the filtering citeria calculation in subsequent BFS rounds
 * 
*/
void DistGraph::CalculateDistanceFromUpdatedGhosts(std::vector<bool>& ghost_updated, std::vector<graph_indexing_t>& frontier_buffer,
                std::vector<bfs_distance_t>& distances_out){
    // TODO: this function assumes at least one ghost vertex is present

    // distances_out.resize(this->local_count + this->ghost_count);
    std::fill(distances_out.begin(), distances_out.end(), DIST_GRAPH_BFS_INFINITY);
    // std::vector<graph_indexing_t> frontier_buffer(distances_out.size());            // to hold BFS queue
    graph_indexing_t curr_frontier_size = 0;
    for (graph_indexing_t ghost_i = 0; ghost_i < this->ghost_count; ghost_i++)
    {
        if (ghost_updated[ghost_i])
        {
            // ghost vertices are in the last part in the distance vector, so we add this->local_count offset
            distances_out[this->local_count + ghost_i] = 0;      

            frontier_buffer[curr_frontier_size++] = this->local_count + ghost_i;     // initial frontier is all the updated ghost vertices  
        }
        

    }
    
    
    graph_indexing_t curr_frontier_start = 0;


    bfs_distance_t curr_distance = 1;

    while (curr_frontier_size != 0) {

        graph_indexing_t next_frontier_size = 0;
        graph_indexing_t next_frontier_slot = curr_frontier_start + curr_frontier_size;
        for (graph_indexing_t frontier_i = curr_frontier_start; frontier_i < curr_frontier_start + curr_frontier_size;
             frontier_i++) {
            graph_indexing_t frontier_vertex = frontier_buffer[frontier_i];

            // looping neighbors
            for (graph_indexing_t neighbor_i = this->local_xdj[frontier_vertex];
                 neighbor_i < this->local_xdj[frontier_vertex + 1]; neighbor_i++) 
            {
                auto neighbor = local_adjncy[neighbor_i];

                if (distances_out[neighbor] == DIST_GRAPH_BFS_INFINITY) {
                    // neighbor is not reached by frontier, visit it now
                    distances_out[neighbor] = curr_distance;

                    // add to next frontier
                    frontier_buffer[next_frontier_slot++] = neighbor;
                    next_frontier_size++;
                }
            }
        }
        curr_distance++;
        curr_frontier_start+= curr_frontier_size;
        curr_frontier_size = next_frontier_size;
    }

}

bool DistGraph::RunLocalMultiBFSToStable(std::vector<BFSValue>& bfs_vector){
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);


    bool changed = false;       // to detect if the BFS incremented

    std::vector<graph_indexing_t> frontier(bfs_vector.size());
    std::iota(frontier.begin(), frontier.end(), 0);


    std::vector<bool> added_to_next_updated(bfs_vector.size());

    while (frontier.size() > 0)
    {
        std::fill(added_to_next_updated.begin(), added_to_next_updated.end(), false);
        std::vector<std::pair<graph_indexing_t, BFSValue>> next_updated;                 // next vertices for potential updates
        for (auto& frontier_vertex : frontier) {
            // looping neighbors, add to next_updated
            for (graph_indexing_t neighbor_i = this->local_xdj[frontier_vertex];
                    neighbor_i < this->local_xdj[frontier_vertex + 1]; neighbor_i++) 
            {
                auto neighbor = local_adjncy[neighbor_i];
                if(neighbor >= this->local_count) continue;      // we dont update the ghost vertices

                if (!added_to_next_updated[neighbor])
                {
                    next_updated.push_back({neighbor, bfs_vector[neighbor]});
                    added_to_next_updated[neighbor] = true;
                }
                
            }
        }

        frontier.clear();
    
        // now considering only the nighborhoods of next potential updates
        for (auto& it: next_updated) {
            auto& vertex = it.first;
            auto& curr_bfs_val = it.second;
            bool value_changed = false;

            for (graph_indexing_t neighbor_i = this->local_xdj[vertex];
                 neighbor_i < this->local_xdj[vertex + 1]; neighbor_i++) 
            {
                auto neighbor = local_adjncy[neighbor_i];

                if (bfs_vector[neighbor].label != DIST_GRAPH_BFS_NO_LABEL &&
                    curr_bfs_val.distance > (bfs_vector[neighbor].distance + 1)) {

                    curr_bfs_val.label = bfs_vector[neighbor].label;
                    curr_bfs_val.distance = bfs_vector[neighbor].distance + 1;

                    if(! value_changed) frontier.push_back(vertex);
                    value_changed = true;
                    changed = true;
                    
                }
            }
        }

        // writing updated values to BFS vector
        for (auto& it: next_updated){
            bfs_vector[it.first].label = it.second.label;
            bfs_vector[it.first].distance = it.second.distance;

        }
    }

    return changed;
}


/**
 * maintains a combined bfs frontier and processes only the neighborhood of the frontier
*/
bool DistGraph::RunLocalMultiBFSToStable2(std::vector<BFSValue>& bfs_vector, std::vector<bool>& ghost_updated, 
    bfs_distance_t ghost_min_update){

    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);

    if(ghost_min_update == DIST_GRAPH_BFS_INFINITY){    // no change in ghost, dont have to run local BFS
        return false;
    }

    bool changed = false;       // to detect if the BFS incremented

    std::vector<graph_indexing_t> frontier;

    for (graph_indexing_t ghost_i = 0; ghost_i < this->ghost_count; ghost_i++)
    {
        if (ghost_updated[ghost_i])
        {
            frontier.push_back(this->local_count + ghost_i);
        }
        
    }

    std::vector<bool> added_to_next_updated(bfs_vector.size());

    while (frontier.size() > 0)
    {
        std::fill(added_to_next_updated.begin(), added_to_next_updated.end(), false);
        std::vector<std::pair<graph_indexing_t, BFSValue>> next_updated;                 // next vertices for potential updates
        for (auto& frontier_vertex : frontier) {
            // looping neighbors, add to next_updated
            for (graph_indexing_t neighbor_i = this->local_xdj[frontier_vertex];
                    neighbor_i < this->local_xdj[frontier_vertex + 1]; neighbor_i++) 
            {
                auto neighbor = local_adjncy[neighbor_i];
                if(neighbor >= this->local_count) continue;      // we dont update the ghost vertices
                if (!added_to_next_updated[neighbor])
                {
                    next_updated.push_back({neighbor, bfs_vector[neighbor]});
                    added_to_next_updated[neighbor] = true;
                }
                
            }
        }

        frontier.clear();
    
        // now considering only the nighborhoods of next potential updates
        for (auto& it: next_updated) {
            auto& vertex = it.first;
            auto& curr_bfs_val = it.second;
            bool value_changed = false;

            for (graph_indexing_t neighbor_i = this->local_xdj[vertex];
                 neighbor_i < this->local_xdj[vertex + 1]; neighbor_i++) 
            {
                auto neighbor = local_adjncy[neighbor_i];

                if (bfs_vector[neighbor].label != DIST_GRAPH_BFS_NO_LABEL &&
                    curr_bfs_val.distance > (bfs_vector[neighbor].distance + 1)) {

                    curr_bfs_val.label = bfs_vector[neighbor].label;
                    curr_bfs_val.distance = bfs_vector[neighbor].distance + 1;

                    if(! value_changed) frontier.push_back(vertex);
                    value_changed = true;
                    changed = true;
                    
                }
            }
        }

        // writing updated values to BFS vector
        for (auto& it: next_updated){
            bfs_vector[it.first].label = it.second.label;
            bfs_vector[it.first].distance = it.second.distance;

        }
    }

    return changed;
    
    

     


}


/**
 * optimized ghost update routine for distributed BFS
 * checks the previous sent buffer and current sending buffer, and sends only updated values
 * updated is defined as: BFS distance got reduced
 * ghost_recv_buffer is updated only for the received values. others are unchanged
*/
void DistGraph::ExchangeUpdatedOnlyBFSGhost(std::vector<BFSValue>& ghost_send_buffer,
                                            std::vector<BFSValue>& ghost_send_buffer_prev,
                                            std::vector<BFSValue>& ghost_recv_buffer)
{
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);

    std::vector<int> updated_only_send_counts(procs_n,0);
    std::vector<int> updated_only_send_counts_scanned(procs_n,0);
    std::vector<GhostBFSValue> updated_only_send_buffer;


    for (int proc_i = 0; proc_i < procs_n; proc_i++)
    {
        for (int send_offset = 0; send_offset < this->send_counts[proc_i]; send_offset++)
        {
            int idx = send_offset + this->send_counts_scanned[proc_i];

            if (ghost_send_buffer[idx].distance < ghost_send_buffer_prev[idx].distance)             // updated only check
            {
                updated_only_send_counts[proc_i]++;
                updated_only_send_buffer.push_back({.label = ghost_send_buffer[idx].label,
                                                    .distance = ghost_send_buffer[idx].distance,
                                                    .offset = static_cast<graph_indexing_t>(send_offset)});
            }
        }
        
    }

    omp_par::scan(&updated_only_send_counts[0], &updated_only_send_counts_scanned[0], procs_n);

    std::vector<int> updated_only_recv_counts(procs_n,0);

    this->EndExchangingUpdatedOnlyGhostCounts(updated_only_send_counts,updated_only_recv_counts);

    // this->ExchangeUpdatedOnlyBFSCounts(updated_only_send_counts,updated_only_recv_counts);
    
    std::vector<int> updated_only_recv_counts_scanned(procs_n,0);
    omp_par::scan(&updated_only_recv_counts[0], &updated_only_recv_counts_scanned[0], procs_n);

    int total_updated_only_recv_count = updated_only_recv_counts_scanned[procs_n-1] + updated_only_recv_counts[procs_n-1];

    // if(!my_rank) print_log("[",my_rank,"] receving_now/total_ghost (%) = ",100 * static_cast<float>(total_updated_only_recv_count)/this->ghost_count, "%");

    std::vector<GhostBFSValue> updated_only_recv_buffer(total_updated_only_recv_count);


    // MPI_Barrier(this->comm);
    par::Mpi_Alltoallv_sparse(updated_only_send_buffer.data(), updated_only_send_counts.data(), updated_only_send_counts_scanned.data(), 
            updated_only_recv_buffer.data(), updated_only_recv_counts.data(), updated_only_recv_counts_scanned.data(), comm);
    // MPI_Barrier(this->comm);


    // now place the received values in correct places in the original receive buffer

    for (int proc_i = 0; proc_i < procs_n; proc_i++)
    {
        // if (updated_only_recv_counts[proc_i] == 0)
        // {
        //     continue;
        // }
        for (int received_i = 0; received_i < updated_only_recv_counts[proc_i]; received_i++)
        {
            int idx_in_changed_only_recv_buffer = updated_only_recv_counts_scanned[proc_i] + received_i;
            GhostBFSValue received_value = updated_only_recv_buffer[idx_in_changed_only_recv_buffer];

            int idx_in_original_recv_buffer = static_cast<int>(received_value.offset) + this->ghost_counts_scanned[proc_i];

            ghost_recv_buffer[idx_in_original_recv_buffer].label = received_value.label;
            ghost_recv_buffer[idx_in_original_recv_buffer].distance = received_value.distance;

        }
        

        
    }
    



    
}


/**
 * MPI routine to exchange the updated only counts before performing the actual updated only exchange.
 * Uses `this->ghost_procs` and `this->send_procs` to perform point to point comm
 * We need this because if instead we use `MPI_Alltoall`, it will be a global lock
 * and the subsequent `Mpi_Alltoallv_sparse` would not get the full benefit
*/
void DistGraph::ExchangeUpdatedOnlyBFSCounts(std::vector<int>& updated_only_send_counts, std::vector<int>& updated_only_recv_counts_out){
    
    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);

    // exchange my own count (this is not actually needed since there are no ghost or sending vertices to self)
    updated_only_recv_counts_out[my_rank] = updated_only_send_counts[my_rank];


    // now exchanging counts only within neighboring processes

    MPI_Request* requests = new MPI_Request[this->ghost_procs.size() + this->send_procs.size()];
    assert(requests);

    MPI_Status* statuses = new MPI_Status[this->ghost_procs.size() + this->send_procs.size()];
    assert(statuses);

    int mpi_idx = 0;

    for (auto& recv_proc : this->ghost_procs)
    {
        par::Mpi_Irecv(&(updated_only_recv_counts_out[recv_proc]), 1, recv_proc, this->GHOST_COUNT_EXCHANGE_TAG, this->comm, &(requests[mpi_idx]));
        mpi_idx++;
    }
    
    for (auto& send_proc : this->send_procs)
    {
        par::Mpi_Isend(&(updated_only_send_counts[send_proc]), 1, send_proc, this->GHOST_COUNT_EXCHANGE_TAG, this->comm, &(requests[mpi_idx]));
        mpi_idx++;
    }

    MPI_Waitall(this->ghost_procs.size() + this->send_procs.size(), requests, statuses);

    delete [] requests;
    delete [] statuses;
    

}


/**
 * this will issue the Mpi_Irecv requests. later we should wait on these requests.
*/
void DistGraph::StartReceivingUpdatedOnlyGhostCounts(){

    std::fill(this->updated_only_recv_counts.begin(), this->updated_only_recv_counts.end(), 0);

    int mpi_idx = 0;
    for (auto& recv_proc : this->ghost_procs)
    {
        par::Mpi_Irecv(&(this->updated_only_recv_counts[recv_proc]), 1, recv_proc, this->GHOST_COUNT_EXCHANGE_TAG, 
                this->comm, &(this->ghost_count_requests[mpi_idx]));
        mpi_idx++;
    }

}

/**
 * Counterpart (i.e. the matching end) procedure of StartReceivingUpdatedOnlyGhostCounts
*/
void DistGraph::EndExchangingUpdatedOnlyGhostCounts(std::vector<int>& updated_only_send_counts, std::vector<int>& updated_only_recv_counts_out){

    int procs_n, my_rank;
    MPI_Comm_size(this->comm, &procs_n);
    MPI_Comm_rank(this->comm, &my_rank);


    // since the first part of this->ghost_count_requests contains the recev requests (from ghost procs), we have to offset that
    int mpi_idx = this->ghost_procs.size();

    for (auto& send_proc : this->send_procs)
    {
        par::Mpi_Isend(&(updated_only_send_counts[send_proc]), 1, send_proc, this->GHOST_COUNT_EXCHANGE_TAG, 
                this->comm, &(this->ghost_count_requests[mpi_idx]));
        mpi_idx++;
    }
    // print_log("[", my_rank, "]: waiting on EndExchangingUpdatedOnlyGhostCounts");

    // auto start_ = std::chrono::high_resolution_clock::now();
    MPI_Waitall(this->ghost_procs.size() + this->send_procs.size(), this->ghost_count_requests, this->ghost_count_statuses);
    // print_log("[", my_rank, "]: done EndExchangingUpdatedOnlyGhostCounts");
    // auto end_ = std::chrono::high_resolution_clock::now();
    // if(!my_rank) print_log("count exchange waitall time: ", std::chrono::duration_cast<std::chrono::microseconds>(end_ - start_).count(), "us");



    std::copy(this->updated_only_recv_counts.begin(), updated_only_recv_counts.end(), updated_only_recv_counts_out.begin());

    // exchange my own count (this is not actually needed since there are no ghost or sending vertices to self)
    updated_only_recv_counts_out[my_rank] = updated_only_send_counts[my_rank];

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

PartitionStatus DistGraph::PartitionPtScotch(std::vector<uint16_t>& partition_labels_out) {
    int procs_n;
    MPI_Comm_size(this->comm, &procs_n);
    std::vector<uint64_t> dist_xadj(this->local_xdj.begin(), this->local_xdj.begin() + (this->local_count + 1));
    return GetPtScotchPartitions(this->vtx_dist, dist_xadj, this->dist_adjncy, this->local_count, this->global_count,
                                 procs_n, partition_labels_out, this-> comm);
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

    if(procs_n > 1){
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
    }

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