#include "bfs-com.hpp"
#include <numeric>


/*
Assumes the last element(s) in bfs_state are ghost layer nodes
*/
bool UpdateGhost(std::vector<unsigned long>& bfs_state, unsigned long graph_part_size, std::vector<int>& receive_counts, std::vector<std::vector<unsigned long>>& send_indices_per_proc, MPI_Comm& com){
    // int total_receive_size = std::accumulate(receive_counts.begin(),receive_counts.end(),0);
    int total_receive_size = 0;
    std::vector<int> receive_displacements;
    for (int receive_count: receive_counts)
    {
        receive_displacements.push_back(total_receive_size);
        total_receive_size+=receive_count;
    }
    auto receive_buf = std::vector<unsigned long>(total_receive_size);


    std::vector<int> send_counts;
    std::vector<int> send_displacements;
    std::vector<unsigned long> send_buf;
    int total_send_size = 0;

    for (std::vector<unsigned long> send_indices : send_indices_per_proc)
    {
        send_counts.push_back(send_indices.size());
        for (unsigned long send_i : send_indices)
        {
            send_buf.push_back(bfs_state[send_i]);      // packing for sending
        }
        
        send_displacements.push_back(total_send_size);
        total_send_size+=send_indices.size();

    }

    MPI_Alltoallv(&send_buf[0],&send_counts[0],&send_displacements[0],MPI_UNSIGNED_LONG,&receive_buf[0],&receive_counts[0],&receive_counts[0],MPI_UNSIGNED_LONG,com);

    //updating bfs state by checking whether the received values are better
    bool is_stable = true;
    for (int recev_i = 0; recev_i < total_receive_size; recev_i++)
    {
        int bfs_index = recev_i + static_cast<int>(graph_part_size) - total_receive_size;
        if (bfs_state[bfs_index] > receive_buf[recev_i])
        {
            bfs_state[bfs_index] = receive_buf[recev_i];
            is_stable = false;      
        }
        
    }

    return is_stable;
    

}


//