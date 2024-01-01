/******************************************************************************
 * FILE: mpi_hello.c
 * DESCRIPTION:
 *   MPI tutorial example code: Simple hello world program
 * AUTHOR: Blaise Barney
 * LAST REVISED: 03/05/10
 ******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <numeric>

#include "graph.hpp"
#include "util.hpp"

#define MASTER 0
#define N 1000

#define BFS_SEED 500000

// #define SYNC_ROUNDS 3


Graph GetMyGraph(int my_rank, int num_tasks, int graph_size){
    Graph my_graph;

    for (vertex_t v = 1; v <= static_cast<vertex_t>(graph_size); v++)
    {
        if (IsMyVertex(v,my_rank,num_tasks,graph_size))
        {
            my_graph.AddVertex(v);
        }
        
    }

    auto vertices = my_graph.GetVertices();

    for (vertex_t v : vertices)
    {
        if (IsValidVertex(v - N,graph_size) && IsMyVertex(v-N,my_rank,num_tasks,graph_size))
        {
            my_graph.AddEdge(v, v - N);
        }
        if (IsValidVertex(v + N,graph_size) && IsMyVertex(v+N,my_rank,num_tasks,graph_size))
        {
            my_graph.AddEdge(v, v + N);
        }
        if (v % N)
        {
            if (IsValidVertex(v + 1,graph_size) && IsMyVertex(v+1,my_rank,num_tasks,graph_size))
            {
                my_graph.AddEdge(v, v + 1);
            }
        }
        if ((v - 1) % N)
        {
            if (IsValidVertex(v - 1,graph_size) && IsMyVertex(v-1,my_rank,num_tasks,graph_size))
            {
                my_graph.AddEdge(v, v - 1);
            }
        }
    }

    return my_graph;
    
    

}

int main(int argc, char *argv[])
{
    int numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
    MPI_Get_processor_name(hostname, &len);

    int graph_size = N*N;

    // if (taskid == 3)
    {
        // std::cout << "task: " << taskid <<" out of " << numtasks << "\nmy vertices\n";
        Graph my_graph = GetMyGraph(taskid,numtasks,graph_size);
        int my_boundary_size = GetMyBoundarySize(my_graph, taskid,numtasks,graph_size,N);
        std::cout << "task: " << taskid <<" out of " << numtasks << "\tboundary size =" << my_boundary_size << "\n";
        std::vector<int> boundary_sizes(numtasks);
        MPI_Allgather(&my_boundary_size, 1, MPI_INT, &boundary_sizes[0], 1, MPI_INT, MPI_COMM_WORLD);
        // std::cout << "task: " << taskid << " received boundary sizes " << VectorToString(boundary_sizes);

        int total_receiving_boundary_size =  0;
        std::vector<int> receive_buffer_displacements;

        for (auto b_size : boundary_sizes)
        {
            receive_buffer_displacements.push_back(total_receiving_boundary_size);
            total_receiving_boundary_size+=b_size;
            
        }
        

        auto my_boundary_vertices = GetMyBoundaryVertices(my_graph, taskid,numtasks,graph_size,N);

        // std::cout << "task: " << taskid << " my boundary vertices " << VectorToString(my_boundary_vertices);
        
        std::vector<vertex_t> received_all_boundary_vertices_buffer(total_receiving_boundary_size);
        MPI_Allgatherv(&my_boundary_vertices[0], my_boundary_size, MPI_UNSIGNED_LONG, &received_all_boundary_vertices_buffer[0], &boundary_sizes[0], &receive_buffer_displacements[0], MPI_UNSIGNED_LONG, MPI_COMM_WORLD);


        // std::cout << "task: " << taskid << " received boundary vertices " << VectorToString<vertex_t>(received_all_boundary_vertices_buffer);

        std::vector<std::vector<int>> my_ghost_vertices_relative_indices =  GetMyGhostVerticesRelativeIndices(my_graph,received_all_boundary_vertices_buffer,boundary_sizes,taskid,numtasks,graph_size, N);

        int ghost_vertex_count = 0;
        std::vector<int> ghost_vertex_count_per_proc;
        std::vector<int> ghost_vertex_count_per_proc_prefix_sum;

        std::vector<vertex_t> my_ghost_vertices;

        // if (taskid == 1){
        
        // std::cout << "task: " << taskid << " my_ghost_vertices_relative_indices = \n";

        int offset = 0;
        int tmp_i = 0;

        std::vector<int> my_ghost_vertices_relative_indices_flat;

        
        for (auto indices : my_ghost_vertices_relative_indices)
        {
            ghost_vertex_count_per_proc_prefix_sum.push_back(ghost_vertex_count);
            ghost_vertex_count+=indices.size();
            ghost_vertex_count_per_proc.push_back(indices.size());

            for (auto index : indices)
            {
                int vertex = received_all_boundary_vertices_buffer[offset+index];
                my_graph.AddVertex(static_cast<vertex_t>(vertex));
                my_ghost_vertices.push_back(vertex);
                my_ghost_vertices_relative_indices_flat.push_back(index);
                
            }

            offset += boundary_sizes[tmp_i];
            tmp_i++;        

            // if (indices.size()==0)
            // {
            //     std::cout << "[]\n";
            // }
            
            // std::cout << VectorToString<int>(indices);
        }
        // std::cout << "task: " << taskid << " my_ghost_vertices = " << VectorToString<vertex_t>(my_ghost_vertices);

        for (auto ghost_v : my_ghost_vertices)
        {
            AddEdgesToGhostVertex(my_graph,ghost_v,taskid,numtasks,graph_size, N);
        }
        
        // std::cout << "task: " << taskid << " my_gaph_with_ghost = \n";
        // my_graph.Print();
        // }
        // std::cout << "task: " << taskid << " ghost_vertex_count_per_proc" << VectorToString(ghost_vertex_count_per_proc);
        std::vector<int> other_ghosts_in_my_graph_counts(numtasks);

        MPI_Alltoall(&ghost_vertex_count_per_proc[0],1,MPI_INT,&other_ghosts_in_my_graph_counts[0],1,MPI_INT,MPI_COMM_WORLD);

        // std::cout << "task: " << taskid << " other_ghosts_in_my_graph_counts " << VectorToString(other_ghosts_in_my_graph_counts);

        int other_ghosts_in_my_graph_total_count = 0;
        std::vector<int> other_ghosts_in_my_graph_indices_receive_buffer_displacements;

        for (auto count : other_ghosts_in_my_graph_counts)
        {
            other_ghosts_in_my_graph_indices_receive_buffer_displacements.push_back(other_ghosts_in_my_graph_total_count);
            other_ghosts_in_my_graph_total_count+=count;
        }

        std::vector<int> other_ghosts_in_my_graph_indices(other_ghosts_in_my_graph_total_count);

        MPI_Alltoallv(&my_ghost_vertices_relative_indices_flat[0],&ghost_vertex_count_per_proc[0],&ghost_vertex_count_per_proc_prefix_sum[0],MPI_INT,
            &other_ghosts_in_my_graph_indices[0],&other_ghosts_in_my_graph_counts[0],&other_ghosts_in_my_graph_indices_receive_buffer_displacements[0],MPI_INT,MPI_COMM_WORLD);


        std::vector<vertex_t> other_ghosts_in_my_graph;

        for (auto other_ghost_v_idx : other_ghosts_in_my_graph_indices)
        {
            other_ghosts_in_my_graph.push_back(my_boundary_vertices[other_ghost_v_idx]);
            
        }

        // std::cout << "task: " << taskid << " other_ghosts_in_my_graph " << VectorToString(other_ghosts_in_my_graph);

        std::vector<int> ghost_send_displacements(other_ghosts_in_my_graph_indices_receive_buffer_displacements);
        std::vector<vertex_t> ghost_send_vertices(other_ghosts_in_my_graph);

        if (my_graph.FindVertex(BFS_SEED))
        {
            my_graph.InitSingleBFS(BFS_SEED);
        }

        std::vector<unsigned long> my_ghost_received_values(my_ghost_vertices.size());
        std::vector<unsigned long> other_ghost_send_values(other_ghosts_in_my_graph_total_count);

        bool global_is_not_stable = true;



        // for (int sync_i = 0; sync_i < SYNC_ROUNDS; sync_i++)
        int sync_steps = 0;
        while(global_is_not_stable)
        {
            global_is_not_stable = false;
            my_graph.RunBFSToStable();


            for (int other_ghost_i=0; other_ghost_i < other_ghosts_in_my_graph_total_count; other_ghost_i++)
            {
                other_ghost_send_values[other_ghost_i] = my_graph.GETBFSValue(other_ghosts_in_my_graph[other_ghost_i]);
            }
            
            MPI_Alltoallv(&other_ghost_send_values[0],&other_ghosts_in_my_graph_counts[0],&ghost_send_displacements[0],MPI_UNSIGNED_LONG,
                &my_ghost_received_values[0],&ghost_vertex_count_per_proc[0],&ghost_vertex_count_per_proc_prefix_sum[0],MPI_UNSIGNED_LONG,MPI_COMM_WORLD);

            // if (taskid==1)
            // {
            //     std::cout << "task: " << taskid << " my_ghost_received_values " <<VectorToString(my_ghost_received_values);
            // }
            
            bool local_is_not_stable = false;

            std::vector<unsigned long> & my_bfs_state = my_graph.GetBFSState();
            for (int recev_i = 0; recev_i < my_ghost_received_values.size(); recev_i++)
            {
                int bfs_index = recev_i + static_cast<int>(my_graph.GetSize()) - my_ghost_received_values.size();
                if (my_bfs_state[bfs_index] > my_ghost_received_values[recev_i])
                {
                    my_bfs_state[bfs_index] = my_ghost_received_values[recev_i];
                    local_is_not_stable = true;
                }
            }

            MPI_Allreduce(&local_is_not_stable, &global_is_not_stable, 1, MPI_CXX_BOOL, MPI_LOR, MPI_COMM_WORLD);
            sync_steps++;

        }   

        // if (taskid==3)
        // {
        //     my_graph.PrintSingleBFS();
        // }

        std::cout << "sync steps = " << sync_steps << "\n";
        

        

        
    }
    



    // my_graph.Print();
    // std::cout << "\n====BFS===\n";
    // my_graph.InitSingleBFS(13);
    // my_graph.RunBFSToStable();
    // auto bfs = my_graph.GetBFSState();
    // my_graph.PrintSingleBFS();
    MPI_Finalize();
}
