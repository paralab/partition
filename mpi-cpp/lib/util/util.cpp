#include "util.hpp"
#include <gmsh.h>

#include <stdexcept>
#include <cassert>

bool IsValidVertex(int x, int graph_size)
{
    return 1 <= x && x <= graph_size;
}

bool IsMyVertex(vertex_t v, int my_rank, int num_tasks, int graph_size){
    int per_proc = graph_size/num_tasks;
    int mod = graph_size%num_tasks;
    int lower, upper;
    if (my_rank < mod)
    {
        lower = my_rank * (per_proc+1) + 1;
        upper = lower + (per_proc + 1) - 1;
    }else
    {
        lower = mod* (per_proc+1) + (my_rank-mod)*per_proc + 1;
        upper = lower + per_proc - 1;
    }   

    return static_cast<vertex_t>(lower) <= v && v<= static_cast<vertex_t>(upper); 

}

bool IsMyGhostVertex(vertex_t v, int my_rank, int num_tasks, int graph_size, int N){
    bool is_a_ghost_vertex = false;
    if (IsValidVertex(v - N,graph_size) && IsMyVertex(v-N,my_rank,num_tasks,graph_size))
    {
        is_a_ghost_vertex = true;
    }
    else if (IsValidVertex(v + N,graph_size) && IsMyVertex(v+N,my_rank,num_tasks,graph_size))
    {
        is_a_ghost_vertex = true;
    }
    else if (v % N)
    {
        if (IsValidVertex(v + 1,graph_size) && IsMyVertex(v+1,my_rank,num_tasks,graph_size))
        {
            is_a_ghost_vertex = true;
        }
    }
    else if ((v - 1) % N)
    {
        if (IsValidVertex(v - 1,graph_size) && IsMyVertex(v-1,my_rank,num_tasks,graph_size))
        {
            is_a_ghost_vertex = true;
        }
    }

    return is_a_ghost_vertex;

}

int GetMyBoundarySize(Graph& graph, int my_rank, int num_tasks, int graph_size, int N){
    int boundary_size = 0;
    // graph.Print();

    auto vertices = graph.GetVertices();
    for (auto v : vertices)
    {
        if (IsValidVertex(v - N,graph_size) && !IsMyVertex(v-N,my_rank,num_tasks,graph_size))
        {
            boundary_size++;
        }
        else if (IsValidVertex(v + N,graph_size) && !IsMyVertex(v+N,my_rank,num_tasks,graph_size))
        {
            boundary_size++;
        }
        else if (v % N)
        {
            if (IsValidVertex(v + 1,graph_size) && !IsMyVertex(v+1,my_rank,num_tasks,graph_size))
            {
                boundary_size++;
            }
        }
        else if ((v - 1) % N)
        {
            if (IsValidVertex(v - 1,graph_size) && !IsMyVertex(v-1,my_rank,num_tasks,graph_size))
            {
                boundary_size++;
            }
        }
    }

    return boundary_size;
    
}


std::vector<vertex_t> GetMyBoundaryVertices(Graph& graph, int my_rank, int num_tasks, int graph_size, int N){
    std::vector<vertex_t> boundary_vertices;
    auto vertices = graph.GetVertices();
    for (auto v : vertices)
    {
        if (IsValidVertex(v - N,graph_size) && !IsMyVertex(v-N,my_rank,num_tasks,graph_size))
        {
            boundary_vertices.push_back(v);
        }
        else if (IsValidVertex(v + N,graph_size) && !IsMyVertex(v+N,my_rank,num_tasks,graph_size))
        {
            boundary_vertices.push_back(v);
        }
        else if (v % N)
        {
            if (IsValidVertex(v + 1,graph_size) && !IsMyVertex(v+1,my_rank,num_tasks,graph_size))
            {
                boundary_vertices.push_back(v);
            }
        }
        else if ((v - 1) % N)
        {
            if (IsValidVertex(v - 1,graph_size) && !IsMyVertex(v-1,my_rank,num_tasks,graph_size))
            {
                boundary_vertices.push_back(v);
            }
        }
    }


    return boundary_vertices;
}


std::vector<std::vector<int>> GetMyGhostVerticesRelativeIndices(Graph& graph, std::vector<vertex_t> external_vertices, std::vector<int> external_vertex_counts, int my_rank, int num_tasks, int graph_size, int N){
    std::vector<std::vector<int>> ghost_vertex_relatives_indices;
    int prefix_sum = 0;
    for (int p_i = 0; p_i < num_tasks; p_i++)
    {
        ghost_vertex_relatives_indices.push_back({});
        if (p_i != my_rank)
        {
            for (int i = 0; i < external_vertex_counts[p_i]; i++)
            {
                int vertex = external_vertices[prefix_sum + i];
                if (IsMyGhostVertex(vertex, my_rank, num_tasks, graph_size, N))
                {
                    ghost_vertex_relatives_indices[p_i].push_back(i);
                }
                
            }
        } else
        {
            ghost_vertex_relatives_indices[p_i].resize(0);
        }
        

        prefix_sum += external_vertex_counts[p_i];   

        
    }
    
    return ghost_vertex_relatives_indices;

    
}

void AddEdgesToGhostVertex(Graph& graph, vertex_t v, int my_rank, int num_tasks, int graph_size, int N){
    if (graph.FindVertex(v-N))
    {
        graph.AddEdge(v, v-N);
    }
    if (graph.FindVertex(v+N))
    {
        graph.AddEdge(v, v+N);
    }
    if (v % N)
    {
        if (graph.FindVertex(v+1))
        {
            graph.AddEdge(v, v+1);
        }
    }
    if ((v - 1) % N)
    {
        if (graph.FindVertex(v-1))
        {
            graph.AddEdge(v, v-1);
        }
    }
}



// template <typename T> std::string VectorToString(std::vector<T> vec){
//     std::ostringstream output;
//     for (auto element : vec)
//     {
//         output << element << " ";
//     }
//     output << "\n";
//     return output.str();
// }

void AssignPartitionLabelsInOrder(std::vector<uint64_t> &ordering, uint64_t count, uint64_t partition_count, std::vector<uint64_t> &labels_out){
    assert(labels_out.size() == count);
    uint64_t partition_size = count/partition_count;
    uint64_t large_partition_count = count%partition_count;

    for (uint64_t p_i = 0; p_i < partition_count; p_i++)
    {
        uint64_t size;
        uint64_t offset;
        if (p_i < large_partition_count)
        {
            size = partition_size+1;
            offset = p_i*(partition_size+1);
        }else
        {
            size = partition_size;
            offset = large_partition_count*(partition_size+1) + (p_i-large_partition_count)*partition_size;
        }

        for (uint64_t i = offset; i < (offset+size); i++)
        {
            labels_out[ordering[i]] = p_i;
        }       
        
        
    }
    return;

}

void GetSamplesFromOrdered(std::vector<uint64_t> &order, std::vector<uint64_t> &input_arr, uint64_t sample_count, std::vector<uint64_t> &samples_out){
    assert(input_arr.size() == order.size());
    assert(samples_out.size() == sample_count);
    for (uint64_t sample_i = 0; sample_i < sample_count; sample_i++)
    {
        samples_out[sample_i] = input_arr[order[sample_i*(input_arr.size()/sample_count)]];
    }
    
}

