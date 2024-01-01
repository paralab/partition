#ifndef UTIL_H
#define UTIL_H

#include "graph.hpp"



bool IsValidVertex(int x, int graph_size);
bool IsMyVertex(vertex_t v, int my_rank, int num_tasks, int graph_size);
int GetMyBoundarySize(Graph& graph, int my_rank, int num_tasks, int graph_size, int N);
std::vector<vertex_t> GetMyBoundaryVertices(Graph& graph, int my_rank, int num_tasks, int graph_size, int N);
std::vector<std::vector<int>> GetMyGhostVerticesRelativeIndices(Graph& graph, std::vector<vertex_t> external_vertices, std::vector<int> external_vertex_counts, int my_rank, int num_tasks, int graph_size, int N);
void AddEdgesToGhostVertex(Graph& graph, vertex_t v, int my_rank, int num_tasks, int graph_size, int N);

template <typename T> std::string VectorToString(std::vector<T> vec){
    std::ostringstream output;
    for (auto element : vec)
    {
        output << element << " ";
    }
    output << "\n";
    return output.str();
}
#endif
