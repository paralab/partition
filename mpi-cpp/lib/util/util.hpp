#ifndef UTIL_H
#define UTIL_H

#include "graph.hpp"

#include <iostream>
#include <sstream>
#include <utility>


bool IsValidVertex(int x, int graph_size);
bool IsMyVertex(vertex_t v, int my_rank, int num_tasks, int graph_size);
int GetMyBoundarySize(Graph& graph, int my_rank, int num_tasks, int graph_size, int N);
std::vector<vertex_t> GetMyBoundaryVertices(Graph& graph, int my_rank, int num_tasks, int graph_size, int N);
std::vector<std::vector<int>> GetMyGhostVerticesRelativeIndices(Graph& graph, std::vector<vertex_t> external_vertices, std::vector<int> external_vertex_counts, int my_rank, int num_tasks, int graph_size, int N);
void AddEdgesToGhostVertex(Graph& graph, vertex_t v, int my_rank, int num_tasks, int graph_size, int N);
void AssignPartitionLabelsInOrder(std::vector<uint64_t> &ordering, uint64_t count, uint64_t partition_count, std::vector<uint64_t> &labels_out);
void GetSamplesFromOrdered(std::vector<uint64_t> &order, std::vector<uint64_t> &input_arr, uint64_t sample_count, std::vector<uint64_t> &samples_out);
template <typename T> std::string VectorToString(std::vector<T> vec){
    std::ostringstream output;
    for (auto element : vec)
    {
        output << element << " ";
    }
    output << "\n";
    return output.str();
}



// template<typename ...Args>
// void print_log(Args && ...args)
// {
//     (std::cout << ... << args);
//     std::cout << std::endl;
// }

// Variadic template function to mimic std::cout with spaces between arguments
template<typename T, typename... Args>
void print_log(const T& first, const Args&... args) {
    std::ostringstream oss;
    oss << first;  // Output the first argument directly
    // Use fold expression to concatenate all arguments with spaces in between
    ((oss << ' ' << args), ...);
    // Output the concatenated string
    std::cout << oss.str() <<std::endl;
}

#endif
