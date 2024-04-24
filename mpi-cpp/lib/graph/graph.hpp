#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>
#include <iostream>

typedef uint64_t vertex_t;

class Graph
{
private:
    std::vector<std::unordered_set<vertex_t>> adj_list;
    std::unordered_map<vertex_t, unsigned long> vertex_to_index;
    std::vector<vertex_t> index_to_vetex;
    unsigned long graph_size;
    std::vector<unsigned long> bfs_status;
    std::vector<uint64_t> multi_bfs_distances;
    std::vector<uint64_t> multi_bfs_labels;
    std::vector<double> multi_pagerank_weights;
    std::vector<uint64_t> multi_pagerank_labels;



    unsigned long infinity;
public:
    Graph(/* args */);
    ~Graph();
    void AddVertex(vertex_t v);
    bool AddEdge(vertex_t v, vertex_t u);
    std::vector<vertex_t>& GetVertices();
    std::unordered_set<vertex_t> GetNeighbors(vertex_t v);
    bool FindVertex(vertex_t v);
    void InitSingleBFS(vertex_t BFS_seed);
    void InitMultiBFS(std::vector<vertex_t>& seeds, uint64_t count);
    void RunBFSToStable();
    void RunMultiBFSToStable();
    std::vector<unsigned long> & GetBFSState();
    std::vector<uint64_t>& GetMultiBFSLabels();
    unsigned long GETBFSValue(vertex_t v);
    unsigned long GetSize();
    void Print();
    void PrintSingleBFS();
    std::vector<uint64_t> GetCSR_xadj();
    std::vector<uint64_t> GetCSR_adjncy();


};


#endif

