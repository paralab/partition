#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <sstream>
#include <iostream>

typedef unsigned long vertex_t;

class Graph
{
private:
    std::vector<std::unordered_set<vertex_t>> adj_list;
    std::unordered_map<vertex_t, unsigned long> vertex_to_index;
    std::vector<vertex_t> index_to_vetex;
    unsigned long graph_size;
    std::vector<unsigned long> bfs_status;
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
    void RunBFSToStable();
    std::vector<unsigned long> & GetBFSState();
    unsigned long GETBFSValue(vertex_t v);
    unsigned long GetSize();
    void Print();
    void PrintSingleBFS();

};


#endif

