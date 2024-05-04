#ifndef DIST_GRAPH_H
#define DIST_GRAPH_H


#include <vector>
#include <cstdint>
#include "../mesh-util/mesh-util.hpp"
#include <stdint.h>


#define DIST_GRAPH_BFS_NO_LABEL UINT16_MAX
#define DIST_GRAPH_BFS_INFINITY UINT32_MAX


// typedef uint64_t distgraph_vertex_t;

struct BFSValue
{
    uint16_t label;
    uint32_t distance;
};

class DistGraph 
{
private:
    MPI_Comm comm;
    uint64_t local_count;
    uint64_t ghost_count;

    /**
     * Local numbering (i.e. local index) of nodes will be in range [0, size(own_elements)+ size(ghost_elements)]
     * Hence ghost elements will have largest local indices
    */
    std::vector<uint64_t> local_degrees;

    //CSR row pointers
    std::vector<uint64_t> local_xdj;
    

    //CSR column offsets
    std::vector<uint64_t> local_adjncy;

    std::vector<int> ghost_counts;
    std::vector<int> ghost_counts_scanned;

    bool RunLocalBFSToStable(std::vector<BFSValue>& bfs_vector);


public:
    DistGraph(const std::vector<ElementWithCoord>& own_elements,const std::vector<ElementWithTag>& ghost_elements,
                     const std::vector<std::pair<ElementWithTag, ElementWithTag>>& local_connectivity,
                     const std::vector<std::pair<ElementWithTag, ElementWithTag>>& boundary_connectivity,
                     const std::vector<uint64_t>& proc_element_counts,
                     const std::vector<uint64_t>& proc_element_counts_scanned, 
                     const std::vector<int>& ghost_element_counts,
                     MPI_Comm comm) ;
    std::string GraphToString();
    void Erase();
    void PartitionBFS(std::vector<uint16_t> partition_labels_out);
    // ~DistGraph();


};

#endif