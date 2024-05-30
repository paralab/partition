#ifndef DIST_GRAPH_H
#define DIST_GRAPH_H


#include <vector>
#include <cstdint>
#include "../mesh-util/mesh-util.hpp"
#include "../usort/dtypes.h"
#include "../util/util.hpp"

#include "../config.h"

#include <stdint.h>
#include <cstddef>

#if GRAPH_INDEXING_TYPE == 64
using graph_indexing_t = uint64_t;
#elif GRAPH_INDEXING_TYPE == 32
using graph_indexing_t = uint32_t;
#elif GRAPH_INDEXING_TYPE == 16
using graph_indexing_t = uint16_t;
#else
#error "Invalid GRAPH_INDEXING_TYPE specified. Allowed values: 16, 32, 64"
#endif






#if BFS_DISTANCE_TYPE == 64
using bfs_distance_t = uint64_t;
#define DIST_GRAPH_BFS_INFINITY UINT64_MAX

#elif BFS_DISTANCE_TYPE == 32
using bfs_distance_t = uint32_t;
#define DIST_GRAPH_BFS_INFINITY UINT32_MAX

#elif BFS_DISTANCE_TYPE == 16
using bfs_distance_t = uint16_t;
#define DIST_GRAPH_BFS_INFINITY UINT16_MAX

#else
#error "Invalid BFS_DISTANCE_TYPE specified. Allowed values: 16, 32, 64"
#endif


#if BFS_LABEL_TYPE == 64
using bfs_label_t = uint64_t;
#define DIST_GRAPH_BFS_NO_LABEL UINT64_MAX

#elif BFS_LABEL_TYPE == 32
using bfs_label_t = uint32_t;
#define DIST_GRAPH_BFS_NO_LABEL UINT32_MAX


#elif BFS_LABEL_TYPE == 16
using bfs_label_t = uint16_t;
#define DIST_GRAPH_BFS_NO_LABEL UINT16_MAX

#else
#error "Invalid BFS_LABEL_TYPE specified. Allowed values: 16, 32, 64"
#endif

struct BFSValue
{
    bfs_label_t label;
    bfs_distance_t distance;
};


template <>
class par::Mpi_datatype<BFSValue> {

    /** 
         @return the MPI_Datatype for the C++ datatype "BFSValue"
        **/
    public:
    static MPI_Datatype value() {
        static bool         first = true;
        static MPI_Datatype custom_mpi_type;

        if (first)
        {
            first = false;
            int block_lengths[2] = {1, 1};
            MPI_Datatype types[2] = {Mpi_datatype<bfs_label_t>::value(), Mpi_datatype<bfs_distance_t>::value()};
            MPI_Aint offsets[2];
            offsets[0] = offsetof(BFSValue, label);
            offsets[1] = offsetof(BFSValue, distance);


            MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
            MPI_Type_commit(&custom_mpi_type);
        }       

        return custom_mpi_type;
    }
};

std::ostream& operator<<(std::ostream& os, const BFSValue& obj);


/**
 * To be used for BFS value exchange for updated only ghost vertices
*/
struct GhostBFSValue
{
    bfs_label_t label;
    bfs_distance_t distance;
    graph_indexing_t offset;            // offset w.r.t original send buffer
};


template <>
class par::Mpi_datatype<GhostBFSValue> {

    /** 
         @return the MPI_Datatype for the C++ datatype "GhostBFSValue"
        **/
    public:
    static MPI_Datatype value() {
        static bool         first = true;
        static MPI_Datatype custom_mpi_type;

        if (first)
        {
            first = false;
            int block_lengths[3] = {1, 1, 1};
            MPI_Datatype types[3] = {Mpi_datatype<bfs_label_t>::value(), Mpi_datatype<bfs_distance_t>::value(), Mpi_datatype<graph_indexing_t>::value()};
            MPI_Aint offsets[3];
            offsets[0] = offsetof(GhostBFSValue, label);
            offsets[1] = offsetof(GhostBFSValue, distance);
            offsets[2] = offsetof(GhostBFSValue, offset);



            MPI_Type_create_struct(3, block_lengths, offsets, types, &custom_mpi_type);
            MPI_Type_commit(&custom_mpi_type);
        }       

        return custom_mpi_type;
    }
};

std::ostream& operator<<(std::ostream& os, const GhostBFSValue& obj);


#define DIST_GRAPH_PAGERANK_NO_LABEL UINT16_MAX

struct PageRankValue
{
    uint16_t label;
    float value;
};

template <>
class par::Mpi_datatype<PageRankValue> {

    /** 
         @return the MPI_Datatype for the C++ datatype "PageRankValue"
        **/
    public:
    static MPI_Datatype value() {
        static bool         first = true;
        static MPI_Datatype custom_mpi_type;

        if (first)
        {
            first = false;
            int block_lengths[2] = {1, 1};
            MPI_Datatype types[2] = {MPI_UINT16_T, MPI_FLOAT};
            MPI_Aint offsets[2];
            offsets[0] = offsetof(PageRankValue, label);
            offsets[1] = offsetof(PageRankValue, value);


            MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
            MPI_Type_commit(&custom_mpi_type);
        }       

        return custom_mpi_type;
    }
};

std::ostream& operator<<(std::ostream& os, const PageRankValue& obj);

class DistGraph 
{
private:
    MPI_Comm comm;
    uint64_t local_count;
    uint64_t ghost_count;
    uint64_t send_count;


    /**
     * Local numbering (i.e. local index) of nodes will be in range [0, size(own_elements)+ size(ghost_elements)]
     * Hence ghost elements will have largest local indices
    */
    std::vector<graph_indexing_t> local_degrees;

    //CSR row pointers
    std::vector<graph_indexing_t> local_xdj;
    

    //CSR column offsets
    std::vector<graph_indexing_t> local_adjncy;

    std::vector<uint64_t> dist_adjncy;

    std::vector<uint64_t> vtx_dist;


    std::vector<int> ghost_counts;
    std::vector<int> ghost_counts_scanned;

    std::vector<int> ghost_procs;

    std::vector<uint64_t> sending_scatter_map;
    std::vector<int> send_counts;
    std::vector<int> send_counts_scanned;
    std::vector<int> send_procs;


    void RunFirstBFSIteration(std::vector<BFSValue>& bfs_vector, graph_indexing_t seed, bfs_label_t label);
    void CalculateDistanceFromUpdatedGhosts(std::vector<bool>& ghost_updated, std::vector<graph_indexing_t>& frontier_buffer,
        std::vector<bfs_distance_t>& distances_out);
    bool RunLocalMultiBFSToStable(std::vector<BFSValue>& bfs_vector, std::vector<BFSValue>& bfs_vector_tmp, 
        std::vector<bool> vector_diff, bfs_distance_t ghost_min_update, std::vector<bfs_distance_t>& distance_from_updated_ghosts);
    void ExchangeUpdatedOnlyBFSGhost(std::vector<BFSValue>& ghost_send_buffer,
                                     std::vector<BFSValue>& ghost_send_buffer_prev,
                                     std::vector<BFSValue>& ghost_recv_buffer);
    void ExchangeUpdatedOnlyBFSCounts(std::vector<int>& updated_only_send_counts, std::vector<int>& updated_only_recv_counts_out);
    bool RunLocalMultiPageRankToStable(std::vector<PageRankValue>& pagerank_vector,
                std::vector<graph_indexing_t> vertex_degrees, const float min_relative_change);

    void GetVertexDegrees(std::vector<graph_indexing_t>& degrees_out);


public:
    DistGraph(const std::vector<ElementWithCoord>& own_elements,const std::vector<ElementWithTag>& ghost_elements,
                     const std::vector<std::pair<ElementWithTag, ElementWithTag>>& local_connectivity,
                     const std::vector<std::pair<ElementWithTag, ElementWithTag>>& boundary_connectivity,
                     const std::vector<uint64_t>& proc_element_counts,
                     const std::vector<uint64_t>& proc_element_counts_scanned, 
                     const std::vector<int>& ghost_element_counts,
                     MPI_Comm comm) ;
    std::string PrintLocal();
    std::string PrintDist();

    void Erase();
    PartitionStatus PartitionBFS(std::vector<uint16_t>& partition_labels_out);
    void PartitionPageRank(std::vector<uint16_t>& partition_labels_out);
    PartitionStatus PartitionParmetis(std::vector<uint16_t>& partition_labels_out);
    void GetPartitionMetrics(std::vector<uint16_t>& local_partition_labels, std::vector<uint32_t>& partition_sizes_out,
                             std::vector<uint32_t>& partition_boundaries_out);
    // ~DistGraph();


};

#endif