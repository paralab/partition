
#include "util.hpp"
#include "scotch-util.hpp"
#include "ptscotch.h"
#include "scotch.h"
#include "mpi.h"
#include <chrono>
#include <cassert>
#include <stdexcept>

PartitionStatus GetPtScotchPartitions(std::vector<uint64_t>& vtxdist, std::vector<uint64_t>& xadj,
                                      std::vector<uint64_t>& adjncy, uint64_t num_vertices_local,
                                      uint64_t num_vertices_global, int partition_count,
                                      std::vector<uint16_t>& partition_labels_out, MPI_Comm comm) {
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    SCOTCH_Strat stradat;
    SCOTCH_Dgraph graph;

    std::vector<SCOTCH_Num> vertloctab(xadj.begin(), xadj.end());
    std::vector<SCOTCH_Num> edgeloctab(adjncy.begin(), adjncy.end());

    SCOTCH_Num local_total_arcs = static_cast<SCOTCH_Num>(edgeloctab.size()); // including arcs to/from ghosts

    SCOTCH_stratInit(&stradat); /* Default strategy will be used */
    SCOTCH_dgraphInit(&graph, comm);

    SCOTCH_dgraphBuild(&graph, 0, static_cast<SCOTCH_Num>(num_vertices_local),
                       static_cast<SCOTCH_Num>(num_vertices_local), vertloctab.data(), NULL, NULL, NULL,
                       local_total_arcs, local_total_arcs, edgeloctab.data(), NULL, NULL);
    int graph_status = SCOTCH_dgraphCheck(&graph);
    if (! (graph_status==0))
    {
        throw std::runtime_error("SOTCH error, graph is not consistent");
    }
    
    std::vector<SCOTCH_Num> partition_labels(num_vertices_local, 0);

    //warmup?
    MPI_Barrier(comm);
    SCOTCH_dgraphPart(&graph, static_cast<SCOTCH_Num>(partition_count), &stradat, partition_labels.data());

    MPI_Barrier(comm);
    auto start = std::chrono::high_resolution_clock::now();
    int return_code =
        SCOTCH_dgraphPart(&graph, static_cast<SCOTCH_Num>(partition_count), &stradat, partition_labels.data());
    MPI_Barrier(comm);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if (!my_rank) {
        print_log("PtScotch time:\t\t\t", duration.count(), " us");
    }
    SCOTCH_dgraphExit(&graph);
    SCOTCH_stratExit(&stradat);
    partition_labels_out.assign(partition_labels.begin(), partition_labels.end());
    return {.return_code = return_code, .time_us = static_cast<int>(duration.count())};
}