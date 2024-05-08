#include "metis.h"
#include "parmetis.h"
#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>

#include "metis-util.hpp"
#include "../util/util.hpp"

std::vector<uint64_t> GetMETISPartitions(std::vector<uint64_t> &xadj, std::vector<uint64_t> &adjncy, int32_t num_vertices,
                                         int32_t partition_count)
{
    std::vector<int32_t> xadj__(xadj.begin(), xadj.end());
    std::vector<int32_t> adjncy__(adjncy.begin(), adjncy.end());

    int options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    int32_t edgecut;
    std::vector<int32_t> partition_labels(num_vertices);
    int32_t ncon=1;
    auto start = std::chrono::high_resolution_clock::now();
    auto result = METIS_PartGraphKway(&num_vertices, &ncon, xadj__.data(), adjncy__.data(), nullptr, nullptr, nullptr,
                        &partition_count, nullptr, nullptr, options, &edgecut, partition_labels.data());
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    print_log("METIS time:\t", duration.count(), " ms");
    assert(result == METIS_OK);
    std::vector<uint64_t> partition_labels__(partition_labels.begin(), partition_labels.end());

    return partition_labels__;
}

int GetParMETISPartitions(std::vector<uint64_t>& vtxdist, std::vector<uint64_t>& xadj, std::vector<uint64_t>& adjncy,
                          uint64_t num_vertices_local, int32_t partition_count,
                          std::vector<uint16_t>& partition_labels_out, MPI_Comm comm) {

    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    std::vector<idx_t> vtxdist__(vtxdist.begin(), vtxdist.end());
    std::vector<idx_t> xadj__(xadj.begin(), xadj.end());
    std::vector<idx_t> adjncy__(adjncy.begin(), adjncy.end());
    int32_t ncon = 1;

    std::vector<idx_t> options = {0, 0, 0};

    idx_t edgecut = 0;
    std::vector<idx_t> partitions_labels(num_vertices_local);

    idx_t zero = 0;

    std::vector<real_t> tpwgts(partition_count,1/(static_cast<real_t>(partition_count)));

    std::vector<real_t> ubvec(1,1.05);

    MPI_Barrier(comm);
    auto start = std::chrono::high_resolution_clock::now();

    int return_code = ParMETIS_V3_PartKway(vtxdist__.data(), xadj__.data(), adjncy__.data(), nullptr, nullptr, &zero, &zero,
                                           &ncon, &partition_count, tpwgts.data(), ubvec.data(), options.data(), &edgecut,
                                           partitions_labels.data(), &comm);
    MPI_Barrier(comm);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    if (!my_rank)
    {
        print_log("ParMetis time:\t", duration.count(), " ms");
    }
    partition_labels_out.assign(partitions_labels.begin(), partitions_labels.end());
    return return_code;
}