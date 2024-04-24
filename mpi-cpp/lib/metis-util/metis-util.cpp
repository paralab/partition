#include <metis.h>
#include <iostream>
#include <vector>
#include <cassert>
#include <chrono>

#include "metis-util.hpp"
#include "util.hpp"

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
