#ifndef METIS_UTIL_H
#define METIS_UTIL_H

#include "../util/util.hpp"

std::vector<uint64_t> GetMETISPartitions(std::vector<uint64_t> &xadj, std::vector<uint64_t> &adjncy, int32_t num_vertices,
                                         int32_t partition_count);
PartitionStatus GetParMETISPartitions(std::vector<uint64_t>& vtxdist, std::vector<uint64_t>& xadj, std::vector<uint64_t>& adjncy,
                          uint64_t num_vertices_local, int32_t partition_count,
                          std::vector<uint16_t>& partition_labels_out, MPI_Comm comm);
#endif