#ifndef METIS_UTIL_H
#define METIS_UTIL_H

std::vector<uint64_t> GetMETISPartitions(std::vector<uint64_t> &xadj, std::vector<uint64_t> &adjncy, int32_t num_vertices,
                                         int32_t partition_count);

#endif