#include <util.hpp>
#include "mpi.h"
#include <vector>

PartitionStatus GetPtScotchPartitions(std::vector<uint64_t>& vtxdist, std::vector<uint64_t>& xadj,
                                      std::vector<uint64_t>& adjncy, uint64_t num_vertices_local,
                                      uint64_t num_vertices_global, int partition_count,
                                      std::vector<uint16_t>& partition_labels_out, MPI_Comm comm);