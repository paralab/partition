#ifndef BFS_COM_H
#define BFS_COM_H

#include <vector>
#include <mpi.h>

bool UpdateGhost(std::vector<unsigned long>& bfs_state, unsigned long graph_part_size, std::vector<int>& receive_counts, std::vector<std::vector<unsigned long>>& send_indices_per_proc, MPI_Comm& com);

#endif