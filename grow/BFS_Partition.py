import networkx as nx
import typing

"""
assumes an undirected, connected graph
"""
def get_BFS_partitions(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:
    # %%

    # BFS from partition centers
    # ==================================

    # BFS_random_seed_indices = [random.randint(0, len(elems)) for _ in range(partition_count)]
    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"

    vertex_to_BFS_partition = {}
    BFS_partition_to_frontier = {}
    BFS_partition_growing_status = [True for _ in range(partition_count)]

    #initial seeds
    for p_i, s in enumerate(seeds):
        vertex_to_BFS_partition[s] = p_i
        BFS_partition_to_frontier[p_i] = set([s])

    while any(BFS_partition_growing_status):
        for p_i in range(partition_count):
            if not BFS_partition_growing_status[p_i]:
                continue
            new_frontier = set()
            for curr_frontier_elem in BFS_partition_to_frontier[p_i]:
                    for neigh in G.neighbors(curr_frontier_elem):
                        if neigh not in vertex_to_BFS_partition:       # if not assigned to any partition
                            vertex_to_BFS_partition[neigh] = p_i       # then add to current growing partition
                            new_frontier.add(neigh)  
            BFS_partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0:
                BFS_partition_growing_status[p_i] = False


    assert len(vertex_to_BFS_partition) == G.number_of_nodes(), "error: some nodes have been left out of partitioning"

    print("BFS partitioning done")
    return vertex_to_BFS_partition