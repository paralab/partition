import networkx as nx
import typing
import random
import math

"""
assumes an undirected, connected graph
"""
def get_grow_partitions(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:
    # %%

    # grow from partition centers
    # ==================================

    # grow_random_seed_indices = [random.randint(0, len(elems)) for _ in range(partition_count)]
    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    optimal_partition_size = G.number_of_nodes()//partition_count

    vertex_to_partition = {}
    partition_to_frontier = {}
    partition_to_vertices = {}
    partition_growing_status = [True for _ in range(partition_count)]
    partition_to_current_edge_cut = [set() for _ in range(partition_count)]

    # grow_steps_n = 16
    # grow_step_size = (optimal_partition_size)*(2-0.5)/grow_steps_n
    # grow_thresholds = [optimal_partition_size*0.5 + grow_step_size*i for i in range(grow_steps_n)]
    # probabilities_upto_threshold = [1 - i*(1/grow_steps_n) for i in range(grow_steps_n)]


    # optimal_cut_size = (2*G.number_of_edges()/(math.log2(G.number_of_nodes()))**2)/partition_count
    # grow_steps_n = 16
    # grow_step_size = (optimal_cut_size)*(8-0.125)/grow_steps_n
    # grow_thresholds = [optimal_cut_size*0.125 + grow_step_size*i for i in range(grow_steps_n)]
    # probabilities_upto_threshold = [1 - i*(1/grow_steps_n) for i in range(grow_steps_n)]

    # print(grow_thresholds)
    # print(probabilities_upto_threshold)



    # random.seed(42)

    #initial seeds
    for p_i, s in enumerate(seeds):
        vertex_to_partition[s] = p_i
        partition_to_frontier[p_i] = set([s])
        partition_to_vertices[p_i] = set([s])


    while any(partition_growing_status):
        for p_i in range(partition_count):
            if not partition_growing_status[p_i]:
                continue
            partition_size = len(partition_to_vertices[p_i])
            # frontier_size = len(partition_to_frontier[p_i])
            # if partition_size <= optimal_partition_size/2:
            #     probability = 1
            # else:
            #     # probability = (optimal_partition_size/2)/partition_size
            # #     probability = 1/(partition_size**0.3 + len(partition_to_current_edge_cut[p_i])**0.5)
            #     probability = 1/(len(partition_to_current_edge_cut[p_i])+1)



            # probability = probabilities_upto_threshold[-1]  #default lowest probability
            # for thres_i in range(grow_steps_n):
            #     if (frontier_size + partition_to_current_edge_cut[p_i]) <= grow_thresholds[thres_i]:
            #         probability = probabilities_upto_threshold[thres_i]
            #         break


            # print(p_i, probability)

            # if not(random.random() <= probability):
            #     continue

            new_frontier = set()
            for curr_frontier_vertex in partition_to_frontier[p_i]:
                    # if random.random() <= probability:      # advance the frontier with a probability
                    for neigh in G.neighbors(curr_frontier_vertex):
                        if neigh not in vertex_to_partition:       # if not assigned to any partition
                            vertex_to_partition[neigh] = p_i       # then add to current growing partition
                            partition_to_vertices[p_i].add(neigh)
                            new_frontier.add(neigh)
                        else:
                            if vertex_to_partition[neigh] != p_i:
                                partition_to_current_edge_cut[p_i].add(neigh)
                    # else:
                    #     new_frontier.add(curr_frontier_vertex)      # if probability failed, keep it for the next frontier
            partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0:
                partition_growing_status[p_i] = False


    assert len(vertex_to_partition) == G.number_of_nodes(), "error: some nodes have been left out of partitioning"

    return vertex_to_partition