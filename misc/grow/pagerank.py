import networkx as nx
import typing
import copy
import random
from datetime import datetime 

"""
assumes an undirected, connected graph
"""
def get_pagerank_partitions(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:
    additional_rounds_limit = 5
    initial_weight_for_seeds = 100000000
    
    vertex_to_partition = {}       
    vertex_to_weight = {}

    for vertex in list(G.nodes):
        vertex_to_weight[vertex] = 0
    
    for BFS_i,root in enumerate(seeds):
        vertex_to_partition[root] = BFS_i
        vertex_to_weight[root] = initial_weight_for_seeds

    is_stable = False
    additional_rounds_count = 0
    while (not is_stable) and additional_rounds_count < additional_rounds_limit:
        print(f"reached\t: {100*len(vertex_to_partition)/G.number_of_nodes()} %")
        if len(vertex_to_partition) == G.number_of_nodes():
            additional_rounds_count+=1
        is_stable = True
        vertex_to_partition_new_copy = copy.deepcopy(vertex_to_partition)
        vertex_to_weight_new_copy = copy.deepcopy(vertex_to_weight)

        for vertex in list(G.nodes):
            partition_to_weight_sum = {}
            for neigh in G.neighbors(vertex):
                if neigh not in vertex_to_partition:
                    continue
                if vertex_to_partition[neigh] not in partition_to_weight_sum:
                    partition_to_weight_sum[vertex_to_partition[neigh]] = vertex_to_weight[neigh]/len(list(G.neighbors(neigh)))
                else:
                    partition_to_weight_sum[vertex_to_partition[neigh]] += vertex_to_weight[neigh]/len(list(G.neighbors(neigh)))
            if len(partition_to_weight_sum) == 0:
                continue
            if vertex_to_weight_new_copy[vertex]==0 or  max(partition_to_weight_sum.values())/vertex_to_weight_new_copy[vertex] > 1.2 :
                # print(f"relative change: {'new' if vertex_to_weight_new_copy[vertex]==0 else max(partition_to_weight_sum.values())/vertex_to_weight_new_copy[vertex] }")
                is_stable = False
                vertex_to_weight_new_copy[vertex] = max(partition_to_weight_sum.values())
                vertex_to_partition_new_copy[vertex]= max(partition_to_weight_sum, key=lambda k: partition_to_weight_sum.get(k))

        vertex_to_partition = vertex_to_partition_new_copy
        vertex_to_weight = vertex_to_weight_new_copy


    return vertex_to_partition
    

def run_pagerank(G: nx.Graph ,initial_weight: int = 1000, rounds: int = 15, damping_factor: float = 1) -> typing.Dict[int,float]:
    vertex_to_weight = {}

    for vertex in list(G.nodes):
        vertex_to_weight[vertex] = initial_weight
        
    for round in range(rounds):
        print(f"round: {round+1} out of {rounds}")
        vertex_to_weight_new_copy = copy.deepcopy(vertex_to_weight)
        for vertex in list(G.nodes):
            weights_from_neighbors = 0
            for neigh in G.neighbors(vertex):
                weights_from_neighbors+= vertex_to_weight[neigh]/len(list(G.neighbors(neigh)))          

            new_weight = damping_factor*weights_from_neighbors + (1-damping_factor)*(1/G.number_of_nodes())
            vertex_to_weight_new_copy[vertex] = new_weight

        vertex_to_weight = vertex_to_weight_new_copy
    return vertex_to_weight


