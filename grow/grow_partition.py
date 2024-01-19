import networkx as nx
import typing
import random
import math
import copy
from collections import defaultdict 

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



"""
assumes partitions are labelled 0,1,2,...., (partition_count-1)
"""
def get_partition_cut_size(G: nx.graph, vertex_to_partition: typing.Dict[int,int], partition_i: int, vertices_in_partition: typing.Set[int]) -> int:
    partition_cuts  = 0

    for u in vertices_in_partition:
        for neigh in G.neighbors(u):
            if neigh in vertex_to_partition and vertex_to_partition[neigh] != partition_i:
                partition_cuts+=1


    # for u,v in G.edges:
    #     if u not in vertex_to_partition or v not in vertex_to_partition:
    #         continue
    #     part_u = vertex_to_partition[u]
    #     part_v = vertex_to_partition[v]
    #     if part_u != part_v:
    #         if part_u == partition_i or part_v == partition_i:
    #             partition_cuts+=1
    
    return partition_cuts


"""
assumes an undirected, connected graph
"""
def get_grow_partitions_stretch_to_frontier(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:

    # grow from partition centers
    # ==================================

    # grow_random_seed_indices = [random.randint(0, len(elems)) for _ in range(partition_count)]
    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    optimal_partition_size = G.number_of_nodes()//partition_count

    vertex_to_partition = {}
    vertex_to_distance = {}
    partition_to_frontier = {}
    partition_to_vertices = {}
    partition_growing_status = [True for _ in range(partition_count)]
    # partition_to_current_edge_cut = [set() for _ in range(partition_count)]
    partition_to_current_edge_cut_size = [0 for _ in range(partition_count)]

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
        vertex_to_distance[s] = 0
        partition_to_frontier[p_i] = set([s])
        partition_to_vertices[p_i] = set([s])


    while any(partition_growing_status):
        # print(partition_to_vertices)
        pass
        for p_i in range(partition_count):
            if not partition_growing_status[p_i]:
                continue
            partition_size = len(partition_to_vertices[p_i])
            frontier_size = len(partition_to_frontier[p_i])
            if frontier_size == 0:
                partition_growing_status[p_i] = False   # other partitions already has acquired this frontier
                continue

            current_cut_size = get_partition_cut_size(G, vertex_to_partition, p_i, partition_to_vertices[p_i])
            distance_delta = current_cut_size + frontier_size


            new_frontier = set()
            for curr_frontier_vertex in partition_to_frontier[p_i]:
                    # if random.random() <= probability:      # advance the frontier with a probability
                    current_distance = vertex_to_distance[curr_frontier_vertex]
                    for neigh in G.neighbors(curr_frontier_vertex):
                        if neigh in partition_to_vertices[p_i]:     # already discovered by this partition
                            continue
                        if (neigh not in vertex_to_partition) or vertex_to_distance[neigh] > (current_distance + distance_delta):       # if not assigned to any partition
                            
                            if neigh in vertex_to_partition:        # removing from other partitions and other frontiers
                                other_partition = vertex_to_partition[neigh]
                                partition_to_vertices[other_partition].remove(neigh)
                                partition_to_frontier[other_partition].discard(neigh)
                            
                            vertex_to_distance[neigh] = current_distance + distance_delta
                            vertex_to_partition[neigh] = p_i       # then add to current growing partition
                            partition_to_vertices[p_i].add(neigh)
                            new_frontier.add(neigh)
                        # else:
                        #     if vertex_to_partition[neigh] != p_i:
                        #         partition_to_current_edge_cut[p_i].add(neigh)
                    # else:
                    #     new_frontier.add(curr_frontier_vertex)      # if probability failed, keep it for the next frontier
            partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0:
                partition_growing_status[p_i] = False


    assert len(vertex_to_partition) == G.number_of_nodes(), "error: some nodes have been left out of partitioning"

    return vertex_to_partition


"""
assumes partitions are labelled as 0,1,2,..... partition_count-1
"""
def all_partitions_reached_threshold(vertex_to_partition_distances: typing.Dict[int,typing.Dict[int,int]], partition_count: int, threshold: int) -> bool:
    partition_sizes = [0 for _ in range(partition_count)]
    for vertex in vertex_to_partition_distances:
        for reached_BFS_i in vertex_to_partition_distances[vertex]:
            partition_sizes[reached_BFS_i]+=1
    
    print(partition_sizes)

    for part_size in partition_sizes:
        if part_size < threshold:
            return False

    return True

"""
assumes partitions are labelled as 0,1,2,..... partition_count-1
"""
def adjust_distances_to_layers(vertex_to_partition_distances: typing.Dict[int,typing.Dict[int,int]], partition_count: int) -> typing.Dict[int,typing.Dict[int,int]]:
    BFS_to_layer_to_size = [{} for _ in range(partition_count)]
    for vertex in vertex_to_partition_distances:
        for BFS_i in vertex_to_partition_distances[vertex]:
            distance = vertex_to_partition_distances[vertex][BFS_i]
            if distance not in BFS_to_layer_to_size[BFS_i]:
                BFS_to_layer_to_size[BFS_i][distance] = 0
            BFS_to_layer_to_size[BFS_i][distance]+=1


    BFS_to_layer_to_cumulative_sizes = [None for _ in range(partition_count)]
    for BFS_i in range(partition_count):


        BFS_layer_to_cumulative_size = {}
        layers = sorted(BFS_to_layer_to_size[BFS_i].keys())
        prev_total = 0
        for layer in layers:
            BFS_layer_to_cumulative_size[layer] = prev_total + BFS_to_layer_to_size[BFS_i][layer]
            prev_total+= BFS_to_layer_to_size[BFS_i][layer]
        # print(BFS_i, BFS_layer_to_cumulative_size, "\n")
        BFS_to_layer_to_cumulative_sizes[BFS_i] = BFS_layer_to_cumulative_size
    
    vertex_to_partition_distance_size_estimate = {}

    for vertex in vertex_to_partition_distances:
        vertex_to_partition_distance_size_estimate[vertex] = {}
        for BFS_i in vertex_to_partition_distances[vertex]:
            vertex_to_partition_distance_size_estimate[vertex][BFS_i] = BFS_to_layer_to_cumulative_sizes[BFS_i][vertex_to_partition_distances[vertex][BFS_i]]
    # print(vertex_to_partition_distances[1])
    # print(vertex_to_partition_distance_size_estimate[1])
    return vertex_to_partition_distance_size_estimate

def get_local_grow_partitions_size_proxy(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:

    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    """
    each vertex is assigned a dict
    vertex_to_partition_distances = {
        vertex_0 = {
            BFS_0 :  distance,
            BFS_4 : distance
        }
    }
    """
    vertex_to_partition_distances = {}       # each vertex will be assigned a tuple (current partition, distance to current partition)
    for vertex in list(G.nodes):
        vertex_to_partition_distances[vertex] = {}
    for BFS_i,root in enumerate(seeds):
        vertex_to_partition_distances[root][BFS_i] = 0
    # is_stable = False
    grow_threshold = 3 * int(G.number_of_nodes()/partition_count)
    while not all_partitions_reached_threshold(vertex_to_partition_distances, partition_count, grow_threshold):
        vertex_to_partition_distances_new_copy = copy.deepcopy(vertex_to_partition_distances)
        for vertex in list(G.nodes):
            for neigh in G.neighbors(vertex):
                for BFS_i in vertex_to_partition_distances[neigh]:
                    # if BFS_i not in vertex_to_partition_distances[neigh]:       # this BFS_i has not reached this neighbor
                    #     continue
                    if (BFS_i not in vertex_to_partition_distances_new_copy[vertex]) or ((vertex_to_partition_distances[neigh][BFS_i] + 1) < vertex_to_partition_distances_new_copy[vertex][BFS_i]):
                        vertex_to_partition_distances_new_copy[vertex][BFS_i] = vertex_to_partition_distances[neigh][BFS_i] + 1

        vertex_to_partition_distances = vertex_to_partition_distances_new_copy

    # print("\n\n")
    # print(vertex_to_partition_distances)
    for vertex in list(G.nodes):
        assert len(vertex_to_partition_distances[vertex]) > 0, f"vertex {vertex} is not reached by any BFS"
    
    
    vertex_to_partition_distances_adjusted_to_sizes = adjust_distances_to_layers(vertex_to_partition_distances,partition_count)
    vertex_to_partition_label = {}

    for vertex in vertex_to_partition_distances_adjusted_to_sizes:
        vertex_to_partition_label[vertex] = min(vertex_to_partition_distances_adjusted_to_sizes[vertex], key=vertex_to_partition_distances_adjusted_to_sizes[vertex].get)

    return vertex_to_partition_label
    
def get_inside_edge_counts(G: nx.Graph,partition_count: int, current_partitioning: typing.Dict[int,int]) -> typing.List[int]:
    inside_edge_counts = [0 for _ in range(partition_count)]

    for u,v in G.edges:
        if not (u in current_partitioning and v in current_partitioning):
            continue
        part_u = current_partitioning[u]
        part_v = current_partitioning[v]
        if part_u == part_v:
            inside_edge_counts[part_u]+=1

    print(inside_edge_counts)

    return inside_edge_counts



def get_local_grow_partitions_size_proxy_with_fennel(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:

    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    """
    each vertex is assigned a dict
    vertex_to_partition_distances = {
        vertex_0 = {
            BFS_0 :  distance,
            BFS_4 : distance
        }
    }
    """
    vertex_to_partition_distances = {}       # each vertex will be assigned a tuple (current partition, distance to current partition)
    for vertex in list(G.nodes):
        vertex_to_partition_distances[vertex] = {}
    for BFS_i,root in enumerate(seeds):
        vertex_to_partition_distances[root][BFS_i] = 0
    # is_stable = False
    grow_threshold = 5 * int(G.number_of_nodes()/partition_count)
    while not all_partitions_reached_threshold(vertex_to_partition_distances, partition_count, grow_threshold):
        vertex_to_partition_distances_new_copy = copy.deepcopy(vertex_to_partition_distances)
        for vertex in list(G.nodes):
            for neigh in G.neighbors(vertex):
                for BFS_i in vertex_to_partition_distances[neigh]:
                    # if BFS_i not in vertex_to_partition_distances[neigh]:       # this BFS_i has not reached this neighbor
                    #     continue
                    if (BFS_i not in vertex_to_partition_distances_new_copy[vertex]) or ((vertex_to_partition_distances[neigh][BFS_i] + 1) < vertex_to_partition_distances_new_copy[vertex][BFS_i]):
                        vertex_to_partition_distances_new_copy[vertex][BFS_i] = vertex_to_partition_distances[neigh][BFS_i] + 1

        vertex_to_partition_distances = vertex_to_partition_distances_new_copy

    # print("\n\n")
    # print(vertex_to_partition_distances)
    for vertex in list(G.nodes):
        assert len(vertex_to_partition_distances[vertex]) > 0, f"vertex {vertex} is not reached by any BFS"

    ideal_partition_size = int(G.number_of_nodes()/partition_count)
    
    
    vertex_to_partition_distances_adjusted_to_sizes = adjust_distances_to_layers(vertex_to_partition_distances,partition_count)
    vertex_to_partition_label = {}

    partition_current_sizes = [0 for _ in range(partition_count)]


    for vertex in vertex_to_partition_distances_adjusted_to_sizes:
        best_partition = None
        best_partition_my_layer_cumulative_size  = float('inf')

        for BFS_i in vertex_to_partition_distances_adjusted_to_sizes[vertex]:
            layer_cumulative_size = vertex_to_partition_distances_adjusted_to_sizes[vertex][BFS_i]
            if layer_cumulative_size > ideal_partition_size:
                continue
            if layer_cumulative_size < best_partition_my_layer_cumulative_size:
                best_partition = BFS_i
        
        if best_partition!=None:
            vertex_to_partition_label[vertex] = best_partition
            partition_current_sizes[best_partition]+=1

    print(f"initially {len(vertex_to_partition_label)} vertices assigned out of {G.number_of_nodes()} \t({100*len(vertex_to_partition_label)/G.number_of_nodes()}%)")
    
    inside_edge_counts = get_inside_edge_counts(G,partition_count,vertex_to_partition_label)
    gamma = 2.0
    alpha = G.number_of_edges() * ((partition_count**(gamma - 1)) / (G.number_of_nodes()**(gamma)))

    
    while len(vertex_to_partition_label) < G.number_of_nodes():
        for vertex in list(G.nodes):
            if vertex in vertex_to_partition_label:
                continue
            # else:
            #     vertex_to_partition_label[vertex] = 0
            #     continue
            a_neighbor_is_assigned = False
            for neigh in G.neighbors(vertex):
                if neigh in vertex_to_partition_label:
                    a_neighbor_is_assigned = True
                    break
            if not a_neighbor_is_assigned:      # if there is no neighbor already assigned, we will process this vertex later
                continue

            best_partition=None
            best_partition_g_value=float("-inf")

            for BFS_i in range(partition_count):
                BFS_i_g_value = 0
                for BFS_j in range(partition_count):
                    inside_edges = inside_edge_counts[BFS_j]
                    part_size = partition_current_sizes[BFS_j]

                    if BFS_i == BFS_j:
                        part_size+=1
                        for neigh in G.neighbors(vertex):
                            if (neigh in vertex_to_partition_label) and (vertex_to_partition_label[neigh] ==BFS_i):
                                inside_edges+=1
                    BFS_i_g_value += (inside_edges - alpha*(part_size**gamma)) 
                if BFS_i_g_value > best_partition_g_value:
                    best_partition = BFS_i
                    best_partition_g_value = BFS_i_g_value
            
            vertex_to_partition_label[vertex] = best_partition

            partition_current_sizes[best_partition]+=1   
            for neigh in G.neighbors(vertex):
                if (neigh in vertex_to_partition_label) and (vertex_to_partition_label[neigh] ==best_partition):
                    inside_edge_counts[best_partition]+=1



    return vertex_to_partition_label