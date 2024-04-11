import networkx as nx
import typing
import random
import math
import copy
from collections import defaultdict
from datetime import datetime

from BFS_Partition import get_BFS_partitions 

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
    print(BFS_to_layer_to_cumulative_sizes[0])
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
            if layer_cumulative_size > (ideal_partition_size):
                continue
            if layer_cumulative_size < best_partition_my_layer_cumulative_size:
                best_partition = BFS_i
                best_partition_my_layer_cumulative_size = layer_cumulative_size
        
        if best_partition!=None:
            vertex_to_partition_label[vertex] = best_partition
            partition_current_sizes[best_partition]+=1

    print(f"initially {len(vertex_to_partition_label)} vertices assigned out of {G.number_of_nodes()} \t({100*len(vertex_to_partition_label)/G.number_of_nodes()} %)")
    
    inside_edge_counts = get_inside_edge_counts(G,partition_count,vertex_to_partition_label)
    gamma = 1.75
    alpha = G.number_of_edges() * ((partition_count**(gamma - 1)) / (G.number_of_nodes()**(gamma)))

    
    while len(vertex_to_partition_label) < G.number_of_nodes():
        for vertex in list(G.nodes):
            if vertex in vertex_to_partition_label:
                continue
            # else:
            #     vertex_to_partition_label[vertex] = -1
            #     continue
            # a_neighbor_is_assigned = False
            # neighbor_count = 0
            assigned_neighbor_count = 0
            for neigh in G.neighbors(vertex):
                # neighbor_count+=1
                if neigh in vertex_to_partition_label:
                    assigned_neighbor_count+=1
                    # a_neighbor_is_assigned = True
                    # break
            # if not a_neighbor_is_assigned:      # if there is no neighbor already assigned, we will process this vertex later
            #     continue
            if assigned_neighbor_count < 1:
                continue
            print(f"processing vertex {vertex}")

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
                    BFS_i_g_value += inside_edges
                    # BFS_i_g_value += (inside_edges - alpha*(part_size**gamma)) 
                if BFS_i_g_value > best_partition_g_value:
                    best_partition = BFS_i
                    best_partition_g_value = BFS_i_g_value
            
            vertex_to_partition_label[vertex] = best_partition

            partition_current_sizes[best_partition]+=1   
            for neigh in G.neighbors(vertex):
                if (neigh in vertex_to_partition_label) and (vertex_to_partition_label[neigh] ==best_partition):
                    inside_edge_counts[best_partition]+=1



    return vertex_to_partition_label




def get_local_BFS_rebalanced_partitions(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:
    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    vertex_to_partition_distance = {}       # each vertex will be assigned a tuple (current partition, distance to current partition)
    for BFS_i,root in enumerate(seeds):
        vertex_to_partition_distance[root] = (BFS_i, 0)

    

    is_stable = False
    while not is_stable:
        is_stable = True
        vertex_to_partition_distance_new_copy = copy.deepcopy(vertex_to_partition_distance)
        for vertex in list(G.nodes):
            for neigh in G.neighbors(vertex):
                if neigh not in vertex_to_partition_distance:
                    continue
                    
                if (vertex not in vertex_to_partition_distance_new_copy) or  ((vertex_to_partition_distance[neigh][1] + 1) < vertex_to_partition_distance_new_copy[vertex][1]):
                    vertex_to_partition_distance_new_copy[vertex] = (vertex_to_partition_distance[neigh][0], vertex_to_partition_distance[neigh][1]+1)
                    is_stable = False
        vertex_to_partition_distance = vertex_to_partition_distance_new_copy

    
    vertex_to_partition_label = {}

    for vertex in vertex_to_partition_distance:
        partititon = vertex_to_partition_distance[vertex][0]
        vertex_to_partition_label[vertex] = partititon


    ideal_partition_size = int(G.number_of_nodes()/partition_count)
    reassignments = 0
    reassignment_rounds = 5
    for r_i in range(reassignment_rounds):
        partition_sizes = [0 for _ in range(partition_count)]
        for vertex in vertex_to_partition_label:
            partititon = vertex_to_partition_label[vertex]
            partition_sizes[partititon]+=1

        vertex_to_partition_label_old_copy = copy.deepcopy(vertex_to_partition_label)

        for vertex in vertex_to_partition_label:
            if partition_sizes[vertex_to_partition_label[vertex]] <= (ideal_partition_size*0.9):
                continue

            best_partition_to_reassign = None
            best_partition_to_reassign_size = float('inf')
            for neigh in G.neighbors(vertex):

                partition = vertex_to_partition_label_old_copy[neigh]
                partition_size = partition_sizes[partititon]

                if vertex_to_partition_label[vertex] == partition:
                    continue        # neighbor in same partition, skip

                if partition_size > (ideal_partition_size*0.8):
                    continue        # do not reassign to already large partition

                if partition_size < best_partition_to_reassign_size:
                    best_partition_to_reassign = partition
                    best_partition_to_reassign_size = partition_size

            if best_partition_to_reassign!=None:
                vertex_to_partition_label[vertex] = best_partition_to_reassign
                reassignments+=1
      

    print(f"reassigned {reassignments} vertices")
    

    return vertex_to_partition_label


def get_grow_partitions_stop_when_ideal(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:
    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"

    vertex_to_BFS_partition = {}
    BFS_partition_to_frontier = {}
    BFS_partition_growing_status = [True for _ in range(partition_count)]

    #initial seeds
    for p_i, s in enumerate(seeds):
        vertex_to_BFS_partition[s] = p_i
        BFS_partition_to_frontier[p_i] = set([s])

    partition_sizes = [1 for _ in range(partition_count)]
    ideal_partition_size = int(G.number_of_nodes()/partition_count)

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
                            partition_sizes[p_i]+=1  
            BFS_partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0 or partition_sizes[p_i]>=(ideal_partition_size*0.9):
                BFS_partition_growing_status[p_i] = False
    print(f"initial pass assigned {len(vertex_to_BFS_partition)} vertices out of {G.number_of_nodes()}\t{100*len(vertex_to_BFS_partition)/G.number_of_nodes()} %")
    print(f"initial partition sizes: {partition_sizes}")
    print(f"initial partition size ratios: {[round(x/ideal_partition_size,2) for x in partition_sizes]}")
    # assigning any left out vertices
    while(len(vertex_to_BFS_partition) < G.number_of_nodes()):
        for vertex in list(G.nodes):
            if vertex in vertex_to_BFS_partition:
                continue        # already assigned
            vertex_to_BFS_partition[vertex] = 0
            for neigh in G.neighbors(vertex):
                if neigh not in vertex_to_BFS_partition:
                    continue
                vertex_to_BFS_partition[vertex] = vertex_to_BFS_partition[neigh]
                partition_sizes[vertex_to_BFS_partition[neigh]]+=1
                break
    print(f"final partition size ratios: {[round(x/ideal_partition_size,2) for x in partition_sizes]}")
    return vertex_to_BFS_partition


def get_grow_partitions_2_passes_rand_seeds(G: nx.Graph ,seeds_given: list[int],partition_count: int) -> typing.Dict[int,int]:
    all_vertices = list(G.nodes)
    random.seed(datetime.now().timestamp())
    seeds= [all_vertices[random.randint(0,G.number_of_nodes())] for _ in range(partition_count)]
    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    
    first_pass_vertex_to_partition_distance = {}       # each vertex will be assigned a tuple (current partition, distance to current partition)
    for BFS_i,root in enumerate(seeds):
        first_pass_vertex_to_partition_distance[root] = (BFS_i, 0)

    is_stable = False
    while not is_stable:
        is_stable = True
        first_pass_vertex_to_partition_distance_new_copy = copy.deepcopy(first_pass_vertex_to_partition_distance)
        for vertex in list(G.nodes):
            for neigh in G.neighbors(vertex):
                if neigh not in first_pass_vertex_to_partition_distance:
                    continue
                    
                if (vertex not in first_pass_vertex_to_partition_distance_new_copy) or  ((first_pass_vertex_to_partition_distance[neigh][1] + 1) < first_pass_vertex_to_partition_distance_new_copy[vertex][1]):
                    first_pass_vertex_to_partition_distance_new_copy[vertex] = (first_pass_vertex_to_partition_distance[neigh][0], first_pass_vertex_to_partition_distance[neigh][1]+1)
                    is_stable = False
        first_pass_vertex_to_partition_distance = first_pass_vertex_to_partition_distance_new_copy

    first_pass_partition_sizes = [0 for _ in range(partition_count)]
    first_pass_vertex_to_partition_label = {}
    for vertex in first_pass_vertex_to_partition_distance:
        first_pass_partition_sizes[first_pass_vertex_to_partition_distance[vertex][0]]+=1
        first_pass_vertex_to_partition_label[vertex] = first_pass_vertex_to_partition_distance[vertex][0]
    
    ideal_partition_size = int(G.number_of_nodes()/partition_count)
    second_pass_step_sizes = [(x/ideal_partition_size)**0.333 for x in first_pass_partition_sizes]
    
    # print(f"first pass sizes\t: {first_pass_partition_sizes}")
    # print(f"second_pass_step_sizes\t: {[round(x,2) for x in second_pass_step_sizes]}")


    # second pass, using new step sizes
    second_pass_vertex_to_partition_distance = {}       # each vertex will be assigned a tuple (current partition, distance to current partition)
    for BFS_i,root in enumerate(seeds):
        second_pass_vertex_to_partition_distance[root] = (BFS_i, 0)

    is_stable = False
    while not is_stable:
        is_stable = True
        second_pass_vertex_to_partition_distance_new_copy = copy.deepcopy(second_pass_vertex_to_partition_distance)
        for vertex in list(G.nodes):
            for neigh in G.neighbors(vertex):
                if neigh not in second_pass_vertex_to_partition_distance:
                    continue

                neighbor_partition = second_pass_vertex_to_partition_distance[neigh][0]
                neighbor_partition_distance = second_pass_vertex_to_partition_distance[neigh][1]
                neighbor_partition_step_size = second_pass_step_sizes[neighbor_partition]

                    
                if (vertex not in second_pass_vertex_to_partition_distance_new_copy) or  ((neighbor_partition_distance + neighbor_partition_step_size) < second_pass_vertex_to_partition_distance_new_copy[vertex][1]):
                    second_pass_vertex_to_partition_distance_new_copy[vertex] = (neighbor_partition, neighbor_partition_distance + neighbor_partition_step_size)
                    is_stable = False
        second_pass_vertex_to_partition_distance = second_pass_vertex_to_partition_distance_new_copy

    second_pass_vertex_to_partition_label = {}
    second_pass_partition_sizes = [0 for _ in range(partition_count)]
    for vertex in second_pass_vertex_to_partition_distance:
        second_pass_partition_sizes[second_pass_vertex_to_partition_distance[vertex][0]]+=1
        second_pass_vertex_to_partition_label[vertex] = second_pass_vertex_to_partition_distance[vertex][0]
    # print(f"second pass sizes\t: {second_pass_partition_sizes}")
    return second_pass_vertex_to_partition_label






def get_grow_partitions_2_passes_for_size_ratio(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:
    # all_vertices = list(G.nodes)
    # random.seed(datetime.now().timestamp())
    # seed_indices = random.sample(range(G.number_of_nodes()), partition_count)
    # seeds= [all_vertices[s_i] for s_i in seed_indices]
    
    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    
    first_pass_vertex_to_partition_distance = {}       # each vertex will be assigned a tuple (current partition, distance to current partition)
    for BFS_i,root in enumerate(seeds):
        first_pass_vertex_to_partition_distance[root] = (BFS_i, 0)

    is_stable = False
    while not is_stable:
        is_stable = True
        first_pass_vertex_to_partition_distance_new_copy = copy.deepcopy(first_pass_vertex_to_partition_distance)
        for vertex in list(G.nodes):
            for neigh in G.neighbors(vertex):
                if neigh not in first_pass_vertex_to_partition_distance:
                    continue
                    
                if (vertex not in first_pass_vertex_to_partition_distance_new_copy) or  ((first_pass_vertex_to_partition_distance[neigh][1] + 1) < first_pass_vertex_to_partition_distance_new_copy[vertex][1]):
                    first_pass_vertex_to_partition_distance_new_copy[vertex] = (first_pass_vertex_to_partition_distance[neigh][0], first_pass_vertex_to_partition_distance[neigh][1]+1)
                    is_stable = False
        first_pass_vertex_to_partition_distance = first_pass_vertex_to_partition_distance_new_copy

    first_pass_partition_sizes = [0 for _ in range(partition_count)]
    first_pass_vertex_to_partition_label = {}
    for vertex in first_pass_vertex_to_partition_distance:
        first_pass_partition_sizes[first_pass_vertex_to_partition_distance[vertex][0]]+=1
        first_pass_vertex_to_partition_label[vertex] = first_pass_vertex_to_partition_distance[vertex][0]
    
    # ideal_partition_size = int(G.number_of_nodes()/partition_count)
    


    first_pass_partition_cuts = [0 for _ in range(partition_count)]

    for u,v in G.edges:

        part_u = first_pass_vertex_to_partition_label[u]
        part_v = first_pass_vertex_to_partition_label[v]
        if part_u != part_v:
            first_pass_partition_cuts[part_u]+=1
            first_pass_partition_cuts[part_v]+=1

    first_pass_partition_size_to_cut_ratios = [s/c for s,c in zip(first_pass_partition_sizes,first_pass_partition_cuts)]

    # second_pass_step_sizes = [(max(first_pass_partition_size_to_cut_ratios)/x)**0.33 for x in first_pass_partition_size_to_cut_ratios]
    second_pass_step_sizes = [(x/max(first_pass_partition_cuts))**0.5 for x in first_pass_partition_cuts]
    
    print(f"first pass cuts\t: {first_pass_partition_cuts}")
    print(f"second_pass_step_sizes\t: {[round(x,2) for x in second_pass_step_sizes]}")


    # second pass, using new step sizes
    second_pass_vertex_to_partition_distance = {}       # each vertex will be assigned a tuple (current partition, distance to current partition)
    for BFS_i,root in enumerate(seeds):
        second_pass_vertex_to_partition_distance[root] = (BFS_i, 0)

    is_stable = False
    while not is_stable:
        is_stable = True
        second_pass_vertex_to_partition_distance_new_copy = copy.deepcopy(second_pass_vertex_to_partition_distance)
        for vertex in list(G.nodes):
            for neigh in G.neighbors(vertex):
                if neigh not in second_pass_vertex_to_partition_distance:
                    continue

                neighbor_partition = second_pass_vertex_to_partition_distance[neigh][0]
                neighbor_partition_distance = second_pass_vertex_to_partition_distance[neigh][1]
                neighbor_partition_step_size = second_pass_step_sizes[neighbor_partition]

                    
                if (vertex not in second_pass_vertex_to_partition_distance_new_copy) or  ((neighbor_partition_distance + neighbor_partition_step_size) < second_pass_vertex_to_partition_distance_new_copy[vertex][1]):
                    second_pass_vertex_to_partition_distance_new_copy[vertex] = (neighbor_partition, neighbor_partition_distance + neighbor_partition_step_size)
                    is_stable = False
        second_pass_vertex_to_partition_distance = second_pass_vertex_to_partition_distance_new_copy

    second_pass_vertex_to_partition_label = {}
    second_pass_partition_sizes = [0 for _ in range(partition_count)]
    for vertex in second_pass_vertex_to_partition_distance:
        second_pass_partition_sizes[second_pass_vertex_to_partition_distance[vertex][0]]+=1
        second_pass_vertex_to_partition_label[vertex] = second_pass_vertex_to_partition_distance[vertex][0]
    # print(f"second pass sizes\t: {second_pass_partition_sizes}")
    return second_pass_vertex_to_partition_label






"""
assumes an undirected, connected graph
"""
def get_grow_partitions_ordered_BFS(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:

    # all_vertices = list(G.nodes)
    # random.seed(datetime.now().timestamp())
    # seed_indices = random.sample(range(G.number_of_nodes()), partition_count)
    # seeds= [all_vertices[s_i] for s_i in seed_indices]

    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    optimal_partition_size = G.number_of_nodes()//partition_count

    vertex_to_partition = {}
    vertex_to_partition_order = {}
    vertex_to_partition_distance = {}
    partition_to_frontier = {}
    partition_to_vertices = {}
    partition_growing_status = [True for _ in range(partition_count)]
    partition_to_current_edge_cut = [set() for _ in range(partition_count)]



    #initial seeds
    for p_i, s in enumerate(seeds):
        vertex_to_partition[s] = p_i
        vertex_to_partition_order[s] = 1
        vertex_to_partition_distance[s] = 0
        partition_to_frontier[p_i] = set([s])
        partition_to_vertices[p_i] = set([s])


    while any(partition_growing_status):
        for p_i in range(partition_count):
            if not partition_growing_status[p_i]:
                continue
            partition_size = len(partition_to_vertices[p_i])


            new_frontier = set()
            ordering = float("-inf")
            for frontier_vertex in partition_to_frontier[p_i]:
                ordering = max(ordering, vertex_to_partition_order[frontier_vertex])
            ordering+=1
            for curr_frontier_vertex in partition_to_frontier[p_i]:
                    for neigh in G.neighbors(curr_frontier_vertex):
                        if neigh not in vertex_to_partition:       # if not assigned to any partition
                            vertex_to_partition[neigh] = p_i       # then add to current growing partition
                            vertex_to_partition_order[neigh] = ordering
                            ordering+=1
                            vertex_to_partition_distance[neigh] = vertex_to_partition_distance[curr_frontier_vertex] + 1
                            partition_to_vertices[p_i].add(neigh)
                            new_frontier.add(neigh)
                        elif vertex_to_partition[neigh] == p_i: # already in current growing partition
                            continue
                        else:
                            if (vertex_to_partition_order[neigh]/optimal_partition_size)*(vertex_to_partition_distance[neigh]) > (vertex_to_partition_order[curr_frontier_vertex]/optimal_partition_size)*(vertex_to_partition_distance[curr_frontier_vertex]) :
                                vertex_to_partition[neigh] = p_i
                                vertex_to_partition_order[neigh] = ordering
                                vertex_to_partition_distance[neigh] = vertex_to_partition_distance[curr_frontier_vertex] + 1
                                ordering+=1
                                partition_to_vertices[p_i].add(neigh)
                                new_frontier.add(neigh)

                                for other_p_i in range(partition_count):
                                    if other_p_i == p_i:
                                        continue
                                    partition_to_frontier[other_p_i].discard(neigh)
                            

            partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0:
                partition_growing_status[p_i] = False


    assert len(vertex_to_partition) == G.number_of_nodes(), "error: some nodes have been left out of partitioning"

    return vertex_to_partition



"""
assumes an undirected, connected graph
"""
def get_grow_partitions_noised_BFS(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:

    # all_vertices = list(G.nodes)
    random.seed(datetime.now().timestamp())
    # seed_indices = random.sample(range(G.number_of_nodes()), partition_count)
    # seeds= [all_vertices[s_i] for s_i in seed_indices]

    assert len(seeds) == partition_count, "len(seeds) should be equal to partition_count"
    optimal_partition_size = G.number_of_nodes()//partition_count

    vertex_to_partition = {}
    # vertex_to_partition_order = {}
    vertex_to_partition_distance = {}
    partition_to_frontier = {}
    partition_to_vertices = {}
    partition_growing_status = [True for _ in range(partition_count)]
    partition_to_current_edge_cut = [set() for _ in range(partition_count)]



    #initial seeds
    for p_i, s in enumerate(seeds):
        vertex_to_partition[s] = p_i
        # vertex_to_partition_order[s] = 1
        vertex_to_partition_distance[s] = 0
        partition_to_frontier[p_i] = set([s])
        partition_to_vertices[p_i] = set([s])


    while any(partition_growing_status):
        for p_i in range(partition_count):
            if not partition_growing_status[p_i]:
                continue
            partition_size = len(partition_to_vertices[p_i])


            new_frontier = set()
            # ordering = float("-inf")
            # for frontier_vertex in partition_to_frontier[p_i]:
            #     ordering = max(ordering, vertex_to_partition_order[frontier_vertex])
            # ordering+=1
            for curr_frontier_vertex in partition_to_frontier[p_i]:
                    increment = 1
                    for neigh in G.neighbors(curr_frontier_vertex):
                        if neigh not in vertex_to_partition:       # if not assigned to any partition
                            vertex_to_partition[neigh] = p_i       # then add to current growing partition
                            # vertex_to_partition_order[neigh] = ordering
                            # ordering+=1
                            vertex_to_partition_distance[neigh] = vertex_to_partition_distance[curr_frontier_vertex] + increment
                            partition_to_vertices[p_i].add(neigh)
                            new_frontier.add(neigh)
                        elif vertex_to_partition[neigh] == p_i: # already in current growing partition
                            continue
                        else:
                            # if (vertex_to_partition_order[neigh]/optimal_partition_size)*(vertex_to_partition_distance[neigh]) > (vertex_to_partition_order[curr_frontier_vertex]/optimal_partition_size)*(vertex_to_partition_distance[curr_frontier_vertex]) :
                            if (vertex_to_partition_distance[neigh]) > (vertex_to_partition_distance[curr_frontier_vertex]) :
                                vertex_to_partition[neigh] = p_i
                                # vertex_to_partition_order[neigh] = ordering
                                vertex_to_partition_distance[neigh] = vertex_to_partition_distance[curr_frontier_vertex] + increment/2
                                # ordering+=1
                                partition_to_vertices[p_i].add(neigh)
                                new_frontier.add(neigh)

                                for other_p_i in range(partition_count):
                                    if other_p_i == p_i:
                                        continue
                                    partition_to_frontier[other_p_i].discard(neigh)
                            

            partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0:
                partition_growing_status[p_i] = False


    assert len(vertex_to_partition) == G.number_of_nodes(), "error: some nodes have been left out of partitioning"

    return vertex_to_partition



def get_grow_partitions_2_passes_oversampled(G: nx.Graph ,seeds_oversampled: list[int],partition_count: int) -> typing.Dict[int,int]:
    
    oversampled_partition_count  = len(seeds_oversampled)
    print(f"oversampled_partition_count: {oversampled_partition_count}")

    vertex_to_BFS_partition = {}
    BFS_partition_to_frontier = {}
    BFS_partition_growing_status = [True for _ in range(oversampled_partition_count)]
    vertex_to_BFS_seed_distance = {}

    #initial seeds
    for p_i, s in enumerate(seeds_oversampled):
        vertex_to_BFS_partition[s] = p_i
        vertex_to_BFS_seed_distance[s] = 0
        BFS_partition_to_frontier[p_i] = set([s])

    while any(BFS_partition_growing_status):
        for p_i in range(oversampled_partition_count):
            if not BFS_partition_growing_status[p_i]:
                continue
            new_frontier = set()
            for curr_frontier_elem in BFS_partition_to_frontier[p_i]:
                    for neigh in G.neighbors(curr_frontier_elem):
                        if neigh not in vertex_to_BFS_partition:       # if not assigned to any partition
                            vertex_to_BFS_partition[neigh] = p_i       # then add to current growing partition
                            vertex_to_BFS_seed_distance[neigh] = vertex_to_BFS_seed_distance[curr_frontier_elem]+1
                            new_frontier.add(neigh)  
            BFS_partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0:
                BFS_partition_growing_status[p_i] = False

    # print(vertex_to_BFS_seed_distance)
    assert len(vertex_to_BFS_partition) == G.number_of_nodes(), "error: some nodes have been left out in initial partitioning"
    
    seeds_oversampled_pairwise_distances = {}
    for u,v in G.edges:
        if vertex_to_BFS_partition[u] != vertex_to_BFS_partition[v]:
            seed1 = seeds_oversampled[vertex_to_BFS_partition[u]]
            seed2 = seeds_oversampled[vertex_to_BFS_partition[v]]

            pair = tuple(sorted([seed1,seed2]))
            if pair not in seeds_oversampled_pairwise_distances:
                seeds_oversampled_pairwise_distances[pair] = float('inf')
            dist = vertex_to_BFS_seed_distance[u] + vertex_to_BFS_seed_distance[v]
            # print(dist)
            seeds_oversampled_pairwise_distances[pair] = min(seeds_oversampled_pairwise_distances[pair], dist)      # sometimes 2 partition may have different boundaries. we want to get the closest distance

    # print(f"seeds_oversampled_pairwise_distances:\n{seeds_oversampled_pairwise_distances}")

    seeds_oversampled_pairwise_distances_list = []

    for key in seeds_oversampled_pairwise_distances:
        seeds_oversampled_pairwise_distances_list.append([key, seeds_oversampled_pairwise_distances[key] ])

    seeds_oversampled_pairwise_distances_sorted = sorted(seeds_oversampled_pairwise_distances_list, key=lambda item : item[1])      # sort pairs using distance
    # print(f"seeds_oversampled_pairwise_distances_sorted:\n{seeds_oversampled_pairwise_distances_sorted}")

    remaining_seeds = [i for i in seeds_oversampled]


    # removing bad seeds
    while len(remaining_seeds) > partition_count and len(seeds_oversampled_pairwise_distances_sorted) > 0:
        # print(f"remaining_seeds: {remaining_seeds}")
        # print(f"seeds_oversampled_pairwise_distances_sorted:\n{seeds_oversampled_pairwise_distances_sorted}")

        closest_pair = seeds_oversampled_pairwise_distances_sorted[0][0]
        # print(f"considering pair: {closest_pair}")
        seed_to_remove = closest_pair[0]        # TODO: have a better way to choose one from the pair
        remaining_seeds.remove(seed_to_remove)
        seeds_oversampled_pairwise_distances_sorted = list(filter(lambda x: x[0][0]!= seed_to_remove and x[0][1]!= seed_to_remove, seeds_oversampled_pairwise_distances_sorted))
    if len(remaining_seeds) > partition_count:
        print("removing still remaning additional seeds, randomly")
        random.seed(datetime.now().timestamp())
        final_seeds = random.sample(remaining_seeds, partition_count)   # TODO: have a better to remove further remining seeds
    else:
        final_seeds = remaining_seeds
    return get_BFS_partitions(G ,final_seeds,partition_count)
