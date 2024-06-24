import networkx as nx
import typing
import random
import math
import itertools
from collections import defaultdict 
import statistics
"""
assumes an undirected, connected graph
"""
def get_parititions_by_oversampling_all_graph_BFS(G: nx.Graph , partition_count: int) -> typing.Dict[int,int]:
    oversample_factor = int(math.log2(partition_count))
    sample_count = oversample_factor + partition_count
    all_vertices = list(G.nodes)
    BFS_roots= [all_vertices[random.randint(0,G.number_of_nodes())] for _ in range(sample_count)]
    vertex_to_distances = {}        # each vertex has an array specifying distance to each BFS root

    for vertex in list(G.nodes):
        vertex_to_distances[vertex] = [None for _ in range(sample_count)]

    finished_count = 0

    BFS_root_to_index = {}          # to store the ordering

    for i,root in enumerate(BFS_roots):
        vertex_to_distances[root][i] = 0
        BFS_root_to_index[root] = i

    while finished_count < G.number_of_nodes():
        for vertex in list(G.nodes):
            if None not in vertex_to_distances[vertex]:     # every BFS has reached this node
                continue
            # print(vertex_to_distances[vertex])
            for BFS_i in range(sample_count):
                if vertex_to_distances[vertex][BFS_i] != None:      # this particular BFS_i has already reached this node
                    continue                                
                distances = []
                for neigh in G.neighbors(vertex):
                    if vertex_to_distances[neigh][BFS_i] != None:
                        distances.append(vertex_to_distances[neigh][BFS_i])
                if len(distances)>0:
                    vertex_to_distances[vertex][BFS_i] = min(distances) + 1

            if None not in vertex_to_distances[vertex]:     # every BFS has reached this node
                finished_count+=1

    for vertex in list(G.nodes):
        assert None not in vertex_to_distances[vertex], f"some BFS did not reach vertex {vertex}"

    BFS_root_pairs = itertools.combinations(BFS_roots,2)
    BFS_root_pairs_with_distances = []

    for pair in BFS_root_pairs:
        distance = vertex_to_distances[pair[0]][BFS_root_to_index[pair[1]]]    # we can choose other way around as well. both should be equal
        BFS_root_pairs_with_distances.append((pair,distance))

    BFS_root_pairs_with_distances_sorted = sorted(BFS_root_pairs_with_distances,key=lambda x: x[1])     # sort pairs by distance, ascending
    # print(BFS_root_pairs_with_distances_sorted)
    filtered_BFS_roots = BFS_roots[:]

    while len(filtered_BFS_roots) > partition_count:
        BFS_root_to_remove = BFS_root_pairs_with_distances_sorted[0][0][0]  # or 0 0 1
        filtered_BFS_roots.remove(BFS_root_to_remove)
        BFS_root_pairs_with_distances_sorted = list(filter(lambda x: x[0][0]!=BFS_root_to_remove and x[0][1]!=BFS_root_to_remove, BFS_root_pairs_with_distances_sorted))


    print(BFS_root_pairs_with_distances_sorted)
    print(filtered_BFS_roots)


    filtered_BFS_roots_indices = list(map(lambda x: BFS_root_to_index[x], filtered_BFS_roots))
    print(filtered_BFS_roots_indices)

    filtered_BFS_roots_indices_to_partition_label = {}

    for label, root_i in enumerate(filtered_BFS_roots_indices):
        filtered_BFS_roots_indices_to_partition_label[root_i] = label

    print(filtered_BFS_roots_indices_to_partition_label)

    vertex_to_partition_label = {}





    for vertex in list(G.nodes):

        # if vertex <30:
        #     print(vertex_to_distances[vertex])

        closest_BFS_root_idx = None
        distance_to_closest_BFS_root = float('inf')

        for BFS_i,BFS_distance in enumerate(vertex_to_distances[vertex]):
            if BFS_i in filtered_BFS_roots_indices:         # considering only filtered BFS s
                if BFS_distance < distance_to_closest_BFS_root:
                    closest_BFS_root_idx = BFS_i
                    distance_to_closest_BFS_root = BFS_distance


        vertex_to_partition_label[vertex] = filtered_BFS_roots_indices_to_partition_label[closest_BFS_root_idx]

    assert len(vertex_to_partition_label) == G.number_of_nodes(), "error: some nodes have been left out of partitioning"

    

    return vertex_to_partition_label


def get_parititions_by_oversampling_early_stop_BFS(G: nx.Graph , partition_count: int) -> typing.Dict[int,int]:
    additional_sample_count = 5 * math.ceil(math.log2(partition_count))
    sample_count = 10 * partition_count

    all_vertices = list(G.nodes)
    BFS_roots= [all_vertices[random.randint(0,G.number_of_nodes())] for _ in range(sample_count)]       # sampling at uniform

    vertex_to_distances = {}        # each vertex has a dict, storing distances to reached BFS s

    for v in all_vertices:
        vertex_to_distances[v] = {}

    BFS_partition_growing_status = [True for _ in range(sample_count)]       # growing stops when reahced threshold

    grow_stop_threshold = 3* int(G.number_of_nodes()/partition_count)

    BFS_partition_to_visited = [set() for _ in range(sample_count)]
    BFS_partition_to_frontier = [set() for _ in range(sample_count)]

    BFS_partition_to_current_distance = [1 for _ in range(sample_count)]

    #initial seeds
    for BFS_i, seed in enumerate(BFS_roots):
        vertex_to_distances[seed][BFS_i] = 0
        BFS_partition_to_frontier[BFS_i].add(seed)
        BFS_partition_to_visited[BFS_i].add(seed)


    while any(BFS_partition_growing_status):
        for BFS_i in range(sample_count):
            if not BFS_partition_growing_status[BFS_i]:
                continue
            new_frontier = set()
            for curr_frontier_elem in BFS_partition_to_frontier[BFS_i]:
                    for neigh in G.neighbors(curr_frontier_elem):
                        if neigh not in BFS_partition_to_visited[BFS_i]:       # if not visited by this BFS_i earlier
                            BFS_partition_to_visited[BFS_i].add(neigh)       # then add to this BFS_i
                            new_frontier.add(neigh)  
                            vertex_to_distances[neigh][BFS_i] = BFS_partition_to_current_distance[BFS_i]
            BFS_partition_to_current_distance[BFS_i]+=1
            BFS_partition_to_frontier[BFS_i] = new_frontier
            if len(BFS_partition_to_visited[BFS_i]) >= grow_stop_threshold or len(new_frontier) == 0:
                BFS_partition_growing_status[BFS_i] = False



    vertex_to_partition = {}        # current best allocation
    partition_to_vertices = [set() for _ in range(sample_count)]

    for v in all_vertices:
        selected_partition = None
        selected_partition_distance = float('inf')

        for reached_partition in vertex_to_distances[v]:
            if vertex_to_distances[v][reached_partition] < selected_partition_distance:
                selected_partition = reached_partition
                selected_partition_distance = vertex_to_distances[v][reached_partition]
        assert selected_partition != None , f"no BFS reached vertex {v}"
        vertex_to_partition[v] = selected_partition
        partition_to_vertices[selected_partition].add(v)


    deleted_partitions = set()

    while len(deleted_partitions) < sample_count - partition_count:
        smallest_partition = None
        smallest_partition_size = float('inf')

        for part_i in range(sample_count):
            if part_i in deleted_partitions:
                continue
            if len(partition_to_vertices[part_i]) < smallest_partition_size:
                smallest_partition_size = len(partition_to_vertices[part_i])
                smallest_partition = part_i

        partition_to_delete = smallest_partition

        deleted_partitions.add(partition_to_delete)

        # reassiging vertices

        for v in partition_to_vertices[partition_to_delete]:
            new_partition = None
            new_partition_distance = float('inf')

            for reached_partition in vertex_to_distances[v]:
                if reached_partition in deleted_partitions:         # not reassiging to already deleted partitions
                    continue
                if vertex_to_distances[v][reached_partition] < new_partition_distance:
                    new_partition = reached_partition
                    new_partition_distance = vertex_to_distances[v][reached_partition]
            assert new_partition != None , f"failed to reassign vertex {v}"
            vertex_to_partition[v] = new_partition
            partition_to_vertices[new_partition].add(v)
            
        partition_to_vertices[partition_to_delete] = set()    


    remaining_partitions = list(set(range(sample_count)).difference(deleted_partitions))

    assert len(remaining_partitions) == partition_count, "remaining parititons count not equal to expected partition count"

    partition_new_labels = {}

    for new_label, part_i in enumerate(remaining_partitions):       # we prefer partitions to be labeled 0,1,....partition_count-1
        partition_new_labels[part_i] = new_label

    for v in all_vertices:
        vertex_to_partition[v] = partition_new_labels[vertex_to_partition[v]]


    assert len(vertex_to_partition) == G.number_of_nodes(), "some nodes are missed in partitioning"

        

    return vertex_to_partition

    

    """
    filtering partitions idea:

    in each vertex, keep current belonged partition

    pick the smallest partition
    remove that and assign its vertices to next best parttion (traverse through its distances list and select next best based on next minimum distance)
    """


"""
assumes partitions are labelled 0,1,2,...., (partition_count-1)
"""
def get_partition_cut_sizes(G: nx.graph, vertex_to_partition: typing.Dict[int,int], partition_count) -> typing.Dict[int,int]:
    partition_cuts  = {}

    for p_i in range(partition_count):
        partition_cuts[p_i] = 0



    for u,v in G.edges:
        part_u = vertex_to_partition[u]
        part_v = vertex_to_partition[v]
        if part_u != part_v:
            partition_cuts[part_u]+=1
            partition_cuts[part_v]+=1
    
    return partition_cuts


def get_ranked_values(input_data: typing.Dict[int,int]) -> typing.Dict[int,int]:
    vals = list(input_data.values())

    vals_sorted = sorted(vals)
    val_to_rank = {}

    for rank, val in enumerate(vals_sorted):
        val_to_rank[val] = rank

    output = {}

    for key in input_data:
        output[key] =  val_to_rank[input_data[key]]

    return output

def get_standardized_values(input_data: typing.Dict[int,int]) -> typing.Dict[int,int]:
    vals = list(input_data.values())
    mean = statistics.mean(vals)
    std_dev = statistics.stdev(vals)

    output = {}

    for key in input_data:
        output[key] =  (input_data[key] - mean)/std_dev

    return output

def get_parititions_by_oversampling_remove_high_cut_partition(G: nx.Graph , partition_count: int) -> typing.Dict[int,int]:
    # additional_sample_count = 3 * math.ceil(math.log2(partition_count))
    sample_count = 2 * partition_count

    all_vertices = list(G.nodes)
    BFS_roots= [all_vertices[random.randint(0,G.number_of_nodes())] for _ in range(sample_count)]       # sampling at uniform

    vertex_to_distances = {}        # each vertex has a dict, storing distances to reached BFS s

    for v in all_vertices:
        vertex_to_distances[v] = {}

    BFS_partition_growing_status = [True for _ in range(sample_count)]       # growing stops when reahced threshold

    # grow_stop_threshold = 10* int(G.number_of_nodes()/partition_count)
    grow_stop_threshold = int(G.number_of_nodes())


    BFS_partition_to_visited = [set() for _ in range(sample_count)]
    BFS_partition_to_frontier = [set() for _ in range(sample_count)]

    BFS_partition_to_current_distance = [1 for _ in range(sample_count)]

    #initial seeds
    for BFS_i, seed in enumerate(BFS_roots):
        vertex_to_distances[seed][BFS_i] = 0
        BFS_partition_to_frontier[BFS_i].add(seed)
        BFS_partition_to_visited[BFS_i].add(seed)


    while any(BFS_partition_growing_status):
        for BFS_i in range(sample_count):
            if not BFS_partition_growing_status[BFS_i]:
                continue
            new_frontier = set()
            for curr_frontier_elem in BFS_partition_to_frontier[BFS_i]:
                    for neigh in G.neighbors(curr_frontier_elem):
                        if neigh not in BFS_partition_to_visited[BFS_i]:       # if not visited by this BFS_i earlier
                            BFS_partition_to_visited[BFS_i].add(neigh)       # then add to this BFS_i
                            new_frontier.add(neigh)  
                            vertex_to_distances[neigh][BFS_i] = BFS_partition_to_current_distance[BFS_i]
            BFS_partition_to_current_distance[BFS_i]+=1
            BFS_partition_to_frontier[BFS_i] = new_frontier
            if len(BFS_partition_to_visited[BFS_i]) >= grow_stop_threshold or len(new_frontier) == 0:
                BFS_partition_growing_status[BFS_i] = False



    vertex_to_partition = {}        # current best allocation
    partition_to_vertices = [set() for _ in range(sample_count)]

    for v in all_vertices:
        selected_partition = None
        selected_partition_distance = float('inf')

        for reached_partition in vertex_to_distances[v]:
            if vertex_to_distances[v][reached_partition] < selected_partition_distance:
                selected_partition = reached_partition
                selected_partition_distance = vertex_to_distances[v][reached_partition]
        assert selected_partition != None , f"no BFS reached vertex {v}"
        vertex_to_partition[v] = selected_partition
        partition_to_vertices[selected_partition].add(v)


    deleted_partitions = set()

    while len(deleted_partitions) < sample_count - partition_count:
        # FIXME: mean calculation can be incorrect with existing 0 values for cuts and sizes
        current_partition_cuts = get_partition_cut_sizes(G,vertex_to_partition, sample_count)
        current_partition_cuts_std = get_standardized_values(current_partition_cuts)


        partition_to_vertices_count = {}
        for p_i, vertex_set in enumerate(partition_to_vertices):
            partition_to_vertices_count[p_i] = len(vertex_set)
        
        partition_to_vertices_count_std = get_standardized_values(partition_to_vertices_count)


        partition_to_delete = None
        partition_to_delete_metric = float('-inf')

        for part_i in range(sample_count):
            if part_i in deleted_partitions:
                continue
            part_i_metric = abs(partition_to_vertices_count_std[part_i]) +current_partition_cuts_std[part_i]
            if part_i_metric > partition_to_delete_metric:
                partition_to_delete_metric = part_i_metric
                partition_to_delete = part_i

        # partition_to_delete = highest_cut_partition

        deleted_partitions.add(partition_to_delete)

        # reassiging vertices

        for v in partition_to_vertices[partition_to_delete]:
            new_partition = None
            new_partition_distance = float('inf')

            for reached_partition in vertex_to_distances[v]:
                if reached_partition in deleted_partitions:         # not reassiging to already deleted partitions
                    continue
                if vertex_to_distances[v][reached_partition] < new_partition_distance:
                    new_partition = reached_partition
                    new_partition_distance = vertex_to_distances[v][reached_partition]
            assert new_partition != None , f"failed to reassign vertex {v}"
            vertex_to_partition[v] = new_partition
            partition_to_vertices[new_partition].add(v)
            
        partition_to_vertices[partition_to_delete] = set()   


    remaining_partitions = list(set(range(sample_count)).difference(deleted_partitions))

    assert len(remaining_partitions) == partition_count, "remaining parititons count not equal to expected partition count"

    partition_new_labels = {}

    for new_label, part_i in enumerate(remaining_partitions):       # we prefer partitions to be labeled 0,1,....,(partition_count-1)
        partition_new_labels[part_i] = new_label

    for v in all_vertices:
        vertex_to_partition[v] = partition_new_labels[vertex_to_partition[v]]


    assert len(vertex_to_partition) == G.number_of_nodes(), "some nodes are missed in partitioning"

        

    return vertex_to_partition



def get_parititions_by_oversampling_remove_extra_once(G: nx.Graph , partition_count: int) -> typing.Dict[int,int]:
    # additional_sample_count = math.ceil(math.log2(partition_count))
    sample_count = 5* partition_count

    all_vertices = list(G.nodes)
    BFS_roots= [all_vertices[random.randint(0,G.number_of_nodes())] for _ in range(sample_count)]       # sampling at uniform

    vertex_to_distances = {}        # each vertex has a dict, storing distances to reached BFS s

    for v in all_vertices:
        vertex_to_distances[v] = {}

    BFS_partition_growing_status = [True for _ in range(sample_count)]       # growing stops when reahced threshold

    # grow_stop_threshold = 10* int(G.number_of_nodes()/partition_count)
    grow_stop_threshold = int(G.number_of_nodes())


    BFS_partition_to_visited = [set() for _ in range(sample_count)]
    BFS_partition_to_frontier = [set() for _ in range(sample_count)]

    BFS_partition_to_current_distance = [1 for _ in range(sample_count)]

    #initial seeds
    for BFS_i, seed in enumerate(BFS_roots):
        vertex_to_distances[seed][BFS_i] = 0
        BFS_partition_to_frontier[BFS_i].add(seed)
        BFS_partition_to_visited[BFS_i].add(seed)


    while any(BFS_partition_growing_status):
        for BFS_i in range(sample_count):
            if not BFS_partition_growing_status[BFS_i]:
                continue
            new_frontier = set()
            for curr_frontier_elem in BFS_partition_to_frontier[BFS_i]:
                    for neigh in G.neighbors(curr_frontier_elem):
                        if neigh not in BFS_partition_to_visited[BFS_i]:       # if not visited by this BFS_i earlier
                            BFS_partition_to_visited[BFS_i].add(neigh)       # then add to this BFS_i
                            new_frontier.add(neigh)  
                            vertex_to_distances[neigh][BFS_i] = BFS_partition_to_current_distance[BFS_i]
            BFS_partition_to_current_distance[BFS_i]+=1
            BFS_partition_to_frontier[BFS_i] = new_frontier
            if len(BFS_partition_to_visited[BFS_i]) >= grow_stop_threshold or len(new_frontier) == 0:
                BFS_partition_growing_status[BFS_i] = False



    vertex_to_partition = {}        # current best allocation
    partition_to_vertices = [set() for _ in range(sample_count)]

    for v in all_vertices:
        selected_partition = None
        selected_partition_distance = float('inf')

        for reached_partition in vertex_to_distances[v]:
            if vertex_to_distances[v][reached_partition] < selected_partition_distance:
                selected_partition = reached_partition
                selected_partition_distance = vertex_to_distances[v][reached_partition]
        assert selected_partition != None , f"no BFS reached vertex {v}"
        vertex_to_partition[v] = selected_partition
        partition_to_vertices[selected_partition].add(v)

    current_partition_cuts = get_partition_cut_sizes(G,vertex_to_partition, sample_count)
    current_partition_cuts_size_ratio_array = [[part_i,part_cuts/len(partition_to_vertices[part_i])] for part_i,part_cuts in current_partition_cuts.items()]
    current_partition_cuts_size_ratio_array_sorted = sorted(current_partition_cuts_size_ratio_array,key=lambda x: x[1], reverse=True) # sort by ratio
    print(current_partition_cuts_size_ratio_array_sorted)
    partitions_to_delete = [part_i_ratio_pair[0] for part_i_ratio_pair in current_partition_cuts_size_ratio_array_sorted[:(sample_count - partition_count)]]

    for part_i_to_delete in partitions_to_delete:
        for v in partition_to_vertices[part_i_to_delete]:
            new_partition = None
            new_partition_distance = float('inf')

            for reached_partition in vertex_to_distances[v]:
                if reached_partition in partitions_to_delete:         # not reassiging to deleted partitions
                    continue
                if vertex_to_distances[v][reached_partition] < new_partition_distance:
                    new_partition = reached_partition
                    new_partition_distance = vertex_to_distances[v][reached_partition]
            assert new_partition != None , f"failed to reassign vertex {v}"
            vertex_to_partition[v] = new_partition
            partition_to_vertices[new_partition].add(v)
            
        partition_to_vertices[part_i_to_delete] = set()   


    remaining_partitions = list(set(range(sample_count)).difference(set(partitions_to_delete)))

    assert len(remaining_partitions) == partition_count, "remaining parititons count not equal to expected partition count"

    partition_new_labels = {}

    for new_label, part_i in enumerate(remaining_partitions):       # we prefer partitions to be labeled 0,1,....,(partition_count-1)
        partition_new_labels[part_i] = new_label

    for v in all_vertices:
        vertex_to_partition[v] = partition_new_labels[vertex_to_partition[v]]


    assert len(vertex_to_partition) == G.number_of_nodes(), "some nodes are missed in partitioning"

        

    return vertex_to_partition


def get_partition_pairs_highest_boundary_sorted(G: nx.graph,vertex_to_partition: typing.Dict[int,int]) -> typing.List[typing.List[int]]:
    partition_pair_to_cut = {}
    for u,v in G.edges:
        part_u = vertex_to_partition[u]
        part_v = vertex_to_partition[v]
        if part_u != part_v:
            sorted_pair = tuple(sorted((part_u,part_v)))
            if sorted_pair not in partition_pair_to_cut:
                partition_pair_to_cut[sorted_pair] = 0
            partition_pair_to_cut[sorted_pair] +=1

    partition_pair_to_cut_list = [[key,pair] for key, pair in partition_pair_to_cut.items()]
    partition_pair_to_cut_list_sorted = sorted(partition_pair_to_cut_list, key=lambda x: x[1], reverse=True)    # pairs sorted in boundary highest to lowest
    return partition_pair_to_cut_list_sorted

def get_parititions_by_oversampling_merge_high_cut_pair(G: nx.Graph , partition_count: int) -> typing.Dict[int,int]:
    optimal_partition_size = G.number_of_nodes()//partition_count
    additional_sample_count = math.ceil(math.log2(partition_count))
    sample_count = additional_sample_count + partition_count

    all_vertices = list(G.nodes)
    BFS_roots= [all_vertices[random.randint(0,G.number_of_nodes())] for _ in range(sample_count)]       # sampling at uniform

    vertex_to_BFS_partition = {}
    BFS_partition_to_frontier = {}
    partition_to_vertices = {}
    BFS_partition_growing_status = [True for _ in range(sample_count)]

    #initial seeds
    for p_i, s in enumerate(BFS_roots):
        vertex_to_BFS_partition[s] = p_i
        BFS_partition_to_frontier[p_i] = set([s])
        partition_to_vertices[p_i] = set([s])

    while any(BFS_partition_growing_status):
        for p_i in range(sample_count):
            if not BFS_partition_growing_status[p_i]:
                continue
            new_frontier = set()
            for curr_frontier_elem in BFS_partition_to_frontier[p_i]:
                    for neigh in G.neighbors(curr_frontier_elem):
                        if neigh not in vertex_to_BFS_partition:       # if not assigned to any partition
                            vertex_to_BFS_partition[neigh] = p_i       # then add to current growing partition
                            partition_to_vertices[p_i].add(neigh)
                            new_frontier.add(neigh)  
            BFS_partition_to_frontier[p_i] = new_frontier
            if len(new_frontier) == 0:
                BFS_partition_growing_status[p_i] = False

    deleted_partitions = set()

    while len(deleted_partitions) < sample_count - partition_count:
        partition_pairs_with_boundary_size = get_partition_pairs_highest_boundary_sorted(G, vertex_to_BFS_partition)
        pair_to_merge = None
        for pair, boundary_size in partition_pairs_with_boundary_size:
            pair_0_size = len(partition_to_vertices[pair[0]])
            pair_1_size = len(partition_to_vertices[pair[1]])

            if (pair_0_size + pair_1_size) >= (2 * optimal_partition_size):
                continue
            else:
                pair_to_merge = pair        # pair with highest boundary
                break
        # pair_to_merge =  get_partition_pair_with_highest_boundary(G,vertex_to_BFS_partition)

        # vertices in pair_to_merge[1] will be allocated into pair_to_merge[0]
        assert pair_to_merge!=None, "no pair could be found to merge with given constraints"

        # reassiging vertices

        for v in partition_to_vertices[pair_to_merge[1]]:
            vertex_to_BFS_partition[v] = pair_to_merge[0]
            partition_to_vertices[pair_to_merge[0]].add(v)
            
        partition_to_vertices[pair_to_merge[1]] = set()   
        deleted_partitions.add(pair_to_merge[1])


    remaining_partitions = list(set(range(sample_count)).difference(deleted_partitions))

    assert len(remaining_partitions) == partition_count, "remaining parititons count not equal to expected partition count"

    partition_new_labels = {}

    for new_label, part_i in enumerate(remaining_partitions):       # we prefer partitions to be labeled 0,1,....,(partition_count-1)
        partition_new_labels[part_i] = new_label

    for v in all_vertices:
        vertex_to_BFS_partition[v] = partition_new_labels[vertex_to_BFS_partition[v]]


    assert len(vertex_to_BFS_partition) == G.number_of_nodes(), "some nodes are missed in partitioning"

        

    return vertex_to_BFS_partition


pass

"""
idea:
start with extra seeds
grow BFS until some size (2N/P) or depth (~ n^(1/dimension))

check if there are unexplored seeds
    then continue BFS until another limit (probably half the limit of the prev limit)
    to implement this have a function like growBFS(G, seed_count, current_frontiers, current_visited)

now consider the current partition allocation (possibly parititon count > p)
remove the partition with
    1. either highest cut
    2. or highest vriance to the distances to its boundary
ressign correspoding vertices to the next best partition
    
After removing check whether any unexplored vertices,

if any, call growBFS again

repeat partition removing until we get 
"""