import typing
import networkx as nx
import copy
import numpy as np
import math

MAX_IMBALANCE = 1.2

def get_diffusion_partitions_from_seeds(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:

    return {}


"""
returns a mapping from vertices to diffusion value
"""
def diffuse_from_seeds(G: nx.Graph ,seeds: list[int],partition_count: int) -> typing.Dict[int,int]:
    vertex_to_val = {}
    for v in G.nodes:
        vertex_to_val[v] = 0
    for s in seeds:
        vertex_to_val[s] = 1
    diffusion_rounds = 100
    diffusion_coefficient = 0.2

    for d_i in range(diffusion_rounds):
        vertex_to_val_new_copy = copy.deepcopy(vertex_to_val)
        for v in G.nodes:
            d_val_dt = -diffusion_coefficient*vertex_to_val[v]*G.degree[v]
            for neigh in G.neighbors(v):
                d_val_dt += diffusion_coefficient*vertex_to_val[neigh]
            vertex_to_val_new_copy[v]+=d_val_dt
        vertex_to_val = vertex_to_val_new_copy
    return vertex_to_val


def get_part_sizes(partition_mapping: typing.Dict[int,int], partition_count: int):
    part_sizes = [0 for _ in range(partition_count)]
    for v in partition_mapping:
        part_sizes[partition_mapping[v]]+=1
    return part_sizes

"""
"""
def get_diffusion_rates(part_sizes: list[int], partition_count: int, vertex_count: int):
    
    print(part_sizes)
    ideal_size = vertex_count//partition_count
    DIFFUSION_MIN = 0.5
    DIFFUSION_MAX = 1

    # normalize = lambda data: [DIFFUSION_MIN + (DIFFUSION_MAX - DIFFUSION_MIN) * (x - min(data)) / (max(data) - min(data)) if max(data) > min(data) else 0.055 for x in data]
    # rates = [(size/ideal_size)/5 for size in part_sizes]

    # rates_ = [math.exp((size/ideal_size)**3) for size in part_sizes]
    # max_rate = max(rates_)
    # rates = [r/max_rate for r in rates_]
    rates = []

    for p in part_sizes:
        if p < 0.5*ideal_size:
            rates.append(0.01)
        elif p <= MAX_IMBALANCE*ideal_size:
            rates.append(0.2)
        else:
            rates.append(0.8)

    print(rates)

    # part_diffusion_rates = normalize(relative_part_sizes)
    return rates

def diffuse_from_partitions(G: nx.Graph , partition_mapping: typing.Dict[int,int], partition_count: int) -> typing.Dict[int,int]:

    part_sizes = get_part_sizes(partition_mapping, partition_count)
    part_diffusion_rates = get_diffusion_rates(part_sizes, partition_count, G.number_of_nodes())
    
    print(part_diffusion_rates)
    vertex_to_val = {}
    for v in G.nodes:
        vertex_to_val[v] = 1

    diffusion_rounds = 10
    # diffusion_coefficient = 0.1

    for d_i in range(diffusion_rounds):
        vertex_to_val_new_copy = copy.deepcopy(vertex_to_val)
        for v in G.nodes:
            d_val_dt = -vertex_to_val[v]*G.degree[v]*part_diffusion_rates[partition_mapping[v]]
            for neigh in G.neighbors(v):
                d_val_dt += vertex_to_val[neigh]*part_diffusion_rates[partition_mapping[neigh]]
            vertex_to_val_new_copy[v]+=d_val_dt
            vertex_to_val_new_copy[v] = min(1.0, vertex_to_val_new_copy[v])
        vertex_to_val = vertex_to_val_new_copy
    return vertex_to_val

def diffuse_one_round(G: nx.Graph,partition_mapping, partition_count, vertex_to_val: typing.Dict[int,float]) -> typing.Dict[int,float]:
    part_sizes = get_part_sizes(partition_mapping, partition_count)
    part_diffusion_rates = get_diffusion_rates(part_sizes, partition_count, G.number_of_nodes())

    vertex_to_val_new_copy = copy.deepcopy(vertex_to_val)
    for v in G.nodes:
        d_val_dt = -vertex_to_val[v]*G.degree[v]*part_diffusion_rates[partition_mapping[v]]
        for neigh in G.neighbors(v):
            d_val_dt += vertex_to_val[neigh]*part_diffusion_rates[partition_mapping[neigh]]
        vertex_to_val_new_copy[v]+=d_val_dt
        vertex_to_val_new_copy[v] = max(min(1.0, vertex_to_val_new_copy[v]), 0)
    return vertex_to_val_new_copy

# vertecies in each level should be disjoint
def diffuse_2_level_compact(G: nx.Graph,partition_mapping, partition_count, vertex_to_val: typing.Dict[int,float],
                            level_0: list[int], level_1: list[int], level_2: list[int]) -> typing.Dict[int,float]:
    part_sizes = get_part_sizes(partition_mapping, partition_count)
    part_diffusion_rates = get_diffusion_rates(part_sizes, partition_count, G.number_of_nodes())
    levels = [level_0, level_0 + level_1, level_0 + level_1 + level_2]

    for li in range(3):
        vertex_set = levels[li]
        vertex_to_val_new_copy = copy.deepcopy(vertex_to_val)
        for v in vertex_set:
            d_val_dt = -vertex_to_val[v]*G.degree[v]*part_diffusion_rates[partition_mapping[v]]
            for neigh in G.neighbors(v):
                d_val_dt += vertex_to_val[neigh]*part_diffusion_rates[partition_mapping[neigh]]
            vertex_to_val_new_copy[v]+=d_val_dt
            vertex_to_val_new_copy[v] = max(min(1.0, vertex_to_val_new_copy[v]), 0)
        vertex_to_val =  vertex_to_val_new_copy
    return vertex_to_val


def refine_by_diffusion(G: nx.Graph , partition_mapping_: typing.Dict[int,int], partition_count: int):
    partition_mapping = copy.deepcopy(partition_mapping_)
    vertex_to_val = {}
    for v in G.nodes:
        vertex_to_val[v] = 1.0
    max_diffusion_rounds = 30

    # while max(part_sizes) > max_part_size:
    for d_i in range(max_diffusion_rounds):

        vertex_to_val = diffuse_one_round(G, partition_mapping, partition_count,vertex_to_val)
        part_sizes = get_part_sizes(partition_mapping, partition_count)
        if(max(part_sizes)/(G.number_of_nodes()/partition_count) < MAX_IMBALANCE):
            break
        part_diffusion_rates = get_diffusion_rates(part_sizes, partition_count, G.number_of_nodes())

        vertex_to_val_new_copy = copy.deepcopy(vertex_to_val)
        partition_mapping_new_copy = copy.deepcopy(partition_mapping)
        # changed_v_set = set()
        # changed_v_set_neigh_ = set()
        for v in G.nodes:
            if (vertex_to_val[v] < 0.5):
                neighborhood = {}       
                incoming_flux = 0
                for neigh in G.neighbors(v):
                    neigh_label = partition_mapping[neigh]
                    if not neigh_label in neighborhood:
                        neighborhood[neigh_label] = 0
                    neighborhood[neigh_label]+=1
                    incoming_flux += vertex_to_val[neigh]*part_diffusion_rates[partition_mapping[neigh]]
                    # changed_v_set_neigh_.add(neigh)
                min_part = min(neighborhood, key=lambda k: part_diffusion_rates[k])
                vertex_to_val_new_copy[v] = incoming_flux /(G.degree(v)*part_diffusion_rates[min_part])
                vertex_to_val_new_copy[v] = max(min(1.0, vertex_to_val_new_copy[v]), 0)
                partition_mapping_new_copy[v] =  min_part
                # print(vertex_to_val[v], vertex_to_val_new_copy[v])
                # changed_v_set.add(v)
        vertex_to_val = vertex_to_val_new_copy
        partition_mapping = partition_mapping_new_copy

        # changed_v_set_neigh = changed_v_set_neigh_.difference(changed_v_set)
        # changed_v_set_neigh_neigh_ = set()
        # for neigh in changed_v_set_neigh:
        #     for neigh_neigh in G.neighbors(neigh):
        #         changed_v_set_neigh_neigh_.add(neigh_neigh)
        # changed_v_set_neigh_neigh = changed_v_set_neigh_neigh_.difference(changed_v_set_neigh).difference(changed_v_set)
        # vertex_to_val = diffuse_2_level_compact(G, partition_mapping, partition_count, vertex_to_val, 
        #                                         list(changed_v_set), list(changed_v_set_neigh), list(changed_v_set_neigh_neigh))
        # for _ in range(5):
        #     vertex_to_val = diffuse_one_round(G, partition_mapping, partition_count,vertex_to_val)

    return [partition_mapping, vertex_to_val]



def partitions_by_diffusion(G: nx.Graph ,seeds: list[int], partition_count: int) -> typing.Dict[int,int]:
    vertex_to_val = {}
    partition_mapping = {}
    for v in G.nodes:
        vertex_to_val[v] = 0.0
    for p_i, s in enumerate(seeds):
        vertex_to_val[s] = 1.0
        partition_mapping[s] = p_i
    
    diffusion_rounds_stop_guess = 200
    visited_all = False
    d_i = 0
    while d_i < diffusion_rounds_stop_guess:
        d_i+=1
        if (not visited_all) and (len(partition_mapping) == G.number_of_nodes()):
            print(f"visited all at round = {d_i}")
            visited_all = True
            # diffusion_rounds_stop_guess = int(1.5*d_i)
        part_sizes = get_part_sizes(partition_mapping, partition_count)
        if visited_all:
            max(part_sizes) <= 1.5 * np.mean(part_sizes)
            print("1.5 criteria met. Stopping...")
            break
        part_diffusion_rates = get_diffusion_rates(part_sizes, partition_count, G.number_of_nodes())
        partition_mapping_new_copy = copy.deepcopy(partition_mapping)

        vertex_to_val_new_copy = copy.deepcopy(vertex_to_val)
        for v in G.nodes:
            incoming_flux = 0
            neigh_parts = set()
            outgoing_flux = vertex_to_val[v]*G.degree[v]*part_diffusion_rates[partition_mapping[v]] if (v in partition_mapping) else 0
            for neigh in G.neighbors(v):
                if neigh in partition_mapping:
                    incoming_flux += vertex_to_val[neigh]*part_diffusion_rates[partition_mapping[neigh]]
                    neigh_parts.add(partition_mapping[neigh])
            vertex_to_val_new_copy[v]+=(incoming_flux - outgoing_flux)
            if len(neigh_parts):
                if ((v not in partition_mapping) or (vertex_to_val_new_copy[v] < 0.5)):      # first diffusion visit or alpha reduces
                    min_part = min(neigh_parts, key=lambda k: part_diffusion_rates[k])
                    vertex_to_val_new_copy[v] = incoming_flux /(G.degree(v)*part_diffusion_rates[min_part])
                    partition_mapping_new_copy[v] =  min_part

            vertex_to_val_new_copy[v] = max(min(1.0, vertex_to_val_new_copy[v]), 0)

        vertex_to_val = vertex_to_val_new_copy
        partition_mapping = partition_mapping_new_copy
    return [partition_mapping, vertex_to_val]