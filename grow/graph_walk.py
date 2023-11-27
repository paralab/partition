import networkx as nx
import typing
import random

def get_seeds_with_walk(G: nx.Graph, seed_count: int):
    all_vertices = list(G.nodes)
    walk_start_seeds= [all_vertices[random.randint(0,G.number_of_nodes())] for _ in range(seed_count)]
    final_seeds = [None for _ in range(seed_count)]

    walk_length = int((G.number_of_nodes()**(1.0/3))/seed_count)

    for walk_i in range(seed_count):
        current = walk_start_seeds[walk_i]
        for step_i in range(walk_length):
            neighbors = list(G.neighbors(current))
            current = random.choice(neighbors)
        final_seeds[walk_i] = current

    assert None not in final_seeds
    return final_seeds