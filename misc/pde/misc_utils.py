import networkx as nx

def part_counts(n, p):
    # Calculate the base size and remainder
    base_size = n // p
    remainder = n % p
    
    parts = []
    
    # First 'remainder' parts get base_size + 1
    for i in range(remainder):
        parts.append(base_size + 1)
    
    # Remaining parts get base_size
    for i in range(p - remainder):
        parts.append(base_size)
    
    return parts

def print_metrics(graph: nx.Graph, partition_labels, elem_to_idx_mapping, p_count):
    partition_sizes = [0 for _ in range(p_count)]
    partition_boundaries = [0 for _ in range(p_count)]
    for pl in partition_labels:
        partition_sizes[pl]+=1

    total_boundaries = 0

    for vertex in graph.nodes:
        partition = partition_labels[elem_to_idx_mapping[vertex]]
        for neigh in graph.neighbors(vertex):
            neigh_partition =  partition_labels[elem_to_idx_mapping[neigh]]
            if neigh_partition != partition:
                partition_boundaries[partition]+=1
                total_boundaries+=1
                break

    rho_max = max(partition_sizes)/(graph.number_of_nodes()/p_count)
    rho_min = min(partition_sizes)/(graph.number_of_nodes()/p_count)
    boundary_ratio = total_boundaries/graph.number_of_nodes()
    bdry_max = max(partition_boundaries)
    output = f"""
        boundary_ratio: {boundary_ratio}\t= {total_boundaries}/{graph.number_of_nodes()}
        bdry_max: {bdry_max}
        rho_max': {rho_max}\t= {max(partition_sizes)}/{int(graph.number_of_nodes()/p_count)}
        rho_min': {rho_min}\t= {min(partition_sizes)}/{int(graph.number_of_nodes()/p_count)}"""
    print(output)
