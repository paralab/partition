import os
os.environ["METIS_DLL"] = "/home/budvin/research/Partitioning/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.so"


import networkx as nx
import matplotlib.pyplot as plt
import metis
import random

import pandas as pd
from collections import defaultdict 


def xy_to_morton_index (coord):
    x, y = coord
    answer = 0
    for i in range(0,64//2):
        answer = answer | ((x & (1 << i)) << (2*i)) | ((y & (1 << i)) << (2*i+1))


G = nx.Graph()
N_2_pow = 6
N = 2**N_2_pow
fig_size = 9
partition_count = 11
results = pd.DataFrame(columns=['method', 'np','lambda_expr', 'lambda','rho_expr', 'rho'])

#NxN grid
for i in range(N*N):
    G.add_node(i,pos=(i//N,i%N))

for i in range(N*N):
    x = i//N
    y = i % N
    # neighbors = [[x-1,y-1],[x-1,y],[x-1,y+1],[x,y-1],[x,y+1],[x+1,y-1],[x+1,y],[x+1,y+1]]
    neighbors = [[x-1,y],[x,y-1],[x,y+1],[x+1,y]]

    for ni in neighbors:
        if 0<=ni[0]<N and 0<=ni[1]<N:
            ni_index = N*ni[0] + ni[1]
            G.add_edge(i,ni_index)

node_count = G.number_of_nodes()
pos = nx.nx_agraph.graphviz_layout(G)
plt.figure(figsize=(fig_size,fig_size))
nx.draw_networkx(G, pos)