import networkx as nx
import matplotlib.pyplot as plt
import fast_matrix_market as fmm
from grave import plot_network

# G = nx.complete_graph(5)



graph_mat = fmm.read_array("dolphins.mtx")

print(graph_mat)

G = nx.Graph()

for i in range(61):
    G.add_node(i)

for i in range(61):
    for j in range(61):
        if graph_mat[i][j] != 0:
            G.add_edge(i,j)
nx.draw_random(G)
plt.show()
# fig, ax = plt.subplots()
# plot_network(G)

# plt.savefig("plot.png")