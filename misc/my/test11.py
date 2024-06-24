# %%
import random
import pprint
from operator import itemgetter
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mayavi import mlab
import random

# %%

pp = pprint.PrettyPrinter(indent=2)
# random.seed(1)
def display(x):
    print("")
    print(x)
    # pp.pprint(x)
    print("")
# display = pp.pprint
figsize=10

# %%

def cluster_to_color(ci):
    random.seed(ci)
    return (random.randint(0,255),random.randint(0,255),random.randint(0,255),255)

def xyz_to_morton_index (coord):
    z,y,x = coord
    answer = 0
    for i in range(0,64//3):
        answer = answer | ((x & (1 << i)) << (2*i)) | ((y & (1 << i)) << (2*i+1)) | ((z & (1 << i)) << (2*i+2))

    return answer


def xy_to_morton_index (coord):
    x, y = coord
    answer = 0
    for i in range(0,32//2):
        answer = answer | ((x & (1 << i)) << (2*i)) | ((y & (1 << i)) << (2*i+1))

    return answer

def get_morton_index_3d_old(data_points):

    all_coord_values = []
    for idx, d in enumerate(data_points):
        # keep pair [value,location]
        all_coord_values.append([
            d[0],[idx,0]
        ])
        all_coord_values.append([
            d[1],[idx,1]
        ])
        all_coord_values.append([
            d[2],[idx,2]
        ])
    all_coord_values_sorted = sorted(all_coord_values, key=itemgetter(0))

    data_points_adjusted = [[-1,-1,-1] for _ in range(len(data_points))]
    for i, x in enumerate(all_coord_values_sorted):
        data_points_adjusted[x[1][0]][x[1][1]] = i
    morton_indices_temp = [[xyz_to_morton_index(data),i] for i,data in enumerate(data_points_adjusted)]


    morton_indices_temp_sorted = sorted(morton_indices_temp, key=itemgetter(0))

    morton_order = [x[1] for x in morton_indices_temp_sorted]

    morton_indices_final = [-1 for _ in range((len(data_points)))]


    for i,x in enumerate(morton_indices_temp_sorted):
        morton_indices_final[x[1]] = i
    return {'morton_index':morton_indices_final,'morton_order':morton_order}

def get_morton_index_3d(data_points):
    all_coord_values = []
    for d in data_points:
        all_coord_values.extend(d)
    all_coord_values_sorted = sorted(all_coord_values)
    min_diff = float('inf')
    for i in range(len(all_coord_values_sorted)-1):
        min_diff = min(min_diff,all_coord_values_sorted[i+1]-all_coord_values_sorted[i])
    min_coord = all_coord_values_sorted[0]
    max_coord = all_coord_values_sorted[-1]
    small_cube_count_per_dim = (max_coord - min_coord)/min_diff
    data_points_adjusted = []
    for d in data_points:
        adjusted = [None,None,None]
        for i in range(3):
            adjusted[i] = int((d[i]-min_coord)/min_diff)
        data_points_adjusted.append(adjusted)

    morton_indices_temp = [[xyz_to_morton_index(data),i] for i,data in enumerate(data_points_adjusted)]


    morton_indices_temp_sorted = sorted(morton_indices_temp, key=itemgetter(0))
    print(morton_indices_temp_sorted[0],morton_indices_temp_sorted[1])

    morton_order = [x[1] for x in morton_indices_temp_sorted]

    morton_indices_final = [-1 for _ in range((len(data_points)))]


    for i,x in enumerate(morton_indices_temp_sorted):
        morton_indices_final[x[1]] = i
    return {'morton_index':morton_indices_final,'morton_order':morton_order}  

# %%

N=10
max_coord_val = 5
n_p=2
data_points = []
# for i in range(N):
#     for j in range(N):
#         for k in range(N):
#             data_points.append([i,j,k])


for i in range(N):
    for j in range(N):
        x = random.uniform(-max_coord_val,max_coord_val)
        y = random.uniform(-max_coord_val,max_coord_val)
        z = random.uniform(-max_coord_val,max_coord_val)
        # data_points.append([x,y])
        data_points.append([x,y,z])



# display(data_points)


# all_coord_values = []
# for idx, d in enumerate(data_points):
#     # keep pair [value,location]
#     all_coord_values.append([
#         d[0],[idx,0]
#     ])
#     all_coord_values.append([
#         d[1],[idx,1]
#     ])
#     all_coord_values.append([
#         d[2],[idx,2]
#     ])
# display(all_coord_values)
# # for i in range(N*N*N):
# #     print(xyz_to_morton_index(data_points[i]))


# all_coord_values_sorted = sorted(all_coord_values, key=itemgetter(0))

# data_points_adjusted = [[-1,-1,-1] for _ in range(len(data_points))]
# for i, x in enumerate(all_coord_values_sorted):
#     # data_points_adjusted[x[1][0]][x[1][1]] = i
#     data_points_adjusted[x[1][0]][x[1][1]] = i


# morton_indices_temp = [[xyz_to_morton_index(data),i] for i,data in enumerate(data_points_adjusted)]



# display(morton_indices_temp)
# morton_indices_temp_sorted = sorted(morton_indices_temp, key=itemgetter(0))
# display(morton_indices_temp_sorted)

# morton_indices_final = [-1 for _ in range((len(data_points)))]

# for i,x in enumerate(morton_indices_temp_sorted):
#     morton_indices_final[x[1]] = i

# display(morton_indices_final)


morton = get_morton_index_3d(data_points)

morton_indices = morton['morton_index']
morton_order = morton['morton_order']

sfc_partition_labels = [-1 for _ in range(len(data_points))]
partition_size = len(data_points)//n_p
for i, elem_idx in enumerate(morton_order):
    sfc_partition_labels[elem_idx] = min(i//partition_size,n_p-1)


G = nx.Graph()

for i,_ in enumerate(data_points):
    G.add_node(i)

for i in range(len(morton_order)-1):
    G.add_edge(morton_order[i],morton_order[i+1])
    pass

morton_graph_labels = {}

for node_idx,label in enumerate(morton_indices):
    morton_graph_labels[node_idx]=label

# plt.figure(figsize=(figsize,figsize))
# nx.draw_networkx(G, data_points,labels=morton_graph_labels)

# %%
# %matplotlib widget
# node_xyz = np.array([data_points[v] for v in sorted(G)])
# edge_xyz = np.array([(data_points[u], data_points[v]) for u, v in G.edges()])

# # Create the 3D figure
# fig = plt.figure()
# # ax = Axes3D(fig)
# ax = fig.add_subplot(projection="3d")

# # Plot the nodes - alpha is scaled by "depth" automatically
# ax.scatter(*node_xyz.T, s=100, ec="w")

# # Plot the edges
# # for vizedge in edge_xyz:
# #     ax.plot(*vizedge.T, color="tab:gray")


# def _format_axes(ax):
#     """Visualization options for the 3D axes."""
#     # Turn gridlines off
#     ax.grid(False)
#     # Suppress tick labels
#     for dim in (ax.xaxis, ax.yaxis, ax.zaxis):
#         dim.set_ticks([])
#     # Set axes labels
#     ax.set_xlabel("x")
#     ax.set_ylabel("y")
#     ax.set_zlabel("z")


# _format_axes(ax)
# fig.tight_layout()
# plt.show()

# %%
# mlab.init_notebook(backend='x3d')
mlab.options.backend = 'envisage'
# s = mlab.test_plot3d()
# mlab.figure()

# %%
xyz = np.array(data_points)
# scalar colors
scalars = np.array(list(G.nodes())) + 1

# mlab.figure()
pts = mlab.points3d(
    xyz[:, 0],
    xyz[:, 1],
    xyz[:, 2],
    scalars,
    scale_factor=0.1,
    scale_mode="none",
    colormap="Blues",
    resolution=20,
)

clrs = [cluster_to_color(label) for label in sfc_partition_labels]
pts.module_manager.scalar_lut_manager.lut.number_of_colors = n_p
pts.module_manager.scalar_lut_manager.lut.table = clrs

pts.mlab_source.dataset.lines = np.array(list(G.edges()))
tube = mlab.pipeline.tube(pts, tube_radius=0.01)
mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))

mlab.orientation_axes()

mlab.show()


