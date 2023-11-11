# %%
import random
import pprint
from operator import itemgetter
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mayavi import mlab
import functools
import math






# %%
N=40
data_points = []
rng = 20
partition_count=5
# for i in range(N):
#     for j in range(N):
#         for k in range(N):
#             data_points.append([i,j,k])


for i in range(N):
    for j in range(N):
        # x = random.uniform(-rng,rng)
        x=1
        # x = random.randint(0,rng)
        # y=1
        y = random.uniform(-rng,14*rng)
        # y = random.randint(0,rng)
        # z = random.randint(0,rng)
        z = random.uniform(-rng,2*rng)
        # z=1
        data_points.append([x,y,z])

# %%


def to_integer_coords(data_points):
    levels = 20

    # # regular cube box method

    # bounding_box = [math.inf,-math.inf]         # -limit to +limit
    # for d in data_points:
    #     bounding_box[0] = min(min(d),bounding_box[0])
    #     bounding_box[1] = max(max(d),bounding_box[1])
    
    # leaf_node_length = (bounding_box[1] - bounding_box[0])/(2**levels)

    # integer_data_points = []

    # for d in data_points:
    #     d_new = []
    #     for dim_i in range(3):
    #         d_new.append(int((d[dim_i]-bounding_box[0])/leaf_node_length) + 1)
    #     integer_data_points.append(d_new)
    # return integer_data_points

    # # regular cube box method end


    # cube djust for each dimension

    bounding_box = [[math.inf,-math.inf],[math.inf,-math.inf],[math.inf,-math.inf]]     # x y z bounding box
    for d in data_points:
        for i in range(3):      # 3 axes
            bounding_box[i][0] = min(bounding_box[i][0],d[i])
            bounding_box[i][1] = max(bounding_box[i][1],d[i])

    leaf_node_lengths = [None,None,None]

    for dim_i in range(3):
        leaf_node_lengths[dim_i] = (bounding_box[i][1] - bounding_box[i][0])/(2**levels)
    
        integer_data_points = []

    for d in data_points:
        d_new = []
        for dim_i in range(3):
            d_new.append(int((d[dim_i]-bounding_box[dim_i][0])/leaf_node_lengths[dim_i]) + 1)
        integer_data_points.append(d_new)
    return integer_data_points


def morton_compare(p1_data, p2_data):
    p1 = p1_data[0]
    p2 = p2_data[0]
    if p1[0]==p2[0] and p1[1]==p2[1] and p1[2]==p2[2]:      # same point
        return 0
    temp_x = p1[0] ^ p2[0]
    temp_y = p1[1] ^ p2[1]
    temp_z = p1[2] ^ p2[2]

    maxC = temp_z
    yOrx = temp_y

    if (yOrx < temp_x):
        if ((temp_x ^ yOrx) >= yOrx):
            yOrx = temp_x

    if (maxC < yOrx):
        if ((maxC ^ yOrx) >= maxC):
            maxC = yOrx
        
    if (maxC == temp_z):
        return -1 if p1[2] < p2[2] else 1
    elif (maxC == temp_y):
        return -1 if p1[1] < p2[1] else 1
    else:
        return -1 if p1[0] < p2[0] else 1

print(data_points)
integer_data_points = to_integer_coords(data_points)
# print(integer_data_points)
integer_data_points_with_index = [ [d,i] for (i,d) in enumerate(integer_data_points)]
print(integer_data_points_with_index)
data_points_morton_ordered_with_index = sorted(integer_data_points_with_index, key=functools.cmp_to_key(morton_compare))
morton_order = [d[1] for d in data_points_morton_ordered_with_index]



figsize=15
G = nx.Graph()
pos={}
for (i,d) in enumerate(data_points):
    G.add_node(i)
    pos[i] = d

for i in range(len(morton_order)-1):
    G.add_edge(morton_order[i],morton_order[i+1])



# %%
def cluster_to_color(ci):
    random.seed(ci*ci)
    return (random.randint(0,255),random.randint(0,255),random.randint(0,255),255)

# %%
partition_size = (len(data_points))//partition_count
large_partition_count = len(data_points) % partition_count
morton_only_partition_labels = [-1 for _ in range(len(data_points))]

centers = []

for partition_idx in range(partition_count):

    if partition_idx < large_partition_count:
        size = partition_size+1
        offset = partition_idx*(partition_size+1)
    else:
        size = partition_size
        offset = large_partition_count*(partition_size+1) + (partition_idx-large_partition_count)*partition_size
    for i in range(offset, offset+size):
        morton_only_partition_labels[morton_order[i]] = partition_idx
        # print(morton_order[i])
# print(morton_only_partition_labels)

# %%
mlab.options.backend = 'envisage'

xyz = np.array(data_points)
scalars = np.array(list(G.nodes()))

pts = mlab.points3d(
    xyz[:, 0],
    xyz[:, 1],
    xyz[:, 2],
    scalars,
    scale_factor=1.5,
    scale_mode="none",
    colormap="Blues",
    resolution=20,
)

clrs = [cluster_to_color(label) for label in morton_only_partition_labels]
pts.module_manager.scalar_lut_manager.lut.number_of_colors = partition_count
pts.module_manager.scalar_lut_manager.lut.table = clrs

# edges_adjusted_for_idx = [(element_to_idx[u],element_to_idx[v]) for u,v in G.edges()]
pts.mlab_source.dataset.lines = np.array(G.edges())
tube = mlab.pipeline.tube(pts, tube_radius=0.02)
mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))

mlab.orientation_axes()
mlab.show()


