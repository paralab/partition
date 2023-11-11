# %%
import os
os.environ["METIS_DLL"] = "/home/budvin/research/Partitioning/metis-5.1.0/build/Linux-x86_64/libmetis/libmetis.so"


# %%
import sys
import glob

import gmsh
import networkx as nx
import metis

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab
from operator import itemgetter
import random
import copy
from collections import defaultdict 

import pprint
from mpl_toolkits.mplot3d import Axes3D
import functools
import math
from datetime import datetime

# %matplotlib widget
all_results_columns = ['mesh_file', 'np','metis_lambda', 'metis_rho_max','metis_rho_min','SFC_lambda','SFC_rho_max','SFC_rho_min','grow_lambda','grow_rho_max','grow_rho_min']
all_results = pd.DataFrame()
folder = r'/home/budvin/research/Partitioning/Meshes/10k_tet/*.mesh'
gmsh.initialize() #sys.argv)

stop_after = 70

partition_count=9


def get_metrics(p_count, parition_labels, graph, elem_to_idx_mapping):
    partition_sizes = [0 for _ in range(p_count)]
    partition_cuts = [0 for _ in range(p_count)]
    for cl in parition_labels:
        partition_sizes[cl]+=1

    edge_cuts = 0

    for u,v in graph.edges:
        part_u = parition_labels[elem_to_idx_mapping[u]]
        part_v = parition_labels[elem_to_idx_mapping[v]]
        if part_u != part_v:
            edge_cuts+=1
            partition_cuts[part_u]+=1
            partition_cuts[part_v]+=1

    rho_max = max(partition_sizes)/(graph.number_of_nodes()/p_count)
    rho_min = min(partition_sizes)/(graph.number_of_nodes()/p_count)
    lmb = edge_cuts/graph.number_of_edges()
    return {
        'lambda': lmb,
        'lambda_expr':f"{edge_cuts}/{graph.number_of_edges()}",
        'rho_max': rho_max,
        'rho_max_expr': f"{max(partition_sizes)}/{int(graph.number_of_nodes()/p_count)}",
        'rho_min': rho_min,
        'rho_min_expr': f"{min(partition_sizes)}/{int(graph.number_of_nodes()/p_count)}",
        'partition_sizes': partition_sizes,
        'partitions_cuts': partition_cuts
    }


def to_integer_coords(data_points):
    levels = 28

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

# %%

# fname = '../../Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh'
# fname = '../../Meshes/10k_tet/136935_sf_hexa.mesh_3592_12718.obj.mesh'
# fname = '../../Meshes/10k_tet/919984_sf_hexa.mesh_3960_15443.obj.mesh'
# fname = '../../Meshes/10k_tet/86233_sf_hexa.mesh_4136_15060.obj.mesh'
# fname = '../../Meshes/10k_tet/40363_sf_hexa.mesh_4618_15439.obj.mesh'
# fname = '../../Meshes/10k_tet/311329_sf_hexa.mesh_3800_14680.obj.mesh'
# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/104512_sf_hexa.mesh_5408_18867.obj.mesh'     # disconnected
# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/57181_sf_hexa.mesh_5006_17194.obj.mesh'
# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/42836_sf_hexa.mesh_4992_16853.obj.mesh'
fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/135214_sf_hexa.mesh_4224_14478.obj.mesh'


list_of_files = filter(os.path.isfile, glob.glob(folder) ) 
  
sorted_list_of_files = sorted( list_of_files, 
                        key =  lambda x: os.stat(x).st_size) 

# for fname in glob.glob(folder):
file_count = 0
# for fname in glob.glob(folder):
for fname in sorted_list_of_files:
# for fname in [fname_]:


    if file_count >= stop_after:
        break


    gmsh.open(fname)



    print(': Model ' + gmsh.model.getCurrent() + ' (' + str(gmsh.model.getDimension()) + 'D)')

    G = nx.Graph()

    entities = gmsh.model.getEntities()
    mesh_type = ''

    for e in entities:
        dim = e[0]
        tag = e[1]

    if dim != 3:
        print("not a 3D mesh\nexiting")
        exit()
        # continue

    type = gmsh.model.getType(e[0], e[1])
    name = gmsh.model.getEntityName(e[0], e[1])
    if len(name): name += ' '
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)

    if elemTypes[0] == 4:
        mesh_type = 'tet'
        fn = 3
        en = 12
        vertices_n = 4 
        elems, _ = gmsh.model.mesh.getElementsByType(4)
        faces    = gmsh.model.mesh.getElementFaceNodes(4, 3)
        _, nodeCoords, _ = gmsh.model.mesh.getNodesByElementType(4,returnParametricCoord=False)
    elif elemTypes[0] == 5:
        # continue;
        # TODO: support hex mesh with correct volume calculation
        mesh_type = 'hex'
        fn = 4
        en = 24
        vertices_n = 8
        elems, _ = gmsh.model.mesh.getElementsByType(5)
        faces    = gmsh.model.mesh.getElementFaceNodes(5, 4)
        _, nodeCoords, _ = gmsh.model.mesh.getNodesByElementType(5,returnParametricCoord=False)
    else:
        print("mesh other than tet or hex")
        # continue;


    print(mesh_type, 'mesh  has ', len(elems), ' elements and ', len(faces), ' faces.')

    f2e = {}
    e2e = {}

    idx_to_element = {}
    element_to_idx = {}
    for i,x in enumerate(elems):
        v = G.add_node(x)
        idx_to_element[i] = x
        element_to_idx[x] = i

    for i in range(0, len(faces), fn):
        f = tuple(sorted(faces[i:i+fn]))
        t = elems[i//en]
        if not f in f2e:
            f2e[f] = [t]
        else:
            f2e[f].append(t)

    # compute neighbors by face
    for i in range(0, len(faces), fn):
        f = tuple(sorted(faces[i:i+fn]))
        t = elems[i//en]
        if not t in e2e:
            e2e[t] = set()
        for tt in f2e[f]:
            if tt != t:
                e2e[t].add(tt)

    for k in e2e:
        for j in e2e[k]:
            G.add_edge(k,j)

    if not nx.is_connected(G):
        print("ignoring disconnected", fname)
        continue



    coord_values_per_elem = vertices_n*3        # for 3d
    elemCenterCoordsXYZ = [[-1,-1,-1] for _ in range(len(elems))]
    elemVolumes = [None for _ in range(len(elems))]
    for i in range(0,len(nodeCoords),coord_values_per_elem):
        elem_idx = i//coord_values_per_elem
        x_tot = 0
        y_tot = 0
        z_tot = 0
        coord_matrix_for_volume = []
        for j in range(coord_values_per_elem):
            if j%3 ==0:
                x_tot+=nodeCoords[i+j]
                coord_matrix_for_volume.append(np.concatenate((nodeCoords[i+j:i+j+3],np.array([1]))))
            elif j%3 ==1:
                y_tot+=nodeCoords[i+j]
            elif j%3 ==2:
                z_tot+=nodeCoords[i+j]
        # setting geometric center as element coordinates
        elemCenterCoordsXYZ[elem_idx][0] = x_tot/vertices_n 
        elemCenterCoordsXYZ[elem_idx][1] = y_tot/vertices_n 
        elemCenterCoordsXYZ[elem_idx][2] = z_tot/vertices_n
        # print(coord_matrix_for_volume)
        volume = np.linalg.det(np.array(coord_matrix_for_volume))/6             # volume calculation
        elemVolumes[elem_idx] = volume


    # %%




    integer_elemCenterCoordsXYZ = to_integer_coords(elemCenterCoordsXYZ)
    integer_elemCenterCoordsXYZ_with_index = [ [d,i] for (i,d) in enumerate(integer_elemCenterCoordsXYZ)]

    integer_elemCenterCoordsXYZ_morton_ordered_with_index = sorted(integer_elemCenterCoordsXYZ_with_index, key=functools.cmp_to_key(morton_compare))
    morton_order = [d[1] for d in integer_elemCenterCoordsXYZ_morton_ordered_with_index]


    # %%
    partition_size = (len(elemCenterCoordsXYZ))//partition_count
    large_partition_count = len(elemCenterCoordsXYZ) % partition_count
    morton_sfc_partition_labels = [-1 for _ in range(len(elemCenterCoordsXYZ))]
    partition_to_elem_idx = {}
    centers = []

    for partition_idx in range(partition_count):

        if partition_idx < large_partition_count:
            size = partition_size+1
            offset = partition_idx*(partition_size+1)
        else:
            size = partition_size
            offset = large_partition_count*(partition_size+1) + (partition_idx-large_partition_count)*partition_size
        partition_to_elem_idx[partition_idx] =  []
        for i in range(offset, offset+size):
            morton_sfc_partition_labels[morton_order[i]] = partition_idx
            partition_to_elem_idx[partition_idx].append(morton_order[i])
    print("SFC partitoning done")
    # calculating center elements for each partition
    # ignoring disconnectedness
    center_element_indices = [None for _ in range(partition_count)]

    for p_i in range(partition_count):
        tot_x = 0
        tot_y = 0
        tot_z = 0
        tot_vol = 0
        for elem_i in partition_to_elem_idx[p_i]:
            tot_vol+=1
            tot_x += (1*elemCenterCoordsXYZ[elem_i][0])
            tot_y += (1*elemCenterCoordsXYZ[elem_i][1])
            tot_z += (1*elemCenterCoordsXYZ[elem_i][2])
        center_x = tot_x/tot_vol
        center_y = tot_y/tot_vol
        center_z = tot_z/tot_vol


        center_elem_idx = None
        center_elem_diff_to_center = math.inf

        # getting the closest element to partition center
        for elem_i in partition_to_elem_idx[p_i]:
            dist_sqr = 0
            dist_sqr += (elemCenterCoordsXYZ[elem_i][0] - center_x)**2
            dist_sqr += (elemCenterCoordsXYZ[elem_i][1] - center_y)**2
            dist_sqr += (elemCenterCoordsXYZ[elem_i][2] - center_z)**2
            if dist_sqr < center_elem_diff_to_center:
                center_elem_idx = elem_i
                center_elem_diff_to_center = dist_sqr
        
        center_element_indices[p_i] = center_elem_idx


    # %%

    # BFS from partition centers
    # ==================================

    element_to_BFS_partition = {}

    #initial seeds
    for i, center_idx in enumerate(center_element_indices):
        element_to_BFS_partition[idx_to_element[center_idx]] = i

    not_visited = set(elems)
    while len(not_visited) > 0:
        element_to_BFS_partition_copy = copy.deepcopy(element_to_BFS_partition)
        visited_this_round = set()
        for vertex in not_visited:
            for neigh in G.neighbors(vertex):
                if neigh in element_to_BFS_partition_copy:
                    element_to_BFS_partition[vertex] = element_to_BFS_partition_copy[neigh]
                    visited_this_round.add(vertex)
                    break


        not_visited = not_visited.difference(visited_this_round)

    BFS_partition_labels = [None for _ in range(len(elems))]

    for elem in element_to_BFS_partition:
        BFS_partition_labels[element_to_idx[elem]] = element_to_BFS_partition[elem]

    print("BFS partitioning done")
    # %%

    # METIS

    (edgecuts, parts) = metis.part_graph(G, partition_count)
    METIS_partition_labels = [lbl for lbl in parts]

    print("METIS done")

    # %%



    SFC_metrics = get_metrics(partition_count,morton_sfc_partition_labels,G,element_to_idx)
    grow_metrics = get_metrics(partition_count,BFS_partition_labels,G,element_to_idx)
    metis_metrics = get_metrics(partition_count,METIS_partition_labels,G,element_to_idx)
    result_row = pd.Series()
    result_row['mesh_idx'] = file_count
    result_row['mesh_file'] = fname
    result_row['np'] = partition_count
    for metric, method_name in zip([SFC_metrics,grow_metrics,metis_metrics],['SFC','grow','metis']):
        for metric_key in metric:
            result_row[f"{method_name}_{metric_key}"] = metric[metric_key]
        pass

    all_results = pd.concat([all_results,pd.DataFrame([result_row])],ignore_index=True)

    file_count+=1
    print(file_count, fname, "done")



# print(result_row)
print(all_results)
out_file_name = datetime.now().strftime('%Y-%m-%d___%H-%M-%S')
# out_file_name = 'new_results'
all_results.to_csv(out_file_name+ '.csv',index=False)
all_results.to_json(out_file_name+'.json',index=False)


exit(0)


### viz part


# %%
def cluster_to_color(ci):
    random.seed(ci*ci)
    return (random.randint(0,255),random.randint(0,255),random.randint(0,255),125)



# %%

selected_partitioning_for_viz = morton_sfc_partition_labels

mlab.options.backend = 'envisage'
# s = mlab.test_plot3d()
# mlab.figure()
xyz = np.array(elemCenterCoordsXYZ)
# scalar colors
scalars = np.array([i for i in range(len(elemCenterCoordsXYZ))])
# scalars = np.array(list(G.nodes()))
print(scalars[:5])
clrs = [cluster_to_color(label) for label in selected_partitioning_for_viz]
print(clrs[:20])


pts = mlab.points3d(
    xyz[:, 0],
    xyz[:, 1],
    xyz[:, 2],
    scalars,
    scale_factor=0.09,
    scale_mode="none",
    colormap="Blues",
    resolution=20,
)


pts.module_manager.scalar_lut_manager.lut.number_of_colors = partition_count+1
pts.module_manager.scalar_lut_manager.lut.table = clrs

# edges_adjusted_for_idx = [(element_to_idx[u],element_to_idx[v]) for u,v in G.edges()]
# pts.mlab_source.dataset.lines = np.array(edges_adjusted_for_idx)
# tube = mlab.pipeline.tube(pts, tube_radius=0.02)
# mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
scalars_centers = np.array([i for i in range(len(center_element_indices))])

center_pts = mlab.points3d(
    xyz[:, 0][[i for i in center_element_indices]],
    xyz[:, 1][[i for i in center_element_indices]],
    xyz[:, 2][[i for i in center_element_indices]],

    scalars_centers,
    scale_factor=0.6,
    scale_mode="none",
    # colormap="Blues",
    color=(0,0,0),
    resolution=20,
)


mlab.draw()
mlab.orientation_axes()
mlab.show()
# %%
