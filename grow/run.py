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
# from mayavi import mlab
from operator import itemgetter
import random
import copy
from collections import defaultdict 

import pprint
# from mpl_toolkits.mplot3d import Axes3D
import functools
import math
from datetime import datetime
import time

from hilbertcurve.hilbertcurve import HilbertCurve


from BFS_Partition import get_BFS_partitions, get_local_BFS_partitions, get_BFS_partitions_rand_seeds
from grow_partition import get_grow_partitions_ordered_BFS, get_grow_partitions_noised_BFS, get_grow_partitions_2_passes_for_size_ratio, get_grow_partitions_2_passes_oversampled
# from graph_walk import get_seeds_with_walk

# from parititions_by_oversampling import get_parititions_by_oversampling_all_graph_BFS, get_parititions_by_oversampling_early_stop_BFS, get_parititions_by_oversampling_remove_high_cut_partition, get_parititions_by_oversampling_remove_extra_once, get_parititions_by_oversampling_merge_high_cut_pair
from parititions_by_oversampling import get_parititions_by_oversampling_merge_high_cut_pair

from vtk_utils import export_points_to_vtk

# %matplotlib widget
method_names = ['SFC_morton','BFS','BFS_grow','METIS']


folder = r'/home/budvin/research/Partitioning/Meshes/10k_tet/*.mesh'
gmsh.initialize() #sys.argv)

stop_after = 70

partition_count=50

delete_vtk_files_after_viewing = True

out_file_name = datetime.now().strftime('%Y-%m-%d___%H-%M-%S-sfc_seed-2pass_2xoversample-np800-70largetetmeshes')

# out_file_name = "largest_mesh"

is_viz_only = True         # set is_viz_only = True for vizualizing partitioning for 1 file

def get_metrics(p_count, parition_labels, graph, elem_to_idx_mapping):
    partition_sizes = [0 for _ in range(p_count)]
    partition_boundaries = [0 for _ in range(p_count)]
    for cl in parition_labels:
        partition_sizes[cl]+=1

    total_boundaries = 0

    for vertex in graph.nodes:
        partition = parition_labels[elem_to_idx_mapping[vertex]]
        for neigh in G.neighbors(vertex):
            neigh_partition =  parition_labels[elem_to_idx_mapping[neigh]]
            if neigh_partition != partition:
                partition_boundaries[partition]+=1
                total_boundaries+=1
                break

    rho_max = max(partition_sizes)/(graph.number_of_nodes()/p_count)
    rho_min = min(partition_sizes)/(graph.number_of_nodes()/p_count)
    boundary_ratio = total_boundaries/graph.number_of_nodes()
    return {
        'boundary_ratio': boundary_ratio,
        'boundary_ratio_expr':f"{total_boundaries}/{graph.number_of_nodes()}",
        'rho_max': rho_max,
        'rho_max_expr': f"{max(partition_sizes)}/{int(graph.number_of_nodes()/p_count)}",
        'rho_min': rho_min,
        'rho_min_expr': f"{min(partition_sizes)}/{int(graph.number_of_nodes()/p_count)}",
        'partition_sizes': partition_sizes,
        'partition_boundaries': partition_boundaries
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

def get_stretched_increment(partition_current_size, frontier_current_size):
    # return math.log2(partition_current_size) + frontier_current_size
    # return partition_current_size
    return math.log2(partition_current_size) + 1
    # return math.sqrt(partition_current_size)

# %%
def cluster_to_color(ci):
    random.seed(ci**3)
    return (random.randint(0,255),random.randint(0,255),random.randint(0,255),125)


# %%
fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh'    # smallest
# fname = '../../Meshes/10k_tet/136935_sf_hexa.mesh_3592_12718.obj.mesh'
# fname = '../../Meshes/10k_tet/919984_sf_hexa.mesh_3960_15443.obj.mesh'
# fname = '../../Meshes/10k_tet/86233_sf_hexa.mesh_4136_15060.obj.mesh'
# fname = '../../Meshes/10k_tet/40363_sf_hexa.mesh_4618_15439.obj.mesh'
# fname_ = '../../Meshes/10k_tet/311329_sf_hexa.mesh_3800_14680.obj.mesh'
# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/104512_sf_hexa.mesh_5408_18867.obj.mesh'     # disconnected

# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/57181_sf_hexa.mesh_5006_17194.obj.mesh'
# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/42836_sf_hexa.mesh_4992_16853.obj.mesh'
# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/135214_sf_hexa.mesh_4224_14478.obj.mesh'
# fname_ = '/home/budvin/research/Partitioning/Meshes/10k_tet/168074_sf_hexa.mesh_5604_19453.obj.mesh'

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/68509_sf_hexa.mesh_4744_20488.obj.mesh"

fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/129930_sf_hexa.mesh_5292_20631.obj.mesh"            # dumbell

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/42836_sf_hexa.mesh_4992_16853.obj.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/coil.msh"       # coil mesh



# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/68645_sf_hexa.mesh_5906_20078.obj.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/97942_sf_hexa.mesh_4696_18792.obj.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/81109_sf_hexa.mesh_35542_125748.obj.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/69987_sf_hexa.mesh_7272_25474.obj.mesh"     # presentation mesh 2

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/40845_sf_hexa.mesh_4936_17372.obj.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/79183_sf_hexa.mesh_4168_15540.obj.mesh"


# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/97596_sf_hexa.mesh_2696_9392.obj.mesh"      #mesh3

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/67923_sf_hexa.mesh_2992_10000.obj.mesh"    #mesh4


fname_ = "/home/budvin/research/Partitioning/Meshes/10k_tet/69220_sf_hexa.mesh_37260_195498.obj.mesh"       # a very large mesh with holes

# fname_ = "/home/budvin/research/Partitioning/mesh_generator/tetMesh3D_lvl0.msh"         # a regular tet mesh of a cube

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_hex/1080516_sf_hexa.mesh"       # 17040 hex elements
# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_hex/55657_sf_hexa.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_hex/97596_sf_hexa.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_hex/395607_sf_hexa.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_hex/40666_sf_hexa.mesh"

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_hex/75651_sf_hexa.mesh"     #largest hex mesh 258038 elements

# fname_ = "/home/budvin/research/Partitioning/Meshes/10k_hex/60222_sf_hexa.mesh"

fname_ = '/home/budvin/research/Partitioning/Meshes/10k_hex/51140_sf_hexa.mesh'

fname_ = '/home/budvin/research/Partitioning/Meshes/10k_hex/472091_sf_hexa.mesh'

# fname_ = "/home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh"

list_of_files = filter(os.path.isfile, glob.glob(folder) ) 
  
sorted_list_of_files = sorted( list_of_files, 
                        key =  lambda x: os.stat(x).st_size) 

file_list = [
    "/home/budvin/research/Partitioning/Meshes/10k_hex/75651_sf_hexa.mesh",
    "/home/budvin/research/Partitioning/Meshes/10k_hex/60222_sf_hexa.mesh",
    '/home/budvin/research/Partitioning/Meshes/10k_hex/1356639_sf_hexa.mesh',
    "/home/budvin/research/Partitioning/Meshes/10k_hex/208721_sf_hexa.mesh",
    "/home/budvin/research/Partitioning/Meshes/10k_hex/69930_sf_hexa.mesh",
    "/home/budvin/research/Partitioning/Meshes/10k_hex/98660_sf_hexa.mesh",
    "/home/budvin/research/Partitioning/Meshes/10k_hex/66485_sf_hexa.mesh"
]

# for fname in glob.glob(folder):
file_count = 0
# for fname in glob.glob(folder):
# for fname in [sorted_list_of_files[-1]]:
# for fname in sorted_list_of_files[2800:]:
# for fname in [file_list[1]]:
# for fname in sorted_list_of_files:
for fname in [fname_]:


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
        print("hex mesh", fname)
        # continue        
        # TODO: support hex mesh with correct volume calculation
        mesh_type = 'hex'
        fn = 4
        en = 24
        vertices_n = 8
        elems, _ = gmsh.model.mesh.getElementsByType(5)
        faces    = gmsh.model.mesh.getElementFaceNodes(5, 4)
        _, nodeCoords, _ = gmsh.model.mesh.getNodesByElementType(5,returnParametricCoord=False)
    else:
        print("ignoring other type mesh", fname)
        continue


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
        if mesh_type=='tet':
            volume = abs(np.linalg.det(np.array(coord_matrix_for_volume)))/6             # volume calculation
        else: # assuming 'hex' mesh

            # Hexahedron:   (taken from https://gmsh.info/doc/texinfo/gmsh.html#Node-ordering)          
            # 
            #        v
            # 3----------2            
            # |\     ^   |\         
            # | \    |   | \         
            # |  \   |   |  \        
            # |   7------+---6       
            # |   |  +-- |-- | -> u  
            # 0---+---\--1   |       
            #  \  |    \  \  |       
            #   \ |     \  \ |       
            #    \|      w  \|        
            #     4----------5
            #     
            # splitting into 5 tetrahedra for volume calculation
            split_tets_nodes = [[0,2,3,7], [0,4,5,7] ,[2,5,6,7], [0,1,2,5], [0,2,5,7]]
            volume = 0
            for split_tet in split_tets_nodes:
                tet_matrix = np.array([coord_matrix_for_volume[split_tet[0]], coord_matrix_for_volume[split_tet[1]], coord_matrix_for_volume[split_tet[2]], coord_matrix_for_volume[split_tet[3]]])
                tet_vol = abs(np.linalg.det(tet_matrix))/6
                volume+=tet_vol

        elemVolumes[elem_idx] = volume


    # %%




    integer_elemCenterCoordsXYZ = to_integer_coords(elemCenterCoordsXYZ)
    integer_elemCenterCoordsXYZ_with_index = [ [d,i] for (i,d) in enumerate(integer_elemCenterCoordsXYZ)]

    integer_elemCenterCoordsXYZ_morton_ordered_with_index = sorted(integer_elemCenterCoordsXYZ_with_index, key=functools.cmp_to_key(morton_compare))
    morton_order = [d[1] for d in integer_elemCenterCoordsXYZ_morton_ordered_with_index]


    # %%
    # # hilbert sort
    # hilbert_curve = HilbertCurve(28, 3)
    # integer_elemCenterCoordsXYZ_hilbert_ordered_with_index = sorted(integer_elemCenterCoordsXYZ_with_index, key=lambda p: hilbert_curve.distance_from_point(p[0]))
    # hilbert_order = [d[1] for d in integer_elemCenterCoordsXYZ_hilbert_ordered_with_index]


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
    assert(None not in morton_sfc_partition_labels)
    print("SFC partitoning done")

    center_element_indices = [None for _ in range(partition_count)]

    # # calculating center elements for each partition
    # # ignoring disconnectedness
    # for p_i in range(partition_count):
    #     tot_x = 0
    #     tot_y = 0
    #     tot_z = 0
    #     tot_vol = 0
    #     for elem_i in partition_to_elem_idx[p_i]:
    #         # vol = elemVolumes[elem_i]         # considering the element volume
    #         vol=1       # not considering the element volume
    #         # vol = 1/elemVolumes[elem_i]
    #         tot_vol+=vol
    #         tot_x += (vol*elemCenterCoordsXYZ[elem_i][0])
    #         tot_y += (vol*elemCenterCoordsXYZ[elem_i][1])
    #         tot_z += (vol*elemCenterCoordsXYZ[elem_i][2])
    #     center_x = tot_x/tot_vol
    #     center_y = tot_y/tot_vol
    #     center_z = tot_z/tot_vol


    #     center_elem_idx = None
    #     center_elem_diff_to_center = math.inf

    #     # getting the closest element to partition center
    #     for elem_i in partition_to_elem_idx[p_i]:
    #         dist_sqr = 0
    #         dist_sqr += (elemCenterCoordsXYZ[elem_i][0] - center_x)**2
    #         dist_sqr += (elemCenterCoordsXYZ[elem_i][1] - center_y)**2
    #         dist_sqr += (elemCenterCoordsXYZ[elem_i][2] - center_z)**2
    #         if dist_sqr < center_elem_diff_to_center:
    #             center_elem_idx = elem_i
    #             center_elem_diff_to_center = dist_sqr
        
    #     center_element_indices[p_i] = center_elem_idx

    # median element in morton order as seed
    for p_i in range(partition_count):
        partition_size = len(partition_to_elem_idx[p_i])
        center_element_indices[p_i] = partition_to_elem_idx[p_i][partition_size//2]


    # %%

    # BFS from partition centers
    # ==================================

    # BFS_random_seed_indices = [random.randint(0, len(elems)) for _ in range(partition_count)]

    start_time = time.perf_counter()
    element_to_BFS_partition = get_BFS_partitions(G,[idx_to_element[c] for c in center_element_indices],partition_count)
    end_time = time.perf_counter()



    BFS_partition_labels = [None for _ in range(len(elems))]

    for elem in element_to_BFS_partition:
        BFS_partition_labels[element_to_idx[elem]] = element_to_BFS_partition[elem]

    assert(None not in BFS_partition_labels)
    print("BFS partitioning done")
    print(f'BFS took\t: {(end_time-start_time):.6f} seconds')

    # %%
    # seeds_by_walk = get_seeds_with_walk(G,partition_count)
    # BFS_random_seed_indices = [random.randint(0, len(elems)) for _ in range(partition_count)]
    # element_to_BFS_grow_partition = get_parititions_by_oversampling_merge_high_cut_pair(G, partition_count)
    
    # element_to_BFS_grow_partition = get_grow_partitions_2_passes_for_size_ratio(G,[idx_to_element[c] for c in center_element_indices],partition_count)
    # element_to_BFS_grow_partition = get_BFS_partitions_rand_seeds(G,partition_count)
    
    oversampled_indices = []
    oversample_rate = 2
    for sample_i in range(oversample_rate * partition_count):
        oversampled_indices.append(morton_order[sample_i* (len(elemCenterCoordsXYZ)//(oversample_rate * partition_count))])

    element_to_BFS_grow_partition = get_grow_partitions_2_passes_oversampled(G,[idx_to_element[c] for c in oversampled_indices],partition_count)

    
    BFS_grow_partition_labels = [None for _ in range(len(elems))]

    for elem in element_to_BFS_grow_partition:
        BFS_grow_partition_labels[element_to_idx[elem]] = element_to_BFS_grow_partition[elem]

    assert None not in BFS_grow_partition_labels
    print("grow partitioning done")

    


    # %%

    # METIS
    start_time = time.perf_counter()
    (edgecuts, parts) = metis.part_graph(G, partition_count)
    end_time = time.perf_counter()
    METIS_partition_labels = [lbl for lbl in parts]

    assert(None not in METIS_partition_labels)

    print("METIS done")
    print(f'METIS took\t: {(end_time-start_time):.6f} seconds')

    # %%



    SFC_metrics = get_metrics(partition_count,morton_sfc_partition_labels,G,element_to_idx)
    grow_metrics = get_metrics(partition_count,BFS_partition_labels,G,element_to_idx)
    grow_stretch_metrics = get_metrics(partition_count,BFS_grow_partition_labels,G,element_to_idx)
    metis_metrics = get_metrics(partition_count,METIS_partition_labels,G,element_to_idx)

    result_row = pd.Series()
    result_row['mesh_idx'] = file_count
    result_row['mesh_file'] = fname
    result_row['np'] = partition_count
    for metric, method_name in zip([SFC_metrics,grow_metrics,grow_stretch_metrics,metis_metrics],method_names):
        for metric_key in metric:
            result_row[f"{method_name}_{metric_key}"] = metric[metric_key]
        pass
    


    # all_results = pd.concat([all_results,pd.DataFrame([result_row])],ignore_index=True)
    if not is_viz_only:
        pd.DataFrame([result_row]).to_json(out_file_name+'.json',index=False,mode='a',lines=True,orient='records')
    else:
        pd.set_option('display.max_columns', None) 
        # with pd.option_context('display.max_columns', 150):
        print(result_row)
        break

    # out_dir_prefix="/home/budvin/research/Partitioning/paralab-partition/grow/vtk_files_3-4/"
    # vtk_file_names = []
    # for labeling, method_name in zip([morton_sfc_partition_labels,BFS_partition_labels,BFS_grow_partition_labels,METIS_partition_labels],method_names):
    #     coloring = [cluster_to_color(label) for label in labeling]
    #     vtk_file_name = f"{out_dir_prefix}{file_count}_{fname.split('/')[-1]}__np-{partition_count}___{method_name}.vtk"
    #     export_points_to_vtk(elemCenterCoordsXYZ,coloring,vtk_file_name)
    #     vtk_file_names.append(vtk_file_name)
    
    # with open(f"{out_dir_prefix}{file_count}_{fname.split('/')[-1]}__np-{partition_count}__vtkfiles.sh", 'w') as f:
    #     for method_name, vtk_file_name in zip(method_names, vtk_file_names):
    #         f.write(f'export {method_name}="{vtk_file_name}"\n')  

    file_count+=1
    print(file_count, fname, "done")



sys.stdout.flush()


if not is_viz_only:
    print(out_file_name)
    exit(0)


### viz part



# %%

vtk_file_names = []

for labeling, method_name in zip([morton_sfc_partition_labels,BFS_partition_labels,BFS_grow_partition_labels,METIS_partition_labels],method_names):
    coloring = [cluster_to_color(label) for label in labeling]
    vtk_file_name = f"{fname_.split('/')[-1]}___{method_name}.vtk"
    export_points_to_vtk(elemCenterCoordsXYZ,coloring,vtk_file_name)
    vtk_file_names.append(vtk_file_name)



# %%

import subprocess

my_env = os.environ.copy()
for i in range(len(method_names)):
    my_env[method_names[i]] = vtk_file_names[i]
subprocess.run(['/home/budvin/bin/ParaView-5.11.2-MPI-Linux-Python3.9-x86_64/bin/paraview','paraview_script.py'],env=my_env)

if delete_vtk_files_after_viewing:
    for vtk_file in vtk_file_names:
        os.remove(vtk_file)

exit(0)

# # %%

# selected_partitioning_for_viz = BFS_partition_labels

# mlab.options.backend = 'envisage'
# # s = mlab.test_plot3d()
# # mlab.figure()
# xyz = np.array(elemCenterCoordsXYZ)
# # scalar colors
# scalars = np.array([i for i in range(len(elemCenterCoordsXYZ))])
# # scalars = np.array(list(G.nodes()))
# print(scalars[:5])
# clrs = [cluster_to_color(label) for label in selected_partitioning_for_viz]
# print(clrs[:20])


# pts = mlab.points3d(
#     xyz[:, 0],
#     xyz[:, 1],
#     xyz[:, 2],
#     scalars,
#     scale_factor=0.09,
#     scale_mode="none",
#     colormap="Blues",
#     resolution=20,
# )


# pts.module_manager.scalar_lut_manager.lut.number_of_colors = partition_count+1
# pts.module_manager.scalar_lut_manager.lut.table = clrs

# # edges_adjusted_for_idx = [(element_to_idx[u],element_to_idx[v]) for u,v in G.edges()]
# # pts.mlab_source.dataset.lines = np.array(edges_adjusted_for_idx)
# # tube = mlab.pipeline.tube(pts, tube_radius=0.02)
# # mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
# scalars_centers = np.array([i for i in range(len(center_element_indices))])

# center_pts = mlab.points3d(
#     xyz[:, 0][[i for i in center_element_indices]],
#     xyz[:, 1][[i for i in center_element_indices]],
#     xyz[:, 2][[i for i in center_element_indices]],

#     scalars_centers,
#     scale_factor=0.6,
#     scale_mode="none",
#     # colormap="Blues",
#     color=(0,0,0),
#     resolution=20,
# )


# mlab.draw()
# mlab.orientation_axes()
# mlab.show()
# # %%
