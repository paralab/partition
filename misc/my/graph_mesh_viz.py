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

# %matplotlib widget


gmsh.initialize() #sys.argv)

# %%

results = pd.DataFrame(columns=['num_elems', 'num_faces', 'type', 'npes', 'metis_lambda', 'metis_rho', 'metis_loads', 'metis_comms'])
index = set()
num_files = 0
n_p = 7
folder = r'Meshes/10k_tet/*.mesh'
stop_after = 200


def xyz_to_morton_index (coord):
    z,y,x = coord
    answer = 0
    for i in range(0,64//3):
        answer = answer | ((x & (1 << i)) << (2*i)) | ((y & (1 << i)) << (2*i+1)) | ((z & (1 << i)) << (2*i+2))

    return answer



def get_morton_index_3d(data_points):

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


def cluster_to_color(ci):
    random.seed(ci*ci)
    return (random.randint(0,255),random.randint(0,255),random.randint(0,255),255)

# Info    : Reading '/home/budvin/research/Partitioning/Meshes/10k_tet/34785_sf_hexa.mesh_11662_49615.obj.mesh'...

# for fname in ['../../Meshes/10k_tet/1582380_sf_hexa.mesh_2368_8512.obj.mesh']:
for fname in ['../../Meshes/10k_tet/69085_sf_hexa.mesh_4408_15032.obj.mesh']:
# for fname in ['../../Meshes/10k_tet/1450323_sf_hexa.mesh_3976_14878.obj.mesh']:

    gmsh.open(fname)

    if len(results) > 0 and gmsh.model.getCurrent() in index:
        gmsh.clear()
        continue
    
    if num_files > stop_after:
        break
    num_files += 1

    print(str(num_files) + ': Model ' + gmsh.model.getCurrent() + ' (' + str(gmsh.model.getDimension()) + 'D)')
    
    G = nx.Graph()

    entities = gmsh.model.getEntities()
    mesh_type = ''

    for e in entities:
        dim = e[0]
        tag = e[1]

    if dim != 3:
        continue
    
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
    else:
        mesh_type = 'hex'
        fn = 4
        en = 24
        vertices_n = 8
        elems, _ = gmsh.model.mesh.getElementsByType(5)
        faces    = gmsh.model.mesh.getElementFaceNodes(5, 4)
        _, nodeCoords, _ = gmsh.model.mesh.getNodesByElementType(5,returnParametricCoord=False)

    print('  Have ', len(elems), ' elements and ', len(faces), ' faces.')

    f2e = {}
    e2e = {}

    V = {}
    element_to_idx = {}
    for i,x in enumerate(elems):
        v = G.add_node(x)
        V[i] = x
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
    print(len(elems))
    print(len(nodeCoords))
    print(faces[:7])
    for a in range(2):
        for b in range(4):
            # for c in range(3):
            base = 12*a + 3*b
            print(nodeCoords[base:base+3])
    print(nodeCoords[:en*2])
    print(elems[:7])
    print(min(elems))
    print(max(elems))
    coord_values_per_elem = vertices_n*3        # for 3d
    elemCenterCoordsXYZ = [[-1,-1,-1] for _ in range(len(elems))]
    for i in range(0,len(nodeCoords),coord_values_per_elem):
        elem_idx = i//coord_values_per_elem
        x_tot = 0
        y_tot = 0
        z_tot = 0
        for j in range(coord_values_per_elem):
            if j%3 ==0:
                x_tot+=nodeCoords[i+j]
            elif j%3 ==1:
                y_tot+=nodeCoords[i+j]
            elif j%3 ==2:
                z_tot+=nodeCoords[i+j]
        elemCenterCoordsXYZ[elem_idx][0] = x_tot/vertices_n 
        elemCenterCoordsXYZ[elem_idx][1] = y_tot/vertices_n 
        elemCenterCoordsXYZ[elem_idx][2] = z_tot/vertices_n         





    morton = get_morton_index_3d(elemCenterCoordsXYZ)
    morton_indices = morton['morton_index']
    morton_order = morton['morton_order']


    sfc_partition_labels = [-1 for _ in range(len(elemCenterCoordsXYZ))]
    curr_partition = 0
    partition_size = len(elemCenterCoordsXYZ)//n_p
    for i, elem_idx in enumerate(morton_order):
        sfc_partition_labels[elem_idx] = min(i//partition_size,n_p-1)

    print(sfc_partition_labels[:10])



    mlab.options.backend = 'envisage'
    # s = mlab.test_plot3d()
    # mlab.figure()
    xyz = np.array(elemCenterCoordsXYZ)
    # scalar colors
    scalars = np.array(list(G.nodes()))

    # mlab.figure()
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

    clrs = [cluster_to_color(label) for label in sfc_partition_labels]
    pts.module_manager.scalar_lut_manager.lut.number_of_colors = n_p
    pts.module_manager.scalar_lut_manager.lut.table = clrs

    edges_adjusted_for_idx = [(element_to_idx[u],element_to_idx[v]) for u,v in G.edges()]
    pts.mlab_source.dataset.lines = np.array(edges_adjusted_for_idx)
    tube = mlab.pipeline.tube(pts, tube_radius=0.02)
    mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))

    mlab.orientation_axes()
    mlab.show()


