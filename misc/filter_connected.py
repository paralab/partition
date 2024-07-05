# %%
import os
import sys


# %%
import glob

import gmsh
import networkx as nx

import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import numpy as np
from mayavi import mlab


# %%
GMSH_TET_ELEMENT_TYPE = 4
GMSH_HEX_ELEMENT_TYPE = 5


GMSH_TRI_FACE_TYPE = 3
GMSH_QUAD_FACE_TYPE = 4



# %%




# %%




folder = r'/home/budvin/research/Partitioning/Meshes/10k_hex/*.mesh'


stop_after = 100


# out_file_name = "largest_mesh"

is_viz_only = True         # set is_viz_only = True for vizualizing partitioning for 1 file



list_of_files = filter(os.path.isfile, glob.glob(folder) ) 
  
sorted_list_of_files = sorted( list_of_files, 
                        key =  lambda x: os.stat(x).st_size) 
# print(sorted_list_of_files[:2])
# sorted_list_of_files.reverse()

reversed_list = list(reversed(sorted_list_of_files))


# for fname in glob.glob(folder):
file_count = 0

# for fname in ['/home/budvin/research/Partitioning/mesh_generator/hex-box-5x5x2.msh']:
# for fname in ['/home/budvin/research/Partitioning/Meshes/10k_hex/1582380_sf_hexa.mesh']:


out_file = open("connected_meshes_hex.txt", "a")  

for fname in reversed_list[6:]:
    file_count+=1
    print("processing file : ", fname)


    if file_count >= stop_after:
        break


    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.open(fname)

    # print('Model ' + gmsh.model.getCurrent() + ' (' + str(gmsh.model.getDimension()) + 'D)')


    element_types = gmsh.model.mesh.getElementTypes(3)

    if len(element_types) == 0:
        print("skipping file: no 3D elements")
        continue
    elif len(element_types) > 1:
        print("skipping file: more than one 3D element type")
        continue

    element_type = element_types[0]

    if element_type == GMSH_TET_ELEMENT_TYPE:
        mesh_type = 'TET'
        face_type = GMSH_TRI_FACE_TYPE
        nodes_per_face = 3
        faces_per_elem = 4
    elif element_type == GMSH_HEX_ELEMENT_TYPE:
        mesh_type = 'HEX'
        face_type = GMSH_QUAD_FACE_TYPE
        nodes_per_face = 4
        faces_per_elem = 6
    else:
        print("skipping file: unknown element type")
        continue

    elem_tags, _ = gmsh.model.mesh.getElementsByType(element_type)
    elem_face_nodes    = gmsh.model.mesh.getElementFaceNodes(element_type, face_type)


    print(mesh_type, 'mesh  has', len(elem_tags), 'elements')


    G = nx.Graph()

    f2e = {}
    e2e = {}


    for i,x in enumerate(elem_tags):
        v = G.add_node(x)


    for i in range(0, len(elem_face_nodes), nodes_per_face):
        f = tuple(sorted(elem_face_nodes[i:i+nodes_per_face]))
        t = elem_tags[i//(nodes_per_face*faces_per_elem)]
        if not f in f2e:
            f2e[f] = [t]
        else:
            f2e[f].append(t)

    # compute neighbors by face
    for i in range(0, len(elem_face_nodes), nodes_per_face):
        f = tuple(sorted(elem_face_nodes[i:i+nodes_per_face]))
        t = elem_tags[i//(nodes_per_face*faces_per_elem)]
        if not t in e2e:
            e2e[t] = set()
        for tt in f2e[f]:
            if tt != t:
                e2e[t].add(tt)

    for k in e2e:
        for j in e2e[k]:
            G.add_edge(k,j)


    if nx.is_connected(G):
        print("connected")
        out_file.write(fname+'\n')
    else:
        print("disconnected")
    gmsh.finalize()

    print()
    sys.stdout.flush()
    out_file.flush()


    continue


    # nx.draw(G, with_labels = True)
    # plt.show()




    vertices_sorted = list(sorted(G))
    vertex_to_idx = {}

    for v_i in range(len(vertices_sorted)):
        vertex_to_idx[vertices_sorted[v_i]] = v_i

    edges_relabled = [(vertex_to_idx[edge[0]], vertex_to_idx[edge[1]]) for edge in G.edges()]

    # mlab.options.backend = 'envisage'

    pos = nx.spring_layout(G, dim=3, seed=1001)
    # numpy array of x,y,z positions in sorted node order
    xyz = np.array([pos[v] for v in vertices_sorted])
    # scalar colors
    scalars = np.array(list(G.nodes())) + 5

    mlab.figure()

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

    pts.mlab_source.dataset.lines = np.array(edges_relabled)
    tube = mlab.pipeline.tube(pts, tube_radius=0.02)
    mlab.pipeline.surface(tube, color=(0.8, 0.8, 0.8))
    mlab.draw()
    mlab.orientation_axes()
    mlab.show()



out_file.close()

