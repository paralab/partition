
import gmsh
import networkx as nx
import os



"""
returns
 - idx to element array
 - element to idx dict
 - e2e graph
 - element centroid coords array : 3xn list
 - valid: True/False
"""
def read_mesh_file(file_name: str):
    TET_ELEMENT_TYPE = 4
    HEX_ELEMENT_TYPE = 5

    TRIANGLE_FACE_TYEPE = 3
    QUAD_FACE_TYPE = 4
    gmsh.initialize()

    gmsh.open(file_name)



    print(': Model ' + gmsh.model.getCurrent() + ' (' + str(gmsh.model.getDimension()) + 'D)')

    G = nx.Graph()

    entities = gmsh.model.getEntities(3)
    mesh_type = ''


    if len(entities) == 0 :
        print(f"{file_name} : not a 3D mesh\n")
        gmsh.finalize()
        return {"valid": False}


    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(3)
    if len(elemTypes) !=1:
        print(f"{file_name} : more than 1 3D element type\n")
        gmsh.finalize()
        return {"valid": False}

    if elemTypes[0] == TET_ELEMENT_TYPE:
        mesh_type = 'tet'
        fn = 3
        en = 12
        vertices_n = 4 
        elems, _ = gmsh.model.mesh.getElementsByType(TET_ELEMENT_TYPE)
        faces    = gmsh.model.mesh.getElementFaceNodes(TET_ELEMENT_TYPE, TRIANGLE_FACE_TYEPE)
        _, nodeCoords, _ = gmsh.model.mesh.getNodesByElementType(4,returnParametricCoord=False)
    elif elemTypes[0] == HEX_ELEMENT_TYPE:
        mesh_type = 'hex'
        fn = 4
        en = 24
        vertices_n = 8
        elems, _ = gmsh.model.mesh.getElementsByType(HEX_ELEMENT_TYPE)
        faces    = gmsh.model.mesh.getElementFaceNodes(HEX_ELEMENT_TYPE, QUAD_FACE_TYPE)
        _, nodeCoords, _ = gmsh.model.mesh.getNodesByElementType(5,returnParametricCoord=False)
    else:
        print(f"{file_name} : element other than hex or tet\n")
        gmsh.finalize()
        return {"valid": False}


    print(mesh_type, 'mesh  has ', len(elems), ' elements and ', len(faces), ' faces.')

    f2e = {}
    e2e = {}

    idx_to_element = [None for i in range(len(elems))]
    element_to_idx = {}
    for i,x in enumerate(elems):
        G.add_node(x)
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
        print(f"{file_name} : not a connected mesh\n")
        gmsh.finalize()
        return {"valid": False}



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
        # setting geometric center as element coordinates
        elemCenterCoordsXYZ[elem_idx][0] = x_tot/vertices_n 
        elemCenterCoordsXYZ[elem_idx][1] = y_tot/vertices_n 
        elemCenterCoordsXYZ[elem_idx][2] = z_tot/vertices_n
        
    gmsh.finalize()
    return {"idx_to_element": idx_to_element, 
            "element_to_idx": element_to_idx, 
            "e2e_graph": G, 
            "element_coords": elemCenterCoordsXYZ, 
            "valid": True  }