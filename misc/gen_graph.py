#!/usr/bin/env python

import gmsh
import sys
import glob
#import graph_tool.all as gt
from networkx.drawing.nx_pydot import write_dot
import networkx as nx
import metis

if len(sys.argv) < 2:
    print("Usage: " + sys.argv[0] + " file")
    exit

gmsh.initialize(sys.argv)


for fname in glob.glob(r'Meshes/10k_hex/*.mesh'):
  # fname = sys.argv[1]
  # gmsh.open(sys.argv[1])
  gmsh.open(fname)
  print('Model ' + gmsh.model.getCurrent() + ' (' + str(gmsh.model.getDimension()) + 'D)')

  # gmsh.model.mesh.refine()

  # G = gt.Graph()
  G = nx.Graph()

  entities = gmsh.model.getEntities()
  mesh_type = ''

  # print('---'*50)
  for e in entities:
    dim = e[0]
    tag = e[1]

    if dim != 3:
      continue
    # print(dim, tag)
    type = gmsh.model.getType(e[0], e[1])
    name = gmsh.model.getEntityName(e[0], e[1])
    if len(name): name += ' '
    #print("Entity " + name + str(e) + " of type " + type)
    elemTypes, elemTags, elemNodeTags = gmsh.model.mesh.getElements(dim, tag)
    
    if elemTypes[0] == 4:
      mesh_type = 'tet'
      fn = 3
      en = 12 
      #print('Tetrahedral Mesh with ', len(elemTags[0]), ' elements')
      elems, _ = gmsh.model.mesh.getElementsByType(4)
      faces    = gmsh.model.mesh.getElementFaceNodes(4, 3)
    else:
      mesh_type = 'hex'
      fn = 4
      en = 24
      #print('Hexahedral Mesh', len(elemTags[0]), ' elements')
      elems, _ = gmsh.model.mesh.getElementsByType(5)
      faces    = gmsh.model.mesh.getElementFaceNodes(5, 4)
    
  # print('---'*50)
  print('  Have ', len(elems), ' elements and ', len(faces), ' faces.')

  # compute face x tet incidence
  f2e = {}
  e2e = {}

  V = {}
  for i,x in enumerate(elems):
    v = G.add_node(x)
    V[i] = x

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



  # print("neighbors by face: ", e2e)

  for k in e2e:
    for j in e2e[k]:
      G.add_edge(k,j)
      

  # gname = fname[:-5] + '.xml.gz'
  # G.save(gname)

  # if len(e2e) == len(elems):
  #   print(fname, ' : built e2e')
  # else:
  #   print(fname, ' : e2e ERROR') 

  # metis test 
  np = 13
  (edgecuts, parts) = metis.part_graph(G, np)


  # for i, p in enumerate(parts):
  #   G.nodes[V[i]]['part'] = p

  parts2 = {}

  for i,p in enumerate(parts):
    parts2[V[i]] = p

  # nx.write_dot(G, 'example.dot') # Requires pydot or pygraphviz

  ###------------ metrics ------------###
  #print(edgecuts)

  part_elems = [0] * np 
  part_sends = [0] * np
  # part_recvs = [0] * np


  for e in elems:
    part_elems[parts2[e]] += 1

  for k in e2e:
    for j in e2e[k]:
      if parts2[k] != parts2[j]:
        # send k - j
        part_sends[parts2[k]] += 1
        #part_recvs[parts2[j]] += 1
      

  print('  METIS Loads: ', part_elems)
  print('  METIS Sends: ', part_sends)
  #print(sum(part_sends))
  # print('Recvs: ', part_recvs)

###------------ Clean up ------------###
gmsh.clear()
gmsh.finalize()
