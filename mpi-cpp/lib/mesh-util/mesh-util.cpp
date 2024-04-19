#include "mesh-util.hpp"
#include "util.hpp"

#include "graph.hpp"
#include <gmsh.h>
#include <stdexcept>
#include <cassert>
#include <unordered_map>
#include <algorithm>

#define NO_ELEMENT SIZE_MAX


Graph GmshGetElementGraph(const std::string &mesh_file_path , std::vector<double>& elem_coordinates_out)
{
    // Initialize Gmsh
    gmsh::initialize();

    // Load the mesh file
    gmsh::open(mesh_file_path);

    // Get the number of nodes and elements
    std::vector<std::pair<int, int>> dimTags;
    gmsh::model::getEntities(dimTags);

    // int dim = dimTags[0].first;
    // int tag = dimTags[0].second;

    // if (dim != 3)
    // {
    //     throw std::invalid_argument("not a 3D mesh\t exiting...");
    // }

    std::vector<int> elementTypes;
    gmsh::model::mesh::getElementTypes(elementTypes, 3);

    if (elementTypes.size() == 0)
    {
        throw std::invalid_argument("no 3D elements\t exiting...");
    }

    if (elementTypes.size() > 1)
    {
        throw std::invalid_argument("more than 1 element type\t exiting...");
    }
    int element_type = elementTypes[0];
    int face_type;
    size_t nodes_per_face;
    size_t faces_per_element;
    size_t nodes_per_element;

    switch (element_type)
    {
    case 4: // linear tet
    {
        face_type = 3; // triangle
        nodes_per_face = 3;
        faces_per_element = 4;
        nodes_per_element=4;
        break;
    }
    case 5: // linear hexahedra
    {
        face_type = 4; // quadtriangle
        nodes_per_face = 4;
        faces_per_element = 6;
        nodes_per_element=8;
        std::cout << "linear hexahedra mesh\n";
        break;
    }

    default:
    {
        throw std::invalid_argument("unknown element type\t exiting...");
        break;
    }
    }

    std::vector<std::size_t> elementTags, elementNodeTags;
    gmsh::model::mesh::getElementsByType(element_type, elementTags, elementNodeTags);
    size_t element_count = elementTags.size();

    std::cout << "element count = " << element_count << std::endl;

    // std::cout << VectorToString(elementTags);



    std::vector<std::size_t> nodeTags;
    elem_coordinates_out.resize(element_count*3);       // 3 for 3D
    std::vector<double> nodeCoords(element_count*nodes_per_element*3);

    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodesByElementType(element_type,
                                          nodeTags,
                                          nodeCoords,
                                          parametricCoord,
                                          -1,
                                          false);
    for (size_t element_i = 0; element_i < element_count; element_i++)
    {
        double x = 0;
        double y=0;
        double z=0;

        for (size_t node_i = 0; node_i < nodes_per_element; node_i++)
        {
            x+= nodeCoords[element_i*nodes_per_element*3 + node_i*3];
            y+= nodeCoords[element_i*nodes_per_element*3 + node_i*3 + 1];
            z+= nodeCoords[element_i*nodes_per_element*3 + node_i*3 + 2];
        }
        x=x/nodes_per_element;
        y=y/nodes_per_element;
        z=z/nodes_per_element;

        elem_coordinates_out[element_i*3] = x;
        elem_coordinates_out[element_i*3+1] = y;
        elem_coordinates_out[element_i*3+2] = z;
        
        
    }
    



    std::vector<std::size_t> faceNodes;
    gmsh::model::mesh::getElementFaceNodes(element_type, face_type, faceNodes);
    assert(faceNodes.size() == (nodes_per_face*faces_per_element*element_count));
    // std::cout << VectorToString(faceNodes);

    std::unordered_map<std::string, std::pair<size_t, size_t>> face_to_element_pair;

    for (size_t node_i = 0; node_i < nodes_per_face*faces_per_element*element_count; node_i+=nodes_per_face)
    {
        size_t element = node_i / (nodes_per_face*faces_per_element);
        std::vector<size_t> nodes_in_face(nodes_per_face);
        for (size_t face_node_i = 0; face_node_i < nodes_per_face; face_node_i++)
        {
            nodes_in_face[face_node_i] = faceNodes[node_i + face_node_i];

        }
        
        std::sort(nodes_in_face.begin(),nodes_in_face.end());
        std::string face="";
        for (size_t face_node_i = 0; face_node_i < nodes_per_face; face_node_i++)
        {
            face+=std::to_string(nodes_in_face[face_node_i]) + " ";
        }
        
        // std::cout << "face hacsh = " << face_hash << std::endl;
        if (face_to_element_pair.find(face) == face_to_element_pair.end())
        {
            face_to_element_pair[face] = std::make_pair(element, NO_ELEMENT);
        }else
        {
            if (face_to_element_pair[face].second == NO_ELEMENT)
            {
                face_to_element_pair[face].second = element;
            }else
            {
                throw std::runtime_error("more than 2 elements found for same face");
            }           
        }      
        


    }
    Graph element_connectivity_graph;

    for (size_t elem_i = 0; elem_i < element_count; elem_i++)
    {
        element_connectivity_graph.AddVertex(elem_i);
    }
    
    for (auto face_mapping_i = face_to_element_pair.begin(); face_mapping_i != face_to_element_pair.end(); face_mapping_i++) {
        // std::cout <<"iterating face " << face_mapping_i->first << " : " << face_mapping_i->second.first <<", " << face_mapping_i->second.second <<"\n";
        if (face_mapping_i->second.second != NO_ELEMENT)
        {
            // std::cout <<"adding edge\n";
            element_connectivity_graph.AddEdge(face_mapping_i->second.first, face_mapping_i->second.second);
        }
        
    }

    // element_connectivity_graph.Print();
    
    

    // std::vector<size_t> elementTags, elementNodeTags;
    // int elementType = gmsh::model::mesh::getElementType("tetrahedron", 1);
    // gmsh::model::mesh::getElementsByType(elementType, elementTags, elementNodeTags);

    // Clean up Gmsh
    gmsh::finalize();

    return element_connectivity_graph;
}