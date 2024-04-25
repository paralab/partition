#include "mesh-util.hpp"
#include "util.hpp"

#include "graph.hpp"
#include <gmsh.h>
#include <stdexcept>
#include <cassert>
#include <unordered_map>
#include <algorithm>

#include <mpi.h>

#define NO_ELEMENT SIZE_MAX


Graph GmshGetElementGraph(const std::string &mesh_file_path , std::vector<double>& elem_coordinates_out, std::vector<size_t>& elem_tags)
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


    std::vector<std::size_t> elementNodeTags;
    gmsh::model::mesh::getElementsByType(element_type, elem_tags, elementNodeTags);
    size_t element_count = elem_tags.size();

    std::cout << "element count = " << element_count << std::endl;

    // std::cout << VectorToString(elem_tags);



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
    gmsh::model::mesh::createFaces();
    std::vector<std::size_t> faceTags;
    std::vector<int> faceOrientations;
    gmsh::model::mesh::getFaces(face_type, faceNodes, faceTags, faceOrientations);
    assert(faceTags.size() == (faces_per_element*element_count));


    // std::unordered_map<std::string, std::pair<size_t, size_t>> face_to_element_pair;
    std::vector<std::pair<size_t,size_t>> element_face_pairs;           // (elemTag, faceTag) pairs

    for (size_t face_i = 0; face_i < faces_per_element*element_count; face_i+=faces_per_element)
    {
        size_t element = elem_tags[face_i / (faces_per_element)];
        // std::vector<size_t> nodes_in_face(nodes_per_face);
        for (size_t elem_face_i = 0; elem_face_i < faces_per_element; elem_face_i++)
        {
            element_face_pairs.push_back({element,faceTags[face_i+elem_face_i]});
        } 
        


    }
    stable_sort(element_face_pairs.begin(), element_face_pairs.end(),
                [&](std::pair<size_t,size_t> pair_1, std::pair<size_t,size_t> pair_2){ 
                    return pair_1.second < pair_2.second;
                });
    Graph element_connectivity_graph;

    for (size_t elem_i = 0; elem_i < element_count; elem_i++)
    {
        element_connectivity_graph.AddVertex(elem_tags[elem_i]);
    }
    bool last_face_added = false;
    for (size_t pair_i = 1; pair_i < faces_per_element*element_count; pair_i++)
    {
        if (element_face_pairs[pair_i-1].second == element_face_pairs[pair_i].second)
        {
            if (last_face_added)
            {
                throw std::runtime_error("more than 2 elements found for same face");
                
            }
            element_connectivity_graph.AddEdge(element_face_pairs[pair_i-1].first, element_face_pairs[pair_i].first);
            last_face_added=true;            
        }else
        {
            last_face_added=false;
        }
        

        
    }
    
    

    // element_connectivity_graph.Print();
    
    

    // std::vector<size_t> elem_tags, elementNodeTags;
    // int elementType = gmsh::model::mesh::getElementType("tetrahedron", 1);
    // gmsh::model::mesh::getElementsByType(elementType, elem_tags, elementNodeTags);

    // Clean up Gmsh
    gmsh::finalize();

    return element_connectivity_graph;
}

ElementType GetElementType(const std::string &mesh_file_path, MPI_Comm comm){
    int  my_rank;
    MPI_Comm_rank(comm, &my_rank);
    // Initialize Gmsh
    gmsh::initialize();
    if (my_rank)
    {
        gmsh::option::setNumber("General.Verbosity", 0);
    }
    

    // Load the mesh file
    gmsh::open(mesh_file_path);

    // Get the number of nodes and elements
    std::vector<std::pair<int, int>> dimTags;
    gmsh::model::getEntities(dimTags);



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
    int gmsh_element_type = elementTypes[0];
    ElementType element_type_out;

    switch (gmsh_element_type)
    {
    case 4: // linear tet
    {
        element_type_out = ElementType::TET;
        std::cout << "linear tetrahedra mesh\n";
        break;
    }
    case 5: // linear hexahedra
    {
        element_type_out = ElementType::HEX;
        std::cout << "linear hexahedra mesh\n";
        break;
    }

    default:
    {
        throw std::invalid_argument("unknown element type\t exiting...");
        break;
    }
    }

    gmsh::finalize();

    return element_type_out;



}

// Overloading the << operator for TetElementWithFaces
std::ostream& operator<<(std::ostream& os, const TetElementWithFaces& obj) {
    os << "(" << obj.element_tag << ",[" <<  obj.x<< "," << obj.y<< "," << obj.z << "], " << obj.morton_encoding << ", [";
    for (size_t i = 0; i < 4; i++)
    {
        os << (i?",": " ") << obj.face_tags[i] ;
    
    }
    os << "])";
    
    return os;
}

// Overloading the << operator for HexElementWithFaces
std::ostream& operator<<(std::ostream& os, const HexElementWithFaces& obj) {
    os << "(" << obj.element_tag << ",[" <<  obj.x<< "," << obj.y<< "," << obj.z << "], " << obj.morton_encoding << ", [";
    for (size_t i = 0; i < 6; i++)
    {
        os << (i?",": " ") << obj.face_tags[i] ;
    
    }
    os << "])";
    
    return os;
}