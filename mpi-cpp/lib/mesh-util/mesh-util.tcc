#include <string>
#include <vector>
#include "mpi.h"
#include "gmsh.h"
#include <cassert>

#include <unordered_map> 
#include <fstream>

#include "../util/util.hpp"

#include "../usort/ompUtils.h"
#include <stdexcept>

template <class T>
void GetElementsWithFacesNodesCentroids(const std::string &part_file_prefix, std::vector<T> &elements_out,
                          ElementType element_type, MPI_Comm comm)
{
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    std::string mesh_part_file_path = part_file_prefix + "_" + std::to_string(my_rank+1) + ".msh";

    // Initialize Gmsh
    gmsh::initialize();
    if (my_rank)
    {
        gmsh::option::setNumber("General.Terminal", 0);
    }
    // Load the mesh file
    gmsh::open(mesh_part_file_path);

    // Get the number of nodes and elements

    int gmsh_element_type;
    int gmsh_face_type;
    size_t nodes_per_face;
    size_t faces_per_element;
    size_t nodes_per_element;

    switch (element_type)
    {
    case ElementType::TET: // linear tet
    {
        gmsh_element_type=4;
        gmsh_face_type = 3; // triangle
        nodes_per_face = 3;
        faces_per_element = 4;
        nodes_per_element = 4;

        break;
    }
    case 5: // linear hexahedra
    {
        gmsh_element_type=5;
        gmsh_face_type = 4; // quadtriangle
        nodes_per_face = 4;
        faces_per_element = 6;
        nodes_per_element = 8;

        break;
    }

    default:
    {
        throw std::invalid_argument("unknown element type\t exiting...");
        break;
    }
    }

    std::vector<uint64_t> localElementNodeTags;
    std::vector<uint64_t> local_elem_tags;


    gmsh::model::mesh::getElementsByType(gmsh_element_type, local_elem_tags, localElementNodeTags, -1);

    // print_log("[", my_rank, "] local_elem_tags:", VectorToString(local_elem_tags));

  


    // uint64_t local_start = (my_rank * global_element_count) / procs_n;
    // uint64_t local_end = ((my_rank+1) * global_element_count) / procs_n;
    // uint64_t local_element_count = local_end-lolocal_startcal_start;

    size_t local_element_count = local_elem_tags.size();
    assert(localElementNodeTags.size() == local_element_count*nodes_per_element);


    std::vector<uint64_t> allNodeTags;
    std::vector<double> allNodeCoords, allNodeParams;
    
    gmsh::model::mesh::getNodes(allNodeTags, allNodeCoords, allNodeParams,-1,-1,true,false);

    std::unordered_map<uint64_t, uint64_t> node_tag_to_index; 
    for (size_t node_i = 0; node_i < allNodeTags.size(); node_i++)
    {
        node_tag_to_index[allNodeTags[node_i]] = node_i;
    }


    std::vector<double> local_elem_coordinates(local_element_count* 3); // 3 for 3D - x,y,z

    for (size_t element_i = 0; element_i < local_element_count; element_i++)
    {
        double x = 0;
        double y = 0;
        double z = 0;

        for (size_t elem_node_i = 0; elem_node_i < nodes_per_element; elem_node_i++)
        {
            size_t nodeTag = localElementNodeTags[element_i*nodes_per_element + elem_node_i];

            x += allNodeCoords[node_tag_to_index[nodeTag] * 3];
            y += allNodeCoords[node_tag_to_index[nodeTag] * 3 + 1];
            z += allNodeCoords[node_tag_to_index[nodeTag] * 3 + 2];
        }
        x = x / nodes_per_element;
        y = y / nodes_per_element;
        z = z / nodes_per_element;

        local_elem_coordinates[element_i * 3] = x;
        local_elem_coordinates[element_i * 3 + 1] = y;
        local_elem_coordinates[element_i * 3 + 2] = z;
    }

    // now we load unique face tag information from partitioned .bin files

    std::vector<uint64_t> elemTagsFromFile(local_element_count);
    std::vector<uint64_t> faceTagsFromFile(local_element_count*faces_per_element);
    std::unordered_map<uint64_t, uint64_t> fileElemToIdx;

    std::string elemtags_part_file_path = part_file_prefix + "_" + std::to_string(my_rank + 1) + "elemTags.bin";
    std::string facetags_part_file_path = part_file_prefix + "_" + std::to_string(my_rank + 1) + "faceTags.bin";

    std::ifstream infile{elemtags_part_file_path, std::ios::binary};
    infile.read(reinterpret_cast<char*>(elemTagsFromFile.data()), local_element_count * sizeof(uint64_t));
    infile.close();

    std::ifstream infile2{facetags_part_file_path, std::ios::binary};
    infile2.read(reinterpret_cast<char*>(faceTagsFromFile.data()), local_element_count * faces_per_element * sizeof(uint64_t));
    infile2.close();

    // print_log("[", my_rank, "] elemTagsFromFile:", VectorToString(elemTagsFromFile));
    // print_log("[", my_rank, "] faceTagsFromFile:", VectorToString(faceTagsFromFile));

    
    std::unordered_map<uint64_t, uint64_t> from_file_elem_tag_to_index; 
    for (size_t elem_i = 0; elem_i < local_element_count; elem_i++)
    {
        from_file_elem_tag_to_index[elemTagsFromFile[elem_i]] = elem_i;
    }


    // print_log("[", my_rank, "] localFaceTags:", VectorToString(localFaceTags));


    elements_out.resize(local_element_count);

    for (size_t element_i = 0; element_i < local_element_count; element_i++)
    {
        elements_out[element_i].element_tag = local_elem_tags[element_i];
        elements_out[element_i].x = local_elem_coordinates[element_i*3];
        elements_out[element_i].y = local_elem_coordinates[element_i*3+1];
        elements_out[element_i].z = local_elem_coordinates[element_i*3+2];

        for (size_t elem_node_i = 0; elem_node_i < nodes_per_element; elem_node_i++)
        {
            size_t nodeTag = localElementNodeTags[element_i*nodes_per_element + elem_node_i];
            elements_out[element_i].node_tags[elem_node_i] = nodeTag;
        }


        size_t element_i_in_file = from_file_elem_tag_to_index[local_elem_tags[element_i]];

        for (size_t elem_face_i = 0; elem_face_i < faces_per_element; elem_face_i++)
        {
            elements_out[element_i].face_tags[elem_face_i] = faceTagsFromFile[element_i_in_file*faces_per_element + elem_face_i];
        }
        


    }



    gmsh::finalize();

    // print_log("[", my_rank, "] elements:", VectorToString(elements_out));

    
}

template <class T>
void ResolveLocalElementConnectivity(const std::vector<T> &elements, ElementType element_type,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &connected_element_pairs_out,
                                std::vector<ElementWithFace> &unconnected_elements_faces_out)
{
    uint64_t faces_per_element;
    switch (element_type)
    {
    case ElementType::TET:
    {
        faces_per_element = 4;
        break;
    }
    case ElementType::HEX:
    {
        faces_per_element = 6;
        break;
    }
    default:
    {
        throw std::invalid_argument("unknown element type");
        break;
    }
    }
    std::vector<ElementWithFace> elements_with_faces;
    for (size_t element_i = 0; element_i < elements.size(); element_i++) {
        for (size_t elem_face_i = 0; elem_face_i < faces_per_element;
             elem_face_i++) {
            elements_with_faces.push_back(
                // ElementWithFace(elements[element_i].element_tag)
                {.element_tag = elements[element_i].element_tag,
                 .global_idx = elements[element_i].global_idx,
                 .face_tag = elements[element_i].face_tags[elem_face_i]});
        }
    }
        // std::vector<int> testvec = {456,34547,56,67,8967,956,85,6867,93,6,3452,3524};
        // omp_par::merge_sort(&testvec[0],&testvec[testvec.size()]);
        // print_log(VectorToString(testvec));


    omp_par::merge_sort(&elements_with_faces[0],&elements_with_faces[elements_with_faces.size()]);
    // print_log(VectorToString(elements_with_faces));

    connected_element_pairs_out.clear();
    unconnected_elements_faces_out.clear();

    {
        bool last_face_added = false;
        for (size_t elem_face_i = 1; elem_face_i < elements_with_faces.size(); elem_face_i++)
        {
            if (elements_with_faces[elem_face_i-1].face_tag == elements_with_faces[elem_face_i].face_tag)
            {
                if (last_face_added)
                {
                    throw std::runtime_error("more than two elements found for face : " + std::to_string(elements_with_faces[elem_face_i].face_tag));
                }else
                {
                    connected_element_pairs_out.push_back(
                        {
                            {
                                .element_tag = elements_with_faces[elem_face_i-1].element_tag,
                                .global_idx = elements_with_faces[elem_face_i-1].global_idx
                            },
                            {
                                .element_tag = elements_with_faces[elem_face_i].element_tag,
                                .global_idx = elements_with_faces[elem_face_i].global_idx

                            }   
                        }                     
                    );
                    last_face_added=true;
                }     
            }else
            {
                if (! last_face_added)
                {
                    unconnected_elements_faces_out.push_back(elements_with_faces[elem_face_i-1]);
                }
                
                last_face_added=false;
            }          
        }
        if (elements_with_faces[elements_with_faces.size()-2].face_tag != elements_with_faces[elements_with_faces.size()-1].face_tag)
        {
            unconnected_elements_faces_out.push_back(elements_with_faces[elements_with_faces.size()-1]);
        }
        
        
    }

    return;
    
    
    // const std::vector<T> elements_cpy(elements);
    // const std::vector<T> elements_sorted(elements);
}