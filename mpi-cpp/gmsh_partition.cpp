#include <iostream>
#include <gmsh.h>
#include <cassert>
#include <fstream>


int main(int argc, char const *argv[])
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << "<input mesh file path> <partition count> <output path>" << std::endl;
        return 1; // indicating an error
    }
    const std::string input_mesh_file_path = argv[1];
    int parts_n = std::stoi(argv[2]);
    const std::string output_path_prefix = argv[3];
    

    std::cout << "starting gmsh simple partitioning" << std::endl;

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0);
    
    gmsh::open(input_mesh_file_path);


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
    int gmsh_face_type;
    size_t nodes_per_face;
    size_t faces_per_element;
    size_t nodes_per_element;


    switch (gmsh_element_type)
    {
    case 4: // linear tet
    {
        gmsh_face_type = 3; // triangle
        nodes_per_face = 3;
        faces_per_element = 4;
        nodes_per_element = 4;
        // std::cout << "linear tetrahedra mesh\n";
        break;
    }
    case 5: // linear hexahedra
    {
        gmsh_face_type = 4; // quadtriangle
        nodes_per_face = 4;
        faces_per_element = 6;
        nodes_per_element = 8;
        // std::cout << "linear hexahedra mesh\n";
        break;
    }

    default:
    {
        throw std::invalid_argument("unknown element type\t exiting...");
        break;
    }
    }

    std::vector<size_t> elementNodeTags;
    std::vector<size_t> elem_tags;

    gmsh::model::mesh::getElementsByType(gmsh_element_type, elem_tags, elementNodeTags, -1);
    size_t element_count = elem_tags.size();

    elementNodeTags.clear();

    std::unordered_map<size_t, size_t> elem_tag_to_index; 
    for (size_t elem_i = 0; elem_i < element_count; elem_i++)
    {
        elem_tag_to_index[elem_tags[elem_i]] = elem_i;
    }

    // now we create all the face tags
    // this can not be done after partitioning because it results in non-unique face tags. 
    // gmsh is not designed to communicate among MPI processes to ensure unique face tags after partitioning.  

    std::vector<std::size_t> faceNodes;
    gmsh::model::mesh::getElementFaceNodes(gmsh_element_type, gmsh_face_type, faceNodes,-1,false);

    assert(faceNodes.size() == element_count*faces_per_element*nodes_per_face);

    gmsh::model::mesh::createFaces();

    std::vector<std::size_t> faceTags;
    std::vector<int> faceOrientations;
    gmsh::model::mesh::getFaces(gmsh_face_type, faceNodes, faceTags, faceOrientations);
    assert(faceTags.size() == (faces_per_element * element_count));

    // face tag creation done


    // now we do a simple checkerboard partitioning of the mesh.

    gmsh::option::setNumber("Mesh.PartitionSplitMeshFiles", 1.0);
    gmsh::option::setNumber("Mesh.PartitionCreateGhostCells", 0.0);
    gmsh::option::setNumber("Mesh.PartitionCreateTopology", 0.0);
    gmsh::option::setNumber("Mesh.PartitionCreatePhysicals", 0.0);


    // gmsh::option::setNumber("Mesh.PartitionOldStyleMsh2", 1.0);
    // gmsh::option::setNumber("Mesh.MshFileVersion", 2.2);

    gmsh::model::mesh::partition(parts_n);

    // gmsh::plugin::setNumber("SimplePartition", "NumSlicesX", parts_n);
    // gmsh::plugin::setNumber("SimplePartition", "NumSlicesY", 1);
    // gmsh::plugin::setNumber("SimplePartition", "NumSlicesZ", 1);
    // gmsh::plugin::run("SimplePartition");




    gmsh::write(output_path_prefix + ".msh");

    // partitioning done

    

    // now we save the face tag information in a partitioned way

    // taken from https://gitlab.onelab.info/gmsh/gmsh/-/issues/2333#note_19656
    std::vector<std::vector<uint64_t>> partition_to_elements(parts_n);
    std::vector<std::pair<int, int> > entities;
    gmsh::model::getEntities(entities,3);
    for (int part_i = 0; part_i < parts_n; part_i++)
    {      
        for(auto e : entities) {
            std::vector<int> entity_partitions;
            gmsh::model::getPartitions(e.first, e.second, entity_partitions);
            if(entity_partitions.size()) {
                assert(1 == entity_partitions.size());
                auto entity_part = entity_partitions[0];

                if ((entity_part - 1) == part_i)        // gmsh uses 1 based indexing for partition tags
                {
                    std::vector<size_t> elems_tmp;
                    std::vector<size_t> elem_nodes_tmp;

                    gmsh::model::mesh::getElementsByType(gmsh_element_type,elems_tmp,elem_nodes_tmp,e.second);
                    partition_to_elements[part_i].insert(partition_to_elements[part_i].end(), elems_tmp.begin(), elems_tmp.end());
                }                   
                                
            }
        }

    }

    for (size_t part_i = 0; part_i < parts_n; part_i++){
        // std::cout << "part: " << part_i << std::endl;
        size_t part_elem_count = partition_to_elements[part_i].size();
        std::vector<uint64_t> face_tags_in_part(part_elem_count * faces_per_element);
        size_t face_tag_idx = 0;
        for (auto part_elem: partition_to_elements[part_i])
        {
            // std::cout << "elem: " << part_elem << ": ";
            size_t elem_idx = elem_tag_to_index[part_elem];
            for (size_t face_i = 0; face_i < faces_per_element; face_i++)
            {
                // std::cout <<faceTags[faces_per_element*elem_idx + face_i] << " ";
                face_tags_in_part[face_tag_idx++] = static_cast<uint64_t>(faceTags[faces_per_element*elem_idx + face_i]);
            }
            // std::cout << std::endl;
            
        }

        std::ofstream outfile{output_path_prefix + "_" + std::to_string(part_i + 1) + "elemTags.bin", std::ios::binary};
        outfile.write(reinterpret_cast<const char*>(partition_to_elements[part_i].data()),
                      part_elem_count * sizeof(size_t));
        outfile.close();

        std::ofstream outfile2{output_path_prefix + "_" + std::to_string(part_i + 1) + "faceTags.bin", std::ios::binary};
        outfile2.write(reinterpret_cast<const char*>(face_tags_in_part.data()),
                      part_elem_count * faces_per_element * sizeof(size_t));
        outfile2.close();

        // above two files can be loaded later to re-construct the globally unique face tags
    }



    std::cout << "gmsh simple partitioning done" << std::endl;
    return 0;
}
