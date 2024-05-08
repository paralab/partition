#include "mesh-util.hpp"
#include "util.hpp"

#include "../graph/graph.hpp"
#include "../usort/parUtils.h"
#include "../usort/ompUtils.h"
#include "../usort/dtypes.h"


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
        // std::cout << "linear tetrahedra mesh\n";
        break;
    }
    case 5: // linear hexahedra
    {
        element_type_out = ElementType::HEX;
        // std::cout << "linear hexahedra mesh\n";
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
    os << "(" << obj.element_tag << "," << obj.global_idx << ",[" <<  obj.x<< "," << obj.y<< "," << obj.z << "], " << obj.morton_encoding << ", [";
    for (size_t i = 0; i < 4; i++)
    {
        os << (i?",": " ") << obj.face_tags[i] ;
    
    }
    os << "])";
    
    return os;
}

// Overloading the << operator for HexElementWithFaces
std::ostream& operator<<(std::ostream& os, const HexElementWithFaces& obj) {
    os << "(" << obj.element_tag << "," << obj.global_idx << ",[" <<  obj.x<< "," << obj.y<< "," << obj.z << "], " << obj.morton_encoding << ", [";
    for (size_t i = 0; i < 6; i++)
    {
        os << (i?",": " ") << obj.face_tags[i] ;
    
    }
    os << "])";
    
    return os;
}

// Overloading the << operator for ElementWithFace
std::ostream& operator<<(std::ostream& os, const ElementWithFace& obj) {
    os << "(" << obj.element_tag << "," << obj.global_idx << "," << obj.face_tag << ")";
    
    return os;
}

// Overloading the << operator for ElementWithCoord
std::ostream& operator<<(std::ostream& os, const ElementWithCoord& obj) {
    os << "(" << obj.element_tag << "," << obj.global_idx  << ",[" << obj.x <<", "<< obj.y << ", " << obj.z << "])";
    
    return os;
}

// Overloading the << operator for ElementWithTag
std::ostream& operator<<(std::ostream& os, const ElementWithTag& obj) {
    os << "(" << obj.element_tag << "," << obj.global_idx  << ")";
    // os << obj.global_idx;
    
    return os;
}

void ResolveBoundaryElementConnectivity(std::vector<ElementWithFace> &unpaired_element_faces,
                                        std::vector<uint64_t> &proc_element_counts,
                                        std::vector<std::pair<ElementWithTag, ElementWithTag>> &boundary_connected_element_pairs_out,
                                        MPI_Comm comm)
{
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    uint64_t local_element_count = proc_element_counts[my_rank];
    std::vector<ElementWithFace> unpaired_element_faces_sorted(unpaired_element_faces.size());
    MPI_Barrier(MPI_COMM_WORLD);

    par::sampleSort<ElementWithFace>(unpaired_element_faces,unpaired_element_faces_sorted,comm);
    MPI_Barrier(MPI_COMM_WORLD);
    // print_log("[", my_rank, "]: unpaired_element_faces_sorted", VectorToString(unpaired_element_faces_sorted));

    std::vector<std::pair<ElementWithTag, ElementWithTag>> connected_boundary_element_pairs;

    {
        bool last_face_added = false;
        for (size_t elem_face_i = 1; elem_face_i < unpaired_element_faces_sorted.size(); elem_face_i++)
        {
            if (unpaired_element_faces_sorted[elem_face_i-1].face_tag == unpaired_element_faces_sorted[elem_face_i].face_tag)
            {
                if (last_face_added)
                {
                    throw std::runtime_error("more than two elements found for a face");
                }else
                {
                    connected_boundary_element_pairs.push_back(
                        {
                            {
                                .element_tag = unpaired_element_faces_sorted[elem_face_i-1].element_tag,
                                .global_idx = unpaired_element_faces_sorted[elem_face_i-1].global_idx
                            },
                            {
                                .element_tag = unpaired_element_faces_sorted[elem_face_i].element_tag,
                                .global_idx = unpaired_element_faces_sorted[elem_face_i].global_idx                            
                            }       
                        }                 
                    );
                    last_face_added=true;
                }     
            }else
            {
                
                last_face_added=false;
            }          
        }        
        
    }
    // print_log("[", my_rank, "]: connected_boundary_element_pairs", VectorToString(connected_boundary_element_pairs));


    /**
     * we need 2 copies of each pair, because we have to send boundary connectivity of a pair to 2 processors
     * it is guranteed that a the 2 elements in a pair belong to 2 different processors
     * if 2 elements a,b are connected the two pairs will be (a,b) and (b,a). 
    */
    std::vector<std::pair<ElementWithTag, ElementWithTag>> connected_boundary_element_pairs_cpy(connected_boundary_element_pairs);
    for (size_t pair_i = 0; pair_i < connected_boundary_element_pairs_cpy.size(); pair_i++)
    {
        std::swap(connected_boundary_element_pairs_cpy[pair_i].first, connected_boundary_element_pairs_cpy[pair_i].second);
    }

    /**
     * Then we sort the array with duplicate pairs according to the first element in the pair. 
     * Then we split it and send this to processors.
     * In this way if a,b are connected (a,b) will be send to owning process of 'a' and (b,a) will be sent to owning process of 'b'
     * 
    */

    std::vector<std::pair<ElementWithTag, ElementWithTag>>  connected_boundary_element_pairs_duplicated(connected_boundary_element_pairs);
    connected_boundary_element_pairs_duplicated.insert(
        connected_boundary_element_pairs_duplicated.end(), connected_boundary_element_pairs_cpy.begin(), connected_boundary_element_pairs_cpy.end());


    omp_par::merge_sort(&connected_boundary_element_pairs_duplicated[0],
                        &connected_boundary_element_pairs_duplicated[connected_boundary_element_pairs_duplicated.size()],
                        [](const auto& a, const auto& b) { 
                            return a.first.global_idx < b.first.global_idx; 
                                
                        });

    // print_log("[", my_rank, "]: connected_boundary_element_pairs_duplicated sorted", VectorToString(connected_boundary_element_pairs_duplicated));

    std::vector<uint64_t> proc_element_counts_scanned(procs_n);
    omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],procs_n);

    std::vector<int> send_counts(procs_n);
    {
        uint64_t current_proc = 0;
        for (size_t pair_i = 0; pair_i < connected_boundary_element_pairs_duplicated.size(); pair_i++)
        {
            if (connected_boundary_element_pairs_duplicated[pair_i].first.global_idx >= proc_element_counts_scanned[current_proc] &&
                connected_boundary_element_pairs_duplicated[pair_i].first.global_idx < (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc]))
            {
                send_counts[current_proc]++;
            } else
            {
                while (1)
                {
                    current_proc++;
                    if (connected_boundary_element_pairs_duplicated[pair_i].first.global_idx >= proc_element_counts_scanned[current_proc] &&
                        connected_boundary_element_pairs_duplicated[pair_i].first.global_idx < (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc]))
                    {
                        send_counts[current_proc]++;
                        break;
                    }                    
                }
                
            }   
        }       

    }
    // print_log("[", my_rank, "]: send_counts", VectorToString(send_counts));

    std::vector<int> recev_counts(procs_n);

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recev_counts.data(), 1, MPI_INT, comm);

    // print_log("[", my_rank, "]: recev_counts", VectorToString(recev_counts));

    std::vector<int> send_displs(procs_n);
    omp_par::scan(&send_counts[0],&send_displs[0],procs_n);
    // print_log("[", my_rank, "]: send_displs", VectorToString(send_displs));



    std::vector<int> recv_displs(procs_n);
    omp_par::scan(&recev_counts[0],&recv_displs[0],procs_n);
    // print_log("[", my_rank, "]: recv_displs", VectorToString(recv_displs));

    int total_receive_count = recev_counts[procs_n-1]+recv_displs[procs_n-1];
    
    boundary_connected_element_pairs_out.resize(total_receive_count);
    MPI_Alltoallv(connected_boundary_element_pairs_duplicated.data(),
                  send_counts.data(), send_displs.data(), par::Mpi_pairtype<ElementWithTag,ElementWithTag>::value(), 
                  boundary_connected_element_pairs_out.data(), recev_counts.data(), recv_displs.data(),
                  par::Mpi_pairtype<ElementWithTag,ElementWithTag>::value(), comm);

    // print_log("[", my_rank, "]: received bdry connectivities", VectorToString(boundary_connected_element_pairs_out));

}

/**
 * Get ghost elements from the boundary edges
 * ghost_element_counts_out is populated with ghost elements, sorted by global_idx
 * ghost_element_counts_out is populated with ghost element count for each process
*/
void ExtractGhostElements(std::vector<std::pair<ElementWithTag, ElementWithTag>>& boundary_connected_element_pairs,
                          std::vector<uint64_t>& proc_element_counts,
                          std::vector<uint64_t>& proc_element_counts_scanned,
                          std::vector<ElementWithTag>& ghost_elements_out, std::vector<int>& ghost_element_counts_out,
                          MPI_Comm comm) {
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    ghost_elements_out.clear();
    ghost_element_counts_out.resize(procs_n);
    std::fill(ghost_element_counts_out.begin(), ghost_element_counts_out.end(), 0);

    if (boundary_connected_element_pairs.size() == 0)
    {
        return;
    }
    

    std::vector<ElementWithTag> ghost_elements_dups(boundary_connected_element_pairs.size());

    // TODO: can be parallelized
    for (size_t boudary_connect_i = 0; boudary_connect_i < boundary_connected_element_pairs.size();
         boudary_connect_i++) {
        ghost_elements_dups[boudary_connect_i].element_tag =
            boundary_connected_element_pairs[boudary_connect_i].second.element_tag;
        ghost_elements_dups[boudary_connect_i].global_idx =
            boundary_connected_element_pairs[boudary_connect_i].second.global_idx;
    }

    omp_par::merge_sort(&ghost_elements_dups[0], &ghost_elements_dups[ghost_elements_dups.size()]);

    // print_log("[", taskid, "]:", "ghost_elements_dups sorted = ", VectorToString(ghost_elements_dups));



    /**
     * removing duplicates
    */
    ghost_elements_out.push_back(ghost_elements_dups[0]);
    for (size_t ghost_element_dup_i = 1; ghost_element_dup_i < ghost_elements_dups.size(); ghost_element_dup_i++) {
        if (ghost_elements_dups[ghost_element_dup_i].global_idx !=
            ghost_elements_dups[ghost_element_dup_i - 1].global_idx) {
            ghost_elements_out.push_back(ghost_elements_dups[ghost_element_dup_i]);
        }
    }

    // print_log("[", my_rank, "]:", "ghost_elements = ", VectorToString(ghost_elements_out));

    {
        uint64_t current_proc = 0;
        for (size_t ghost_element_i = 0; ghost_element_i < ghost_elements_out.size(); ghost_element_i++) {
            if (ghost_elements_out[ghost_element_i].global_idx >= proc_element_counts_scanned[current_proc] &&
                ghost_elements_out[ghost_element_i].global_idx <
                    (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc])) {
                ghost_element_counts_out[current_proc]++;
            } else {
                while (1) {
                    current_proc++;
                    if (ghost_elements_out[ghost_element_i].global_idx >= proc_element_counts_scanned[current_proc] &&
                        ghost_elements_out[ghost_element_i].global_idx <
                            (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc])) {
                        ghost_element_counts_out[current_proc]++;
                        break;
                    }
                }
            }
        }
    }
    // print_log("[", my_rank, "]:", "ghost_element_counts = ", VectorToString(ghost_element_counts_out));
}