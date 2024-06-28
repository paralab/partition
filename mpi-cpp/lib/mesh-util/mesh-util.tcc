#include <string>
#include <vector>
#include "mpi.h"
#include "gmsh.h"
#include <cassert>

#include <unordered_map> 
#include <unordered_set> 

#include <fstream>

#include "../util/util.hpp"

#include "../usort/ompUtils.h"
#include "../usort/parUtils.h"

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


/**
 * Given a new labeling (i.e. a new partitioning) this function redestributes elements_in.
 * Output is in elements_out
 * Final local elements are further sorted according to morton encoding
 * Also adjusts global_idx of elements to new partitioned range
 * 
*/
template <class T>
void Redestribute(std::vector<T> &elements_in, std::vector<uint16_t>& labeling, std::vector<T> &elements_out, MPI_Comm comm){
    
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    if(!my_rank) print_log("starting redistribution");

    assert(elements_in.size() == labeling.size());
    
    std::vector<std::pair<uint64_t, uint16_t>> idx_label_pairs(elements_in.size());    // temp array for sorting

    for (uint64_t i = 0; i < elements_in.size(); i++)
    {        
        idx_label_pairs[i] = {i, labeling[i]};
    }

    std::sort(idx_label_pairs.begin(), idx_label_pairs.end(), 
        [](auto &left, auto &right) {
            return left.second < right.second;
        }
    );

    const std::vector<T> elements_ordered(elements_in.size());

    for (uint64_t i = 0; i < elements_in.size(); i++)
    {        
        elements_ordered[i] = elements_in[idx_label_pairs[i].first];
    }




    std::vector<int> send_counts(procs_n);
    {
        uint64_t current_proc = 0;
        for (uint64_t i = 0; i < elements_ordered.size(); i++)
        {
            if (idx_label_pairs[i].second == current_proc)
            {
                send_counts[current_proc]++;
            } else
            {
                while (1)
                {
                    current_proc++;
                    if (idx_label_pairs[i].second == current_proc)
                    {
                        send_counts[current_proc]++;
                        break;
                    }                    
                }
                
            }   
        }       

    }
    // print_log("[", my_rank, "]: elements_in", VectorToString(elements_in));
    // print_log("[", my_rank, "]: labeling", VectorToString(labeling));


    std::vector<int> recev_counts(procs_n);

    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recev_counts.data(), 1, MPI_INT, comm);

    // print_log("[", my_rank, "]: recev_counts", VectorToString(recev_counts));

    std::vector<int> send_displs(procs_n);
    omp_par::scan(&send_counts[0],&send_displs[0],procs_n);
    // print_log("[", my_rank, "]: send_displs", VectorToString(send_displs));



    std::vector<int> recv_displs(procs_n);
    omp_par::scan(&recev_counts[0],&recv_displs[0],procs_n);
    // print_log("[", my_rank, "]: recv_displs", VectorToString(recv_displs));

    uint64_t new_element_count = recev_counts[procs_n-1]+recv_displs[procs_n-1];
    
    elements_out.resize(new_element_count);
    MPI_Alltoallv(elements_ordered.data(),
                  send_counts.data(), send_displs.data(), par::Mpi_datatype<T>::value(), 
                  elements_out.data(), recev_counts.data(), recv_displs.data(),
                  par::Mpi_datatype<T>::value(), comm);

    // local sorting w.r.t. morton encoding
    omp_par::merge_sort(&elements_out[0], &elements_out[elements_out.size()]);

    // adjuting global_idx
    std::vector<uint64_t> proc_element_counts(procs_n);
    std::vector<uint64_t> proc_element_counts_scanned(procs_n);

    MPI_Allgather(&new_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,MPI_COMM_WORLD);


    omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],procs_n);

    uint64_t global_idx_start = proc_element_counts_scanned[my_rank];

    for (size_t local_elem_i = 0; local_elem_i < new_element_count; local_elem_i++)
    {
        elements_out[local_elem_i].global_idx = global_idx_start + local_elem_i;
    }

    // print_log("[", my_rank, "]: new elements", VectorToString(elements_out));

    if(!my_rank) print_log("redistribution done");


    
}


/**
 * Given a vector of tet or hex elements (with node tags), computes a global numbering for node tags to be used with petsc matvec routines.
 * node tags are renumbered in a way that each mpi rank gets a contiguous section in the numbering.
 * When a node tag is present in multiple MPI processes, its owenership goes to the lowest ranked MPI process
 * 
 * @param[in] local_elements vector of `TetElementWithFacesNodes` or `HexElementWithFacesNodes`
 * @param[in] element_type either `ElementType::TET` or `ElementType::HEX`
 * @param[out] mapping_out output mapping from local node tags to global unique index
 * @param[out] global_count_out number of globally unique nodes
 * @param[out] local_start_out start of the owned range in new numbering index
 * @param[out] local_start_end end of the owned range in new numbering index
 * 
 *  
*/
template <class T>
void GetNodetagToGlobalIdx(const std::vector<T> &local_elements, ElementType element_type, 
        std::unordered_map<uint64_t, uint64_t>& mapping_out, uint64_t& global_count_out, uint64_t& local_start_out, uint64_t& local_end_out, MPI_Comm comm)
{
    /**
     * steps followed
     * 
     * 1. get all local unique nodes to a set
     * 2. make pairs (node_tag, my_rank)
     * 3. global sample sort by `node_tag`
     * 4. for each pair (node_tag, my_rank) make struct (node_tag=node_tag, location=my_rank, owner, global_idx)
     * 5. since structs are already sorted by node_tag populate `owner` in structs by selecting the lowest rank for each node_tag
     * 6. sort by owner
     * 7. reditribute to owners
     * 8. now owner has nodes, do a local sort by node tag, then assign global_idx,
     * 9. local sort by location and reditribute to location
     * 10. using the received node new information, populate unordered map (node, global_idx)
    */

    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);


    /**
     * step 1. get all local unique nodes to a set
    */
    std::unordered_set<uint64_t> all_local_nodes_set;

    for (auto & element : local_elements)
    {
        for (uint64_t node : element.node_tags)
        {
            all_local_nodes_set.insert(node);            
        }
    }

    std::vector<uint64_t> all_local_nodes(all_local_nodes_set.begin(), all_local_nodes_set.end());

    // omp_par::merge_sort(&all_local_nodes[0],&all_local_nodes[all_local_nodes.size()]);
    // print_log("[", my_rank, "]: all_local_nodes", VectorToString(all_local_nodes));

    /**
     * step 2. make pairs (node_tag, my_rank)
    */
    std::vector<NodeLocationPair> node_location_pairs(all_local_nodes.size());

    for (size_t node_i = 0; node_i < all_local_nodes.size(); node_i++)
    {
        node_location_pairs[node_i].node_tag = all_local_nodes[node_i];
        node_location_pairs[node_i].location_rank = my_rank;
    }

    /**
     * step 3. global sample sort by `node_tag`
    */
    std::vector<NodeLocationPair> node_location_pairs_sorted;
    MPI_Barrier(comm);
    par::sampleSort<NodeLocationPair>(node_location_pairs, node_location_pairs_sorted, comm); 
    MPI_Barrier(comm);

    /**
     * step 4. for each pair (node_tag, my_rank) make struct (node_tag=node_tag, location=my_rank, owner, global_idx)
    */
    std::vector<NodeNewInfo> node_new_info(node_location_pairs_sorted.size());
    {
        size_t pair_i = 0;
        while (pair_i < node_location_pairs_sorted.size())
        {
            size_t start_ = pair_i;
            size_t end_ = pair_i;
            uint64_t current_node = node_location_pairs_sorted[pair_i].node_tag;
            int current_owned_rank = node_location_pairs_sorted[pair_i].location_rank;

            pair_i++;
            end_++;

            while (pair_i < node_location_pairs_sorted.size() && (node_location_pairs_sorted[pair_i].node_tag == current_node))
            {
                current_owned_rank = std::min(current_owned_rank, node_location_pairs_sorted[pair_i].location_rank);
                pair_i++;
                end_++;
            }
            /**
             *  step 5. since structs are already sorted by node_tag populate `owner` in structs by selecting the lowest rank for each node_tag
             *  TODO: if required, change the logic for node ownership. currently ownership goes to the lowest rank MPI process. 
            */
            for (size_t temp_i = start_; temp_i < end_; temp_i++)
            {
                node_new_info[temp_i].node_tag = node_location_pairs_sorted[temp_i].node_tag;
                node_new_info[temp_i].location_rank = node_location_pairs_sorted[temp_i].location_rank;
                node_new_info[temp_i].owner_rank = current_owned_rank;
                // global_idx will be populated later
            }
        }
    }

    /**
     * step 6. sort by owner
    */
    std::sort(node_new_info.begin(), node_new_info.end(), 
        [](auto &left, auto &right) {
            return left.owner_rank < right.owner_rank;
        }
    );    


    // for mpi comm
    std::vector<int> send_counts(procs_n, 0);
    std::vector<int> recev_counts(procs_n, 0);
    std::vector<int> send_displs(procs_n, 0);
    std::vector<int> recv_displs(procs_n, 0);
    int total_receive_count;



    /**
     * step 7. reditribute to owners
    */
    {
        int current_proc = 0;
        for (size_t node_info_i = 0; node_info_i < node_new_info.size(); node_info_i++)
        {
            if (node_new_info[node_info_i].owner_rank == current_proc)
            {
                send_counts[current_proc]++;
            } else
            {
                while (1)
                {
                    current_proc++;
                    if (node_new_info[node_info_i].owner_rank == current_proc)
                    {
                        send_counts[current_proc]++;
                        break;
                    }                    
                }
                
            }   
        }       
    }



    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recev_counts.data(), 1, MPI_INT, comm);

    omp_par::scan(&send_counts[0],&send_displs[0],procs_n);

    omp_par::scan(&recev_counts[0],&recv_displs[0],procs_n);

    total_receive_count = recev_counts[procs_n-1]+recv_displs[procs_n-1];

    std::vector<NodeNewInfo> nodes_for_glob_idx_assign(total_receive_count);
   

    MPI_Alltoallv(node_new_info.data(),
                  send_counts.data(), send_displs.data(), par::Mpi_datatype<NodeNewInfo>::value(), 
                  nodes_for_glob_idx_assign.data(), recev_counts.data(), recv_displs.data(),
                  par::Mpi_datatype<NodeNewInfo>::value(), comm);

    /**
     * step 8. now owner has nodes, do a local sort by node tag, then assign global_idx,
    */
    omp_par::merge_sort(&nodes_for_glob_idx_assign[0], &nodes_for_glob_idx_assign[nodes_for_glob_idx_assign.size()]);

    //counting unique node tags
    uint64_t local_uniq_node_count = 0;

    if (nodes_for_glob_idx_assign.size() > 0)
    {
        local_uniq_node_count = 1;

        for (size_t node_i = 1; node_i < nodes_for_glob_idx_assign.size(); node_i++)
        {
            if (nodes_for_glob_idx_assign[node_i].node_tag != nodes_for_glob_idx_assign[node_i-1].node_tag)
            {
                local_uniq_node_count++;
            }
            
        }
    }

    std::vector<uint64_t> proc_uniq_node_counts(procs_n, 0);
    std::vector<uint64_t> proc_uniq_node_counts_scanned(procs_n, 0);
    uint64_t global_uniq_node_count  = 0;
    

    MPI_Allgather(&local_uniq_node_count, 1, MPI_UINT64_T,proc_uniq_node_counts.data(),1,MPI_UINT64_T,comm);

    MPI_Allreduce(&local_uniq_node_count,&global_uniq_node_count,1,MPI_UINT64_T,MPI_SUM,comm);

    omp_par::scan(&proc_uniq_node_counts[0],&proc_uniq_node_counts_scanned[0],procs_n);

    //now we can calculate the contigous range for global_idx to assign to nodes
    uint64_t node_global_idx_start = proc_uniq_node_counts_scanned[my_rank];
    uint64_t node_global_idx_end =  node_global_idx_start + proc_uniq_node_counts[my_rank];


    if (nodes_for_glob_idx_assign.size() > 0){
        uint64_t current_global_idx = node_global_idx_start;
        nodes_for_glob_idx_assign[0].global_idx = current_global_idx;

        for (size_t node_i = 1; node_i < nodes_for_glob_idx_assign.size(); node_i++)
        {
            if (nodes_for_glob_idx_assign[node_i].node_tag != nodes_for_glob_idx_assign[node_i-1].node_tag)
            {
                current_global_idx++;
            }
            nodes_for_glob_idx_assign[node_i].global_idx = current_global_idx;
            
        }
        assert(current_global_idx == (node_global_idx_end - 1));

    }   // global_idx assignment done.

    // print_log("[", my_rank, "]: own nodes", VectorToString(nodes_for_glob_idx_assign));

   

    /**
     * 9. local sort by location and redistribute to location
    */
    std::sort(nodes_for_glob_idx_assign.begin(), nodes_for_glob_idx_assign.end(), 
        [](auto &left, auto &right) {
            return left.location_rank < right.location_rank;
        }
    );  

    std::fill(send_counts.begin(), send_counts.end(), 0);
    std::fill(recev_counts.begin(), recev_counts.end(), 0);
    std::fill(send_displs.begin(), send_displs.end(), 0);
    std::fill(recv_displs.begin(), recv_displs.end(), 0);


    total_receive_count = 0;

    {
        int current_proc = 0;
        for (size_t node_info_i = 0; node_info_i < nodes_for_glob_idx_assign.size(); node_info_i++)
        {
            if (nodes_for_glob_idx_assign[node_info_i].location_rank == current_proc)
            {
                send_counts[current_proc]++;
            } else
            {
                while (1)
                {
                    current_proc++;
                    if (nodes_for_glob_idx_assign[node_info_i].location_rank == current_proc)
                    {
                        send_counts[current_proc]++;
                        break;
                    }                    
                }
                
            }   
        }       
    }



    MPI_Alltoall(send_counts.data(), 1, MPI_INT, recev_counts.data(), 1, MPI_INT, comm);

    omp_par::scan(&send_counts[0],&send_displs[0],procs_n);

    omp_par::scan(&recev_counts[0],&recv_displs[0],procs_n);

    total_receive_count = recev_counts[procs_n-1]+recv_displs[procs_n-1];

    std::vector<NodeNewInfo> nodes_with_glob_idx(total_receive_count);
   

    MPI_Alltoallv(nodes_for_glob_idx_assign.data(),
                  send_counts.data(), send_displs.data(), par::Mpi_datatype<NodeNewInfo>::value(), 
                  nodes_with_glob_idx.data(), recev_counts.data(), recv_displs.data(),
                  par::Mpi_datatype<NodeNewInfo>::value(), comm);

    assert(nodes_with_glob_idx.size() == all_local_nodes.size());

    /**
     * step 10. using the received node new information, populate unordered map (node, global_idx)
    */

    mapping_out.clear();


    for (auto & node : nodes_with_glob_idx)
    {
        mapping_out[node.node_tag] = node.global_idx;        
    }


    global_count_out = global_uniq_node_count;
    local_start_out = node_global_idx_start;
    local_end_out = node_global_idx_end;

    return;

    

}