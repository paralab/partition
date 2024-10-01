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
#include <chrono>

template <class T>
void GetInitialElementsDistribution(const std::string &mesh_file_path, std::vector<T> &local_elements_out,
                          ElementType element_type, MPI_Comm comm)
{
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    int local_element_count;
    std::vector<int> proc_element_counts(procs_n);          // populated only at root
    std::vector<int> proc_element_counts_scanned(procs_n);  // populated only at root
    std::vector<T> all_elements;                            // populated only at root

    // root process will read the mesh and do a simple initial element scatter distribution
    if(!my_rank)
    {

        gmsh::initialize();
        // gmsh::option::setNumber("General.Terminal", 0);
        gmsh::open(mesh_file_path);

        int gmsh_element_type;
        int gmsh_face_type;
        int faces_per_element;
        int nodes_per_element;
        int nodes_per_face;

        switch (element_type)
        {
        case ElementType::TET: // linear tet
        {
            gmsh_element_type   = 4;
            gmsh_face_type      = 3; // triangle
            nodes_per_face      = 3;
            faces_per_element   = 4;
            nodes_per_element   = 4;

            break;
        }
        case 5: // linear hexahedra
        {
            gmsh_element_type   = 5;
            gmsh_face_type      = 4; // quadtriangle
            nodes_per_face      = 4;
            faces_per_element   = 6;
            nodes_per_element   = 8;

            break;
        }

        default:
        {
            throw std::invalid_argument("unknown element type\t exiting...");
            break;
        }
        }

        std::vector<uint64_t> element_node_tags;
        std::vector<uint64_t> element_tags;


        gmsh::model::mesh::getElementsByType(gmsh_element_type, element_tags, element_node_tags, -1);

        // print_log("[", my_rank, "] element_tags:", VectorToString(element_tags));


        size_t total_element_count = element_tags.size();
        assert(element_node_tags.size() == total_element_count*nodes_per_element);


        std::vector<uint64_t> allNodeTags;
        std::vector<double> allNodeCoords, allNodeParams;
        
        gmsh::model::mesh::getNodes(allNodeTags, allNodeCoords, allNodeParams,-1,-1,true,false);

        std::unordered_map<uint64_t, uint64_t> node_tag_to_index; 
        for (size_t node_i = 0; node_i < allNodeTags.size(); node_i++)
        {
            node_tag_to_index[allNodeTags[node_i]] = node_i;
        }


        std::vector<double> element_coordinates(total_element_count* 3); // 3 for 3D - x,y,z

        for (size_t element_i = 0; element_i < total_element_count; element_i++)
        {
            double x = 0;
            double y = 0;
            double z = 0;

            for (size_t elem_node_i = 0; elem_node_i < nodes_per_element; elem_node_i++)
            {
                size_t nodeTag = element_node_tags[element_i*nodes_per_element + elem_node_i];

                x += allNodeCoords[node_tag_to_index[nodeTag] * 3];
                y += allNodeCoords[node_tag_to_index[nodeTag] * 3 + 1];
                z += allNodeCoords[node_tag_to_index[nodeTag] * 3 + 2];
            }
            x = x / nodes_per_element;
            y = y / nodes_per_element;
            z = z / nodes_per_element;

            element_coordinates[element_i * 3       ] = x;
            element_coordinates[element_i * 3 + 1   ] = y;
            element_coordinates[element_i * 3 + 2   ] = z;
        }

        std::vector<std::size_t> faceNodes;
        gmsh::model::mesh::getElementFaceNodes(gmsh_element_type, gmsh_face_type, faceNodes,-1,false);

        assert(faceNodes.size() == total_element_count*faces_per_element*nodes_per_face);

        gmsh::model::mesh::createFaces();

        std::vector<std::size_t> faceTags;
        std::vector<int> faceOrientations;
        gmsh::model::mesh::getFaces(gmsh_face_type, faceNodes, faceTags, faceOrientations);
        assert(faceTags.size() == (faces_per_element * total_element_count));
        gmsh::finalize();

        all_elements.resize(total_element_count);

        for (size_t element_i = 0; element_i < total_element_count; element_i++)
        {
            all_elements[element_i].element_tag = element_tags[element_i];
            all_elements[element_i].x = element_coordinates[element_i*3];
            all_elements[element_i].y = element_coordinates[element_i*3+1];
            all_elements[element_i].z = element_coordinates[element_i*3+2];

            for (size_t elem_node_i = 0; elem_node_i < nodes_per_element; elem_node_i++)
            {
                size_t nodeTag = element_node_tags[element_i*nodes_per_element + elem_node_i];
                all_elements[element_i].node_tags[elem_node_i] = nodeTag;
            }


            for (size_t elem_face_i = 0; elem_face_i < faces_per_element; elem_face_i++)
            {
                all_elements[element_i].face_tags[elem_face_i] = faceTags[element_i*faces_per_element + elem_face_i];
            }    
        }
        std::fill(proc_element_counts.begin(), proc_element_counts.end(), total_element_count/procs_n);

        // Distribute the remainder element counts
        int rem = total_element_count % procs_n;
        for (int proc_i = 0; proc_i < rem; proc_i++) {
            proc_element_counts[proc_i]++;
        }
        omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],procs_n);
        print_log("mesh reading done");
    }   // mesh reading done

    MPI_Scatter(proc_element_counts.data(), 1, MPI_INT, &local_element_count, 1, MPI_INT, 0, comm);
    local_elements_out.resize(local_element_count);
    MPI_Scatterv(all_elements.data(),
        proc_element_counts.data(), proc_element_counts_scanned.data(), par::Mpi_datatype<T>::value(),
        local_elements_out.data(),local_element_count,par::Mpi_datatype<T>::value(), 
        0, comm);

    // print_log("[", my_rank, "] elements:", VectorToString(local_elements_out));

    
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
                {/* .element_tag = */ elements[element_i].element_tag,
                 /* .global_idx = */ elements[element_i].global_idx,
                 /* .face_tag = */ elements[element_i].face_tags[elem_face_i]});
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
                                /* .element_tag = */ elements_with_faces[elem_face_i-1].element_tag,
                                /* .global_idx = */ elements_with_faces[elem_face_i-1].global_idx
                            },
                            {
                                /* .element_tag = */ elements_with_faces[elem_face_i].element_tag,
                                /* .global_idx = */ elements_with_faces[elem_face_i].global_idx

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
 * Creates element connectivity pairs based on node sharing logic
 */
template <class T>
void ResolveElementConnectivityByNodes(const std::vector<T> &elements, ElementType element_type,
                                std::vector<uint64_t>& proc_element_counts,
                                std::vector<uint64_t>& proc_element_counts_scanned,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &local_connected_element_pairs_out,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &boundary_connected_element_pairs_out,
                                MPI_Comm comm)
{
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);

    uint64_t nodes_per_element;
    switch (element_type)
    {
    case ElementType::TET:
    {
        nodes_per_element = 4;
        break;
    }
    case ElementType::HEX:
    {
        nodes_per_element = 8;
        break;
    }
    default:
    {
        throw std::invalid_argument("unknown element type");
        break;
    }
    }
    std::vector<ElementWithNode> elements_with_nodes;
    for (size_t element_i = 0; element_i < elements.size(); element_i++) {
        for (size_t elem_node_i = 0; elem_node_i < nodes_per_element;
             elem_node_i++) {
            elements_with_nodes.push_back(
                {/* .element_tag = */ elements[element_i].element_tag,
                 /* .global_idx = */ elements[element_i].global_idx,
                 /* .node_tag = */ elements[element_i].node_tags[elem_node_i]});
        }
    }


    std::vector<ElementWithNode> elements_with_nodes_sorted;
    MPI_Barrier(comm);

    par::sampleSort<ElementWithNode>(elements_with_nodes,elements_with_nodes_sorted,comm);
    MPI_Barrier(comm);

    std::vector<std::pair<ElementWithTag, ElementWithTag>> connected_element_pairs;


    {
        size_t current_node;
        size_t current_node_start_i = 0;
        size_t current_node_end_i   = 0;
        size_t elem_node_i = 0;

        while (elem_node_i < elements_with_nodes_sorted.size())
        {
            current_node = elements_with_nodes_sorted[elem_node_i].node_tag;
            current_node_start_i = elem_node_i;
            elem_node_i++;
            

            while (elem_node_i < elements_with_nodes_sorted.size() &&
                    (elements_with_nodes_sorted[elem_node_i].node_tag == elements_with_nodes_sorted[elem_node_i-1].node_tag))
            {
                elem_node_i++;
            }

            current_node_end_i = elem_node_i;

            // nested loop in all elements sharing a node. this nested loop will create pairs for connectivity.
            // this will generate both pairs (a,b) and (b,a) if elements 'a' and 'b' are connected
            for (size_t elem_node_j = current_node_start_i; elem_node_j < current_node_end_i; elem_node_j++)
            {
                for (size_t elem_node_k = current_node_start_i; elem_node_k < current_node_end_i; elem_node_k++)
                {
                    if (elem_node_j == elem_node_k)
                    {
                        continue;
                    }
                    
                    connected_element_pairs.push_back(
                        {
                            {
                                .element_tag = elements_with_nodes_sorted[elem_node_j].element_tag,
                                .global_idx = elements_with_nodes_sorted[elem_node_j].global_idx
                            },
                            {
                                .element_tag = elements_with_nodes_sorted[elem_node_k].element_tag,
                                .global_idx = elements_with_nodes_sorted[elem_node_k].global_idx                            
                            }       
                        }                 
                    );
                }
            }              
        }
        
    }

    // now sort the pairs by the first element
    omp_par::merge_sort(&connected_element_pairs[0],
                        &connected_element_pairs[connected_element_pairs.size()],
                        [](const auto& a, const auto& b) { 
                            return a.first.global_idx < b.first.global_idx; 
                                
                        });

    // now send pairs to respective processes.
    // if elements 'a' anb 'b' are connected, pair (a,b) will be sent to the owning process of 'a' and the pair (b,a) will be sent to the owning process of 'b'
    // note: since the connectivity is defined by node sharing logic, potentially there can be multiples of the same pair
    //       in a subsequent step, we get remove these duplicates

    std::vector<int> send_counts(procs_n);
    {
        uint64_t current_proc = 0;
        for (size_t pair_i = 0; pair_i < connected_element_pairs.size(); pair_i++)
        {
            if (connected_element_pairs[pair_i].first.global_idx >= proc_element_counts_scanned[current_proc] &&
                connected_element_pairs[pair_i].first.global_idx < (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc]))
            {
                send_counts[current_proc]++;
            } else
            {
                while (1)
                {
                    current_proc++;
                    if (connected_element_pairs[pair_i].first.global_idx >= proc_element_counts_scanned[current_proc] &&
                        connected_element_pairs[pair_i].first.global_idx < (proc_element_counts_scanned[current_proc] + proc_element_counts[current_proc]))
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


    // this will contain element connectvity for both local-local and local-ghost
    std::vector<std::pair<ElementWithTag, ElementWithTag>> connected_element_pairs_local_and_ghost(total_receive_count);

    
    MPI_Alltoallv(connected_element_pairs.data(),
                  send_counts.data(), send_displs.data(), par::Mpi_pairtype<ElementWithTag,ElementWithTag>::value(), 
                  connected_element_pairs_local_and_ghost.data(), recev_counts.data(), recv_displs.data(),
                  par::Mpi_pairtype<ElementWithTag,ElementWithTag>::value(), comm);

    // now we have to remove duplicates in pairs.
    // if a,b are connected there can be multiple copies of (a,b) and (b,a) and maybe multiples of both.
    // however if a pair (a,b) corresponds to a ghost edge, it is guranteed that we do not have (b,a) in connected_element_pairs_local_and_ghost

    // preprocess step: if both a,b in a pair (a,b) or (b,a) are owned locally, rearrange the pair s.t. first.global_idx < second.global_idx
    
    auto local_start = proc_element_counts_scanned[my_rank];
    auto local_end = proc_element_counts_scanned[my_rank] + proc_element_counts[my_rank];

    for (size_t pair_i = 0; pair_i < connected_element_pairs_local_and_ghost.size(); pair_i++)
    {
        if (    connected_element_pairs_local_and_ghost[pair_i].first.global_idx >= local_start &&
                connected_element_pairs_local_and_ghost[pair_i].first.global_idx < local_end    &&
                connected_element_pairs_local_and_ghost[pair_i].second.global_idx >= local_start &&
                connected_element_pairs_local_and_ghost[pair_i].second.global_idx < local_end)
        {
            if (connected_element_pairs_local_and_ghost[pair_i].second.global_idx < connected_element_pairs_local_and_ghost[pair_i].first.global_idx ) 
            {
                std::swap(connected_element_pairs_local_and_ghost[pair_i].first, connected_element_pairs_local_and_ghost[pair_i].second);
            }
        }         
    }
    

    // normal sort by second element
    std::sort(  connected_element_pairs_local_and_ghost.begin(), 
                connected_element_pairs_local_and_ghost.end(),
                [](const auto& a, const auto& b) { 
                    return a.second.global_idx < b.second.global_idx; 
                                
                });
    

    // stable sort by first element
    std::stable_sort(   connected_element_pairs_local_and_ghost.begin(), 
                        connected_element_pairs_local_and_ghost.end(),
                        [](const auto& a, const auto& b) { 
                            return a.first.global_idx < b.first.global_idx; 
                                
                        });

    // now duplicates are in contigous range


    local_connected_element_pairs_out.clear();
    boundary_connected_element_pairs_out.clear();

    if (connected_element_pairs_local_and_ghost.size() == 0)
    {
        return;
    }

    if (connected_element_pairs_local_and_ghost[0].second.global_idx >= local_start &&
        connected_element_pairs_local_and_ghost[0].second.global_idx < local_end)
    {
        local_connected_element_pairs_out.push_back(connected_element_pairs_local_and_ghost[0]);
    }else
    {
        boundary_connected_element_pairs_out.push_back(connected_element_pairs_local_and_ghost[0]);
    }

    for (size_t pair_i = 1; pair_i < connected_element_pairs_local_and_ghost.size(); pair_i++)
    {

        if (connected_element_pairs_local_and_ghost[pair_i].first.global_idx == connected_element_pairs_local_and_ghost[pair_i-1].first.global_idx &&
            connected_element_pairs_local_and_ghost[pair_i].second.global_idx == connected_element_pairs_local_and_ghost[pair_i-1].second.global_idx)
        {
            // duplicate detected, skipping
            continue;
        }
        if (connected_element_pairs_local_and_ghost[pair_i].second.global_idx >= local_start &&
            connected_element_pairs_local_and_ghost[pair_i].second.global_idx < local_end)
        {
            local_connected_element_pairs_out.push_back(connected_element_pairs_local_and_ghost[pair_i]);
        }else
        {
            boundary_connected_element_pairs_out.push_back(connected_element_pairs_local_and_ghost[pair_i]);
        }
        
    }   
    // if (!my_rank)
    // {
    //     print_log("[", my_rank, "]: local_connected_element_pairs", VectorToString(local_connected_element_pairs_out));
    //     print_log("[", my_rank, "]: boundary_connected_element_pairs", VectorToString(boundary_connected_element_pairs_out));
    // }
    



    
}

/**
 * Performs sample sort on TET or HEX elements on morton_encoding field, using lightweight proxy objects for distributed sorting step.
 */
template <class T>
void SampleSortMorton(std::vector<T> &elements_in, std::vector<T> &elements_out, MPI_Comm comm){
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    if(!my_rank) print_log("using new sample sort");

    /**
     * gather all pre-sort local count information
     */

    int pre_sort_local_count = elements_in.size();

    std::vector<int> pre_sort_all_local_counts(procs_n);
    std::vector<int> pre_sort_all_local_counts_scanned(procs_n);

    MPI_Allgather(&pre_sort_local_count, 1, MPI_INT,pre_sort_all_local_counts.data(),1,MPI_INT,comm);

    omp_par::scan(&pre_sort_all_local_counts[0],&pre_sort_all_local_counts_scanned[0],procs_n);

    /**
     * form a new aux array with (global_idx, morton_encoding information)
     */
    std::vector<SortingElement> aux_array_to_sort(pre_sort_local_count);
    for (int i = 0; i < pre_sort_local_count; i++)
    {
        aux_array_to_sort[i].morton_encoding = elements_in[i].morton_encoding;
        aux_array_to_sort[i].global_idx = pre_sort_all_local_counts_scanned[my_rank] + i;        // global contiguous indexing before sorting
    }
    
    /**
     * global sample sort the aux array by morton_encoding
    */
    std::vector<SortingElement> aux_array_sorted;
    MPI_Barrier(comm);
    par::sampleSort<SortingElement>(aux_array_to_sort,aux_array_sorted,comm);


    int sorted_local_count = aux_array_sorted.size();

    std::vector<int> sorted_local_counts(procs_n);
    std::vector<int> sorted_local_counts_scanned(procs_n);

    MPI_Allgather(&sorted_local_count, 1, MPI_INT,sorted_local_counts.data(),1,MPI_INT,comm);

    omp_par::scan(&sorted_local_counts[0],&sorted_local_counts_scanned[0],procs_n);

    /**
     * temp array to store pairs<original global index, sorted local index>
     */
    std::vector<std::pair<uint64_t, uint64_t>> indices(sorted_local_count);
    for (int i = 0; i < sorted_local_count; i++)
    {
        indices[i].first = aux_array_sorted[i].global_idx;
        indices[i].second = i;
    }

    /**
     * sort by original global index
     */
    std::sort(indices.begin(), indices.end(), 
        [](auto &left, auto &right) {
            return left.first < right.first;
        }
    );

    /**
     * now send the global indices to owner processes to request the actual data objects
     */
    
    std::vector<int> idx_send_counts(procs_n,0);
    {
        uint64_t current_proc = 0;
        for (uint64_t i = 0; i < indices.size(); i++)
        {
            if (indices[i].first < static_cast<uint64_t>(pre_sort_all_local_counts[current_proc] + pre_sort_all_local_counts_scanned[current_proc]))
            {
                idx_send_counts[current_proc]++;
            } else
            {
                while (1)
                {
                    current_proc++;
                    if (indices[i].first < static_cast<uint64_t>(pre_sort_all_local_counts[current_proc] + pre_sort_all_local_counts_scanned[current_proc]))
                    {
                        idx_send_counts[current_proc]++;
                        break;
                    }                    
                }
                
            }   
        }       

    }

    std::vector<int> idx_send_displs(procs_n);
    omp_par::scan(&idx_send_counts[0],&idx_send_displs[0],procs_n);


    std::vector<int> idx_recev_counts(procs_n);
    MPI_Alltoall(idx_send_counts.data(), 1, MPI_INT, idx_recev_counts.data(), 1, MPI_INT, comm);

    std::vector<int> idx_recv_displs(procs_n);
    omp_par::scan(&idx_recev_counts[0],&idx_recv_displs[0],procs_n);

    /**
     * prepare buffer to send indices
     */
    std::vector<uint64_t> sending_indices(sorted_local_count);
    for (int i = 0; i < sorted_local_count; i++)
    {
        sending_indices[i] = indices[i].first;

    }

    int idx_recev_total_count = idx_recev_counts[procs_n-1]+idx_recv_displs[procs_n-1];


    assert(idx_recev_total_count == pre_sort_local_count);

    std::vector<uint64_t> recev_indices(idx_recev_total_count);
    

    MPI_Alltoallv(sending_indices.data(),
                  idx_send_counts.data(), idx_send_displs.data(), MPI_UINT64_T, 
                  recev_indices.data(), idx_recev_counts.data(), idx_recv_displs.data(),
                  MPI_UINT64_T, comm);

    /**
     * now each process have the requested global indices.
     * rearrange the elements according to this index order and prepare for sending
     */

    std::vector<T> sending_elements(pre_sort_local_count);

    for (int i = 0; i < pre_sort_local_count; i++)
    {
        //recev_indices contains the pre-sort global indices. to get the local pre sort local index, we have to offset the scanned pre sort count
        sending_elements[i] = elements_in[recev_indices[i] - pre_sort_all_local_counts_scanned[my_rank]];
    }

    /**
     * now when determining element send recev counts, it is essentially the flipped version of idx send recv counts
     */

    std::vector<int> elem_send_counts = idx_recev_counts;
    std::vector<int> elem_send_displas = idx_recv_displs;

    std::vector<int> elem_recv_counts = idx_send_counts;
    std::vector<int> elem_recv_displs = idx_send_displs;

    std::vector<T> recv_elements(sorted_local_count);


    MPI_Alltoallv(sending_elements.data(),
                  elem_send_counts.data(), elem_send_displas.data(), par::Mpi_datatype<T>::value(), 
                  recv_elements.data(), elem_recv_counts.data(), elem_recv_displs.data(),
                  par::Mpi_datatype<T>::value(), comm);

    /**
     * now we have received the elements, rearrange them to the morton sorted order
     */

    elements_out.resize(sorted_local_count);

    for (int i = 0; i < sorted_local_count; i++)
    {
        elements_out[indices[i].second] = recv_elements[i];
    }

    
    
}


/**
 * Given a new labeling (i.e. a new partitioning) this function redistributes elements_in.
 * Output is in elements_out
 * Final local elements are further sorted according to morton encoding
 * Also adjusts global_idx of elements to new partitioned range
 * 
*/
template <class T>
DistributionStatus Redistribute(std::vector<T> &elements_in, std::vector<uint16_t>& labeling, std::vector<T> &elements_out, MPI_Comm comm){
    
    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    if(!my_rank) print_log("starting redistribution");
    MPI_Barrier(comm);
    auto start = std::chrono::high_resolution_clock::now();

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

    std::vector<T> elements_ordered(elements_in.size());

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

    MPI_Allgather(&new_element_count, 1, MPI_UINT64_T,proc_element_counts.data(),1,MPI_UINT64_T,comm);


    omp_par::scan(&proc_element_counts[0],&proc_element_counts_scanned[0],procs_n);

    uint64_t global_idx_start = proc_element_counts_scanned[my_rank];

    for (size_t local_elem_i = 0; local_elem_i < new_element_count; local_elem_i++)
    {
        elements_out[local_elem_i].global_idx = global_idx_start + local_elem_i;
    }

    // print_log("[", my_rank, "]: new elements", VectorToString(elements_out));
    MPI_Barrier(comm);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    if(!my_rank)
    {
        print_log("redistribution done");
        print_log("element redistribution time:\t", duration.count(), " us");
    } 
    DistributionStatus status;
    status.return_code = 0;
    status.time_us = duration.count();

    return status;


    
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