#include <vector>
#include "../mesh-util/mesh-util.hpp"
#include <unordered_map>
#include <chrono>
#include <petscksp.h>
#include <petscerror.h>




template <class T>
SpMVStatus TestSpMV(const std::vector<T> &elements, ElementType element_type, bool viz_flag, MPI_Comm comm)
{

    int procs_n, my_rank;
    MPI_Comm_size(comm, &procs_n);
    MPI_Comm_rank(comm, &my_rank);
    if(!my_rank) print_log("starting SpMV");

    uint64_t    local_node_idx_start;
    uint64_t    local_node_idx_end;
    uint64_t    global_node_count;
    int         nodes_per_element;
    std::unordered_map<uint64_t, uint64_t> node_tag_to_global_idx_map;
    
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
        throw std::runtime_error("unknown element type");
        break;
    }
    }


    GetNodetagToGlobalIdx(elements,element_type,node_tag_to_global_idx_map,global_node_count, local_node_idx_start, local_node_idx_end,comm);


    uint64_t        local_count = local_node_idx_end - local_node_idx_start;
    Vec             x, y;       
    Mat             A;        

    // PetscFunctionBeginUser;
    PETSC_COMM_WORLD = comm;
    PetscCallAbort(comm, PetscInitializeNoArguments());


    // setting matrix
    PetscCallAbort(comm, MatCreate(PETSC_COMM_WORLD, &A));
    PetscCallAbort(comm, MatSetSizes(A, local_count, local_count, global_node_count, global_node_count));
    PetscCallAbort(comm, MatSetFromOptions(A));

    // populating the matrix with elemental matrices
    // for (auto & local_element : elements)
    // {
    //     for (uint64_t node_x : local_element.node_tags)
    //     {
    //         PetscInt global_x = static_cast<PetscInt>(node_tag_to_global_idx_map[node_x]);
    //         for (uint64_t node_y : local_element.node_tags)
    //         {
    //             PetscInt global_y = static_cast<PetscInt>(node_tag_to_global_idx_map[node_y]);
    //             PetscScalar val = 1.0;
    //             PetscCallAbort(comm, MatSetValue(A, global_x, global_y, val, ADD_VALUES));
    //         }
    //     }
    // }


    // populating the matrix using batch method
    std::vector<PetscInt> rows(elements.size()*nodes_per_element);
    
    {
        uint64_t idx = 0;
        for (auto & local_element : elements)
        {
            for (uint64_t node : local_element.node_tags)
            {
                PetscInt node_global_idx = static_cast<PetscInt>(node_tag_to_global_idx_map[node]);
                rows[idx++] = node_global_idx;
            }
        }
    }

    std::vector<PetscScalar> row_values(elements.size()*(nodes_per_element*nodes_per_element), 1.125);      // TODO: add some logic to populate values
    
    PetscCallAbort(comm, MatSetValuesBatch(A,elements.size(),nodes_per_element,rows.data(),row_values.data()));

    // matrix assembly
    MPI_Barrier(comm);
    auto mat_assemble_start = std::chrono::high_resolution_clock::now();
    PetscCallAbort(comm, MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY));
    PetscCallAbort(comm, MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY));
    MPI_Barrier(comm);
    auto mat_assemble_end = std::chrono::high_resolution_clock::now();
    auto mat_assemble_duration = std::chrono::duration_cast<std::chrono::microseconds>(mat_assemble_end - mat_assemble_start);

    if(!my_rank) print_log("matrix assembly time: \t", mat_assemble_duration.count(), "us");



    #ifdef ENABLE_PETSC_DRAW_FEATURES
    if (viz_flag)
    {
        PetscViewer     viewer;
        PetscCallAbort(comm, PetscViewerDrawOpen(PETSC_COMM_WORLD, NULL, "Matrix Structure", PETSC_DECIDE, PETSC_DECIDE, 500, 500, &viewer));
        PetscCallAbort(comm, PetscViewerDrawSetPause(viewer, -1));
        MatView(A, viewer);
    }
    #endif


    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // setup vectors for y = Ax
    PetscCallAbort(comm, MatCreateVecs(A, &x, &y));
    for (auto & local_element : elements)
    {
        for (uint64_t node : local_element.node_tags)
        {

            PetscInt global_idx= static_cast<PetscInt>(node_tag_to_global_idx_map[node]);
            PetscScalar val = 1.0;
            PetscCallAbort(comm, VecSetValue(x, global_idx, val, ADD_VALUES));            
        }     
    }

    // Assemble the vector x
    MPI_Barrier(comm);
    auto vecx_assemble_start = std::chrono::high_resolution_clock::now();
    VecAssemblyBegin(x);
    VecAssemblyEnd(x);
    MPI_Barrier(comm);
    auto vecx_assemble_end = std::chrono::high_resolution_clock::now();
    auto vecx_assemble_duration = std::chrono::duration_cast<std::chrono::microseconds>(vecx_assemble_end - vecx_assemble_start);
    if(!my_rank) print_log("vector assembly time: \t", vecx_assemble_duration.count(), "us");


    // PetscPrintf(PETSC_COMM_WORLD, "Vec x:\n");
    // VecView(x, PETSC_VIEWER_STDOUT_WORLD);

    // Perform the matrix-vector multiplication
    MPI_Barrier(comm);
    auto matvec_start = std::chrono::high_resolution_clock::now();
    MatMult(A, x, y);
    MPI_Barrier(comm);
    auto matvec_end = std::chrono::high_resolution_clock::now();
    auto matvec_duration = std::chrono::duration_cast<std::chrono::microseconds>(matvec_end - matvec_start);
    if(!my_rank) print_log("matvec time: \t\t", matvec_duration.count(), "us");

    // View the result
    // PetscPrintf(PETSC_COMM_WORLD, "Vec y:\n");
    // VecView(y, PETSC_VIEWER_STDOUT_WORLD);

    // Clean up
    VecDestroy(&x);
    VecDestroy(&y);
    MatDestroy(&A);

    PetscFinalize();
    if(!my_rank) print_log("SpMV done");

    SpMVStatus status;

    status.return_code = 0;
    status.mat_assembly_time_us = mat_assemble_duration.count();
    status.matvec_time_us = matvec_duration.count();

    return status;

}