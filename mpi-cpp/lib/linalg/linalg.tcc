#include <vector>
#include "../mesh-util/mesh-util.hpp"
#include <unordered_map>
#include <chrono>
#include <petscksp.h>
#include <petscerror.h>
#include <petscsys.h>
#include "petscviewer.h"
#include "util.hpp"



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
    // KSP             ksp; 

    const int spmv_repititions = 200;

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

    std::vector<PetscScalar> mat_row_values(elements.size()*(nodes_per_element*nodes_per_element));     
    
    {
        uint64_t idx = 0;
        for (auto & local_element : elements)
        {
            for (uint64_t node_x : local_element.node_tags)
            {
                auto node_x_global_idx = node_tag_to_global_idx_map[node_x];
                for (uint64_t node_y : local_element.node_tags)
                {
                    auto node_y_global_idx = node_tag_to_global_idx_map[node_y];
                    mat_row_values[idx++] = sin(node_x_global_idx+node_y_global_idx + 1);
                }
            }
        }
    }    
    PetscCallAbort(comm, MatSetValuesBatch(A,elements.size(),nodes_per_element,rows.data(),mat_row_values.data()));

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
    std::vector<PetscScalar> vec_row_values(elements.size()*nodes_per_element);
    {
        uint64_t idx = 0;
        for (auto & local_element : elements)
        {
            for (uint64_t node : local_element.node_tags)
            {

                PetscInt global_idx= static_cast<PetscInt>(node_tag_to_global_idx_map[node]);
                PetscScalar val = static_cast<PetscScalar>(cos(global_idx));
                vec_row_values[idx++] = val;         
            }     
        }
    }
    PetscCallAbort(comm, VecSetValues(x, elements.size()*nodes_per_element, rows.data(), vec_row_values.data(), ADD_VALUES));
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
    // MatView(A, PETSC_VIEWER_STDOUT_WORLD);

    // KSPCreate(PETSC_COMM_WORLD, &ksp);
    // KSPSetOperators(ksp, A, A);
    // KSPSetFromOptions(ksp);

    // warmup?
    MatMult(A, y, x);
    MatMult(A, x, y);

    // Perform the matrix-vector multiplication
    MPI_Barrier(comm);
    // PetscLogDefaultBegin();
    auto matvec_start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < spmv_repititions; i++)
    {
        if(i%2)
        {
            MatMult(A, y, x);
        }else
        {
            MatMult(A, x, y);
        }
        
    }
    MPI_Barrier(comm);
    auto matvec_end = std::chrono::high_resolution_clock::now();
    auto matvec_duration = std::chrono::duration_cast<std::chrono::microseconds>(matvec_end - matvec_start);
    if(!my_rank) print_log("matvec time: \t\t", matvec_duration.count(), "us");
    
    
    // PetscViewer viewer;
    // const char* file_name = "profiling.txt";
    // PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    // PetscViewerSetType(viewer, PETSCVIEWERASCII);
    // PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);
    // PetscViewerFileSetName(viewer, file_name);
    // PetscCallAbort(comm, PetscLogView(viewer));
    // PetscViewerDestroy(&viewer);


    // PetscLogView(PETSC_VIEWER_STDOUT_WORLD);

    // View the result
    // PetscPrintf(PETSC_COMM_WORLD, "Vec y:\n");
    // VecView(y, PETSC_VIEWER_STDOUT_WORLD);


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