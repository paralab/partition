#include "util.hpp"
#include "gmsh.h"

#include <stdexcept>
#include <cassert>

#include <Python.h>

#include <vector>





template <> 
std::string VectorToString(std::vector<uint8_t> vec){
    std::ostringstream output;
    output << "[ ";
    for (auto element : vec)
    {
        output << +element << ", ";
    }
    output << "]\n";
    return output.str();
}

void AssignPartitionLabelsInOrder(std::vector<uint64_t> &ordering, uint64_t count, uint64_t partition_count, std::vector<uint64_t> &labels_out){
    assert(labels_out.size() == count);
    uint64_t partition_size = count/partition_count;
    uint64_t large_partition_count = count%partition_count;

    for (uint64_t p_i = 0; p_i < partition_count; p_i++)
    {
        uint64_t size;
        uint64_t offset;
        if (p_i < large_partition_count)
        {
            size = partition_size+1;
            offset = p_i*(partition_size+1);
        }else
        {
            size = partition_size;
            offset = large_partition_count*(partition_size+1) + (p_i-large_partition_count)*partition_size;
        }

        for (uint64_t i = offset; i < (offset+size); i++)
        {
            labels_out[ordering[i]] = p_i;
        }       
        
        
    }
    return;

}

void GetSamplesFromOrdered(std::vector<uint64_t> &order, std::vector<uint64_t> &input_arr, uint64_t sample_count, std::vector<uint64_t> &samples_out){
    assert(input_arr.size() == order.size());
    assert(samples_out.size() == sample_count);
    for (uint64_t sample_i = 0; sample_i < sample_count; sample_i++)
    {
        samples_out[sample_i] = input_arr[order[sample_i*(input_arr.size()/sample_count)]];
    }
    
}

// TODO: template if needed
void CopyToPythonList(std::vector<uint32_t>& in_vector, PyObject* out_list, size_t count) {
    for (size_t i = 0; i < count; i++) {
        PyList_SetItem(out_list, i, PyLong_FromUnsignedLong(in_vector[i]));
    }
}


void ExportMetricsToPandasJson(
    std::string mesh_file, int file_idx, int run_idx, int partition_count, uint64_t global_vertex_count,
    std::vector<uint32_t>& sfc_partition_sizes, std::vector<uint32_t>& sfc_partition_boundaries, int sfc_partition_time, int sfc_mat_assembly_time, int sfc_matvec_time,
    std::vector<uint32_t>& bfs_partition_sizes, std::vector<uint32_t>& bfs_partition_boundaries, int bfs_labeling_time, int bfs_redistribution_time, int bfs_mat_assembly_time, int bfs_matvec_time,
    std::vector<uint32_t>& parmetis_partition_sizes, std::vector<uint32_t>& parmetis_partition_boundaries, int parmetis_labeling_time, int parmetis_redistribution_time, int parmetis_mat_assembly_time, int parmetis_matvec_time,
    std::string metrics_out_file_path) {


    
    Py_Initialize();

    
    PyObject* p_module = PyImport_ImportModule("metric_export");
    assert(p_module!=NULL);


    PyObject* p_func = PyObject_GetAttrString(p_module, "export_metrics");

    

    PyObject* py_mesh_file = PyUnicode_FromString(mesh_file.c_str());
    PyObject* py_file_idx = PyLong_FromLong(file_idx);
    PyObject* py_run_idx = PyLong_FromLong(run_idx);

    PyObject* py_partition_count = PyLong_FromLong(partition_count);
    PyObject* py_global_vertex_count = PyLong_FromUnsignedLong(global_vertex_count);



    PyObject* py_sfc_partition_sizes;
    PyObject* py_sfc_partition_boundaries;

    PyObject* py_bfs_partition_sizes;
    PyObject* py_bfs_partition_boundaries;

    PyObject* py_parmetis_partition_sizes;
    PyObject* py_parmetis_partition_boundaries;

    std::vector<PyObject**> py_lists = {&py_sfc_partition_sizes, &py_sfc_partition_boundaries, &py_bfs_partition_sizes, &py_bfs_partition_boundaries, &py_parmetis_partition_sizes, &py_parmetis_partition_boundaries};
    std::vector<std::vector<uint32_t>*> cpp_vectors = {&sfc_partition_sizes, &sfc_partition_boundaries, &bfs_partition_sizes, &bfs_partition_boundaries, &parmetis_partition_sizes, &parmetis_partition_boundaries};



    for (int list_i = 0; list_i < py_lists.size(); list_i++)
    {
        // print_log(py_lists[list_i]);
        // print_log(VectorToString(*cpp_vectors[list_i]));
        *py_lists[list_i] = PyList_New(cpp_vectors[list_i]->size());
        CopyToPythonList(*cpp_vectors[list_i], *py_lists[list_i], cpp_vectors[list_i]->size());
        // print_log(py_lists[list_i]);


    }


    PyObject* py_sfc_partition_time = PyLong_FromLong(sfc_partition_time);
    PyObject* py_sfc_mat_assembly_time = PyLong_FromLong(sfc_mat_assembly_time);
    PyObject* py_sfc_matvec_time = PyLong_FromLong(sfc_matvec_time);

    PyObject* py_bfs_labeling_time = PyLong_FromLong(bfs_labeling_time);
    PyObject* py_bfs_redistribution_time = PyLong_FromLong(bfs_redistribution_time);
    PyObject* py_bfs_mat_assembly_time = PyLong_FromLong(bfs_mat_assembly_time);
    PyObject* py_bfs_matvec_time = PyLong_FromLong(bfs_matvec_time);


    PyObject* py_parmetis_labeling_time = PyLong_FromLong(parmetis_labeling_time);
    PyObject* py_parmetis_redistribution_time = PyLong_FromLong(parmetis_redistribution_time);
    PyObject* py_parmetis_mat_assembly_time = PyLong_FromLong(parmetis_mat_assembly_time);
    PyObject* py_parmetis_matvec_time = PyLong_FromLong(parmetis_matvec_time);


    PyObject* py_metrics_out_file_path = PyUnicode_FromString(metrics_out_file_path.c_str());
    PyObject* all_args =
        PyTuple_Pack(23, py_mesh_file, py_file_idx, py_run_idx, py_partition_count, py_global_vertex_count, 
                        py_sfc_partition_sizes, py_sfc_partition_boundaries, py_sfc_partition_time, py_sfc_mat_assembly_time, py_sfc_matvec_time,
                        py_bfs_partition_sizes, py_bfs_partition_boundaries, py_bfs_labeling_time, py_bfs_redistribution_time, py_bfs_mat_assembly_time, py_bfs_matvec_time,
                        py_parmetis_partition_sizes, py_parmetis_partition_boundaries, py_parmetis_labeling_time, py_parmetis_redistribution_time, py_parmetis_mat_assembly_time, py_parmetis_matvec_time,
                        py_metrics_out_file_path);

    PyObject_CallObject(p_func, all_args);

    PyErr_Print();
    Py_Finalize();
}


