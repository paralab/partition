#include "util.hpp"
#include "gmsh.h"

#include <stdexcept>
#include <cassert>

#include <Python.h>





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


// TODO: populate missing SFC_time
void ExportMetricsToPandasJson(
    std::string mesh_file, int file_idx, int run_idx, int partition_count, uint64_t global_vertex_count,
    std::vector<uint32_t>& sfc_partition_sizes, std::vector<uint32_t>& sfc_partition_boundaries,
    std::vector<uint32_t>& bfs_partition_sizes, std::vector<uint32_t>& bfs_partition_boundaries, int bfs_time,
    std::vector<uint32_t>& grow_partition_sizes, std::vector<uint32_t>& grow_partition_boundaries, int grow_time,
    std::vector<uint32_t>& parmetis_partition_sizes, std::vector<uint32_t>& parmetis_partition_boundaries, int parmetis_time,
    std::string metrics_out_file_path) {

    // const std::string pythonpath = "/home/budvin/research/Partitioning/paralab-partition/.venv";
    // setenv("PYTHONHOME", pythonpath.c_str(), 1);
    
    Py_Initialize();

    // const char *script_path = "/home/budvin/research/Partitioning/paralab-partition/mpi-cpp";

    // // Append the module directory to Python's sys.path
    // PyObject *sysPath = PySys_GetObject("path");
    // PyList_Append(sysPath, PyUnicode_FromString(script_path));
    
    PyObject* p_module = PyImport_ImportModule("metric_export");
    assert(p_module!=NULL);


    PyObject* p_func = PyObject_GetAttrString(p_module, "export_metrics");

    

    PyObject* py_mesh_file = PyUnicode_FromString(mesh_file.c_str());
    PyObject* py_file_idx = PyLong_FromLong(file_idx);
    PyObject* py_run_idx = PyLong_FromLong(run_idx);

    PyObject* py_partition_count = PyLong_FromLong(partition_count);
    PyObject* py_global_vertex_count = PyLong_FromUnsignedLong(global_vertex_count);
    

    PyObject* py_sfc_partition_sizes = PyList_New(sfc_partition_sizes.size());
    CopyToPythonList(sfc_partition_sizes, py_sfc_partition_sizes, sfc_partition_sizes.size());

    PyObject* py_sfc_partition_boundaries = PyList_New(sfc_partition_boundaries.size());
    CopyToPythonList(sfc_partition_boundaries, py_sfc_partition_boundaries, sfc_partition_boundaries.size());


    PyObject* py_bfs_partition_sizes = PyList_New(bfs_partition_sizes.size());
    CopyToPythonList(bfs_partition_sizes, py_bfs_partition_sizes, bfs_partition_sizes.size());

    PyObject* py_bfs_partition_boundaries = PyList_New(bfs_partition_boundaries.size());
    CopyToPythonList(bfs_partition_boundaries, py_bfs_partition_boundaries, bfs_partition_boundaries.size());

    PyObject* py_bfs_time = PyLong_FromLong(bfs_time);



    PyObject* py_grow_partition_sizes = PyList_New(grow_partition_sizes.size());
    CopyToPythonList(grow_partition_sizes, py_grow_partition_sizes, grow_partition_sizes.size());

    PyObject* py_grow_partition_boundaries = PyList_New(grow_partition_boundaries.size());
    CopyToPythonList(grow_partition_boundaries, py_grow_partition_boundaries, grow_partition_boundaries.size());

    PyObject* py_grow_time = PyLong_FromLong(grow_time);



    PyObject* py_parmetis_partition_sizes = PyList_New(parmetis_partition_sizes.size());
    CopyToPythonList(parmetis_partition_sizes, py_parmetis_partition_sizes, parmetis_partition_sizes.size());

    PyObject* py_parmetis_partition_boundaries = PyList_New(parmetis_partition_boundaries.size());
    CopyToPythonList(parmetis_partition_boundaries, py_parmetis_partition_boundaries,
                     parmetis_partition_boundaries.size());
    PyObject* py_parmetis_time = PyLong_FromLong(parmetis_time);


    PyObject* py_metrics_out_file_path = PyUnicode_FromString(metrics_out_file_path.c_str());

    PyObject* all_args =
        PyTuple_Pack(17, py_mesh_file, py_file_idx, py_run_idx, py_partition_count, py_global_vertex_count, py_sfc_partition_sizes,
                     py_sfc_partition_boundaries, py_bfs_partition_sizes, py_bfs_partition_boundaries, py_bfs_time,
                     py_grow_partition_sizes, py_grow_partition_boundaries, py_grow_time, py_parmetis_partition_sizes,
                     py_parmetis_partition_boundaries, py_parmetis_time, py_metrics_out_file_path);

    PyObject_CallObject(p_func, all_args);
    Py_Finalize();
}


