#include <vector>
#include "../mesh-util/mesh-util.hpp"


struct SpMVStatus
{
    int return_code;
    int mat_assembly_time_us;
    int matvec_time_us;

};

template <class T>
SpMVStatus TestSpMV(const std::vector<T> &elements, ElementType element_type, bool viz_flag, MPI_Comm comm);


#include "linalg.tcc"