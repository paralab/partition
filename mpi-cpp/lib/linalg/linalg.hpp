#include <vector>
#include "../mesh-util/mesh-util.hpp"


template <class T>
void TestSpMV(const std::vector<T> &elements, ElementType element_type, bool viz_flag, MPI_Comm comm);


#include "linalg.tcc"