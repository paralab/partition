#ifndef SFC_H
#define SFC_H

#include <mesh-util.hpp>
#include <mpi.h>

std::vector<uint64_t> SortMorton(std::vector<double> &coords, uint64_t count);
void ConvertToIntegerCoords(std::vector<double> &coords,uint64_t count, std::vector<uint64_t>& coords_integer_out);
bool MortonCompare(uint64_t x1, uint64_t y1, uint64_t z1, uint64_t x2, uint64_t y2, uint64_t z2 );

/**
 * morton encoding with magic bits
 * taken from https://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
*/
// method to seperate bits from a given integer 3 positions apart
inline uint64_t splitBy3(uint64_t a){
    uint64_t x = a & 0x1fffff; // we only look at the first 21 bits
    x = (x | x << 32) & 0x1f00000000ffff; // shift left 32 bits, OR with self, and 00011111000000000000000000000000000000001111111111111111
    x = (x | x << 16) & 0x1f0000ff0000ff; // shift left 32 bits, OR with self, and 00011111000000000000000011111111000000000000000011111111
    x = (x | x << 8) & 0x100f00f00f00f00f; // shift left 32 bits, OR with self, and 0001000000001111000000001111000000001111000000001111000000000000
    x = (x | x << 4) & 0x10c30c30c30c30c3; // shift left 32 bits, OR with self, and 0001000011000011000011000011000011000011000011000011000100000000
    x = (x | x << 2) & 0x1249249249249249;
    return x;
}

inline uint64_t mortonEncode_magicbits(uint64_t x, uint64_t y, uint64_t z){
    uint64_t answer = 0;
    answer |= splitBy3(x) | splitBy3(y) << 1 | splitBy3(z) << 2;
    return answer;
}
/* * * **/

template <class T>
void SetMortonEncoding(std::vector<T> &elements, ElementType element_type, MPI_Comm comm){
    int my_rank;
    MPI_Comm_rank(comm, &my_rank);
    
    double bounding_box_local_min = elements[0].x;
    double bounding_box_local_max = elements[0].x;
    for (size_t i = 0; i < elements.size(); i++)
    {
        bounding_box_local_min = std::min({bounding_box_local_min, elements[i].x, elements[i].y, elements[i].z});
        bounding_box_local_max = std::max({bounding_box_local_max, elements[i].x, elements[i].y, elements[i].z});

    }
    double bounding_box_global_min;
    double bounding_box_global_max;
    MPI_Allreduce(&bounding_box_local_min, &bounding_box_global_min, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&bounding_box_local_max, &bounding_box_global_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    uint64_t levels = 20;
    double spacing = (bounding_box_global_max - bounding_box_global_min)/((double)(1<<levels));
    print_log("[", my_rank, "]:", "bounding_box_global_min = ", bounding_box_global_min);
    print_log("[", my_rank, "]:", "bounding_box_global_max = ", bounding_box_global_max);



    for (size_t element_i = 0; element_i < elements.size(); element_i++)
    {
        uint64_t x = (uint64_t)((elements[element_i].x - bounding_box_global_min)/spacing);
        uint64_t y = (uint64_t)((elements[element_i].y - bounding_box_global_min)/spacing);
        uint64_t z = (uint64_t)((elements[element_i].z - bounding_box_global_min)/spacing);

        elements[element_i].morton_encoding = mortonEncode_magicbits(x,y,z);
    }
    

    // MPI_Reduce()
}
#endif