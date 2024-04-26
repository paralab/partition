#include <vector>
#include <cstdint>
#include <cassert>
#include <numeric>
#include <algorithm>

#include "../util/util.hpp"
#include "sfc.hpp"

std::vector<uint64_t> SortMorton(std::vector<double> &coords, uint64_t count)
{
    assert(coords.size() == count * 3);

    std::vector<uint64_t> coords_integer(count * 3);
    ConvertToIntegerCoords(coords, count, coords_integer);
    // print_log(VectorToString(coords_integer));
    std::vector<uint64_t> sorted_indices(count);
    std::iota(std::begin(sorted_indices), std::end(sorted_indices), 0); // Fill with 0, 1, ..., count-1.
    stable_sort(sorted_indices.begin(), sorted_indices.end(),
                [&coords_integer](size_t i1, size_t i2){ 
                    return MortonCompare(coords_integer[i1*3], coords_integer[i1*3+1], coords_integer[i1*3+2],
                                            coords_integer[i2*3], coords_integer[i2*3+1], coords_integer[i2*3+2]);
                });
    
    return sorted_indices;
}

void ConvertToIntegerCoords(std::vector<double> &coords,uint64_t count, std::vector<uint64_t>& coords_integer_out){
    assert(coords.size() == count*3);
    assert(coords_integer_out.size() == count*3);

    uint64_t levels = 19;

    /**
     * regular method
    */
    std::vector<double> bounding_box = {coords[0], coords[0]};  //min max
    for (uint64_t i = 0; i < count*3; i++)
    {
        bounding_box[0] = std::min(coords[i], bounding_box[0]);
        bounding_box[1] = std::max(coords[i], bounding_box[1]);
       
    }
    print_log("bounding box done");
    print_log("bounding box:", VectorToString(bounding_box));
    double spacing = (bounding_box[1] - bounding_box[0])/((double)(1<<levels));

    print_log("spacing done");
    print_log("spacing:", spacing);
    for (uint64_t i = 0; i < count*3; i++)
    {
        coords_integer_out[i] = (uint64_t)((coords[i] - bounding_box[0])/spacing) + 1;     // (coordinate - min)/spacing
      
    }

    return;


    /**
     * scaled method
    */

    // std::vector<double> bounding_box = {coords[0], coords[0], coords[1], coords[1], coords[2], coords[2]};  //xmin xmax ymin ymax zmin zmax
    // for (uint64_t i = 0; i < count; i++)
    // {
    //     for (int dim_i = 0; dim_i < 3; dim_i++)
    //     {
    //         bounding_box[dim_i*2] = std::min(bounding_box[dim_i*2],coords[i*3 + dim_i]);
    //         bounding_box[dim_i*2+1] = std::max(bounding_box[dim_i*2+1],coords[i*3 + dim_i]);
    //     }        
    // }
    // print_log("bounding box done");
    // print_log("bounding box:", VectorToString(bounding_box));
    // std::vector<double> spacing = {0,0,0};

    // for (int dim_i = 0; dim_i < 3; dim_i++)
    // {
    //     spacing[dim_i] = (bounding_box[dim_i*2 + 1] - bounding_box[dim_i*2])/((double)(1<<levels));
    // }
    // print_log("spacing done");
    // print_log("spacing:", VectorToString(spacing));
    // for (uint64_t i = 0; i < count; i++)
    // {
    //     for (int dim_i = 0; dim_i < 3; dim_i++)
    //     {
    //         coords_integer_out[i*3+dim_i] = (uint64_t)((coords[i*3+dim_i] - bounding_box[dim_i*2])/spacing[dim_i]) + 1;     // (coordinate - min)/spacing
    //     }        
    // }
    // return;
    
}

bool MortonCompare(uint64_t x1, uint64_t y1, uint64_t z1, uint64_t x2, uint64_t y2, uint64_t z2 ){
    if ((x1==x2) && (y1==y2) && (z1==z2))
    {
        return false;
    }

    uint64_t temp_x = x1^x2;
    uint64_t temp_y = y1^y2;
    uint64_t temp_z = z1^z2;

    uint64_t max_c = temp_z;
    uint64_t y_or_x = temp_y;

    if (y_or_x < temp_x)
    {
        if ((temp_x ^ y_or_x) >= y_or_x)
        {
            y_or_x = temp_x;
        }
    }
    if (max_c < y_or_x)
    {
        if ((max_c ^ y_or_x) >= max_c)
        {
            max_c = y_or_x;
        }
    }
    if (max_c == temp_z)
    {
        return (z1 < z2);
    }
    else if (max_c == temp_y)
    {
        return (y1 < y2);
    }
    else
    {
        return (x1 < x2);
    }
}