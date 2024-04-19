#ifndef VTK_UTIL_H
#define VTK_UTIL_H

#include <vector>
#include <string>


void PointsWithScalarsToVtk(std::vector<double>& point_coords, std::vector<double>& scalars, uint64_t count, std::string out_file_name);
void PointsToVtk(std::vector<double>& point_coords, uint64_t count, std::string out_file_name);

#endif