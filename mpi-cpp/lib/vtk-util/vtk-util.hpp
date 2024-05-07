#ifndef VTK_UTIL_H
#define VTK_UTIL_H

#include <vector>
#include <string>

#include "../mesh-util/mesh-util.hpp"


void PointsWithPartitionsToVtk(std::vector<double>& point_coords, std::vector<uint64_t>& partitions, uint64_t count, std::string out_file_name);
void ElementsWithPartitionsToVtk(std::vector<ElementWithCoord>& elements, std::vector<uint16_t>& partitions, uint64_t count, std::string out_file_name);
void PointsToVtk(std::vector<double>& point_coords, uint64_t count, std::string out_file_name);
#endif