#ifndef MESH_UTIL_H
#define MESH_UTIL_H

#include <vector>
#include <string>
#include <graph.hpp>

Graph GmshGetElementGraph(const std::string &mesh_file_path , std::vector<double>& elem_coordinates_out, std::vector<size_t>& elem_tags);

#endif