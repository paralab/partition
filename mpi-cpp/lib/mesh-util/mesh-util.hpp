#ifndef MESH_UTIL_H
#define MESH_UTIL_H

#include <vector>
#include <string>
#include <graph.hpp>
#include <mpi.h>
#include <mesh-util.hpp>

struct TetElementWithFaces {
    uint64_t element_tag;
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[4];
    bool operator<(const TetElementWithFaces& other) const {
        return morton_encoding < other.morton_encoding;
    }
    bool operator<=(const TetElementWithFaces& other) const {
        return morton_encoding <= other.morton_encoding;
    }
    bool operator>=(const TetElementWithFaces& other) const {
        return morton_encoding >= other.morton_encoding;
    }
    bool operator>(const TetElementWithFaces& other) const {
        return morton_encoding > other.morton_encoding;
    }
};
// Overloading the << operator for TetElementWithFaces
std::ostream& operator<<(std::ostream& os, const TetElementWithFaces& obj);

struct HexElementWithFaces {
    uint64_t element_tag;
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[6];
    bool operator<(const HexElementWithFaces& other) const {
        return morton_encoding < other.morton_encoding;
    }
    bool operator<=(const HexElementWithFaces& other) const {
        return morton_encoding <= other.morton_encoding;
    }
    bool operator>=(const HexElementWithFaces& other) const {
        return morton_encoding >= other.morton_encoding;
    }
    bool operator>(const HexElementWithFaces& other) const {
        return morton_encoding > other.morton_encoding;
    }
};

// Overloading the << operator for HexElementWithFaces
std::ostream& operator<<(std::ostream& os, const HexElementWithFaces& obj);

enum ElementType { TET=4, HEX=5 };


Graph GmshGetElementGraph(const std::string &mesh_file_path , std::vector<double>& elem_coordinates_out, std::vector<size_t>& elem_tags);

ElementType GetElementType(const std::string &mesh_file_path, MPI_Comm comm);

template <class T>
void GetElementsWithFacesCentroids(const std::string &mesh_file_path, std::vector<T> &elements_out,
                          ElementType element_type, MPI_Comm comm);

#include "mesh-util.tcc"

#endif