#ifndef MESH_UTIL_H
#define MESH_UTIL_H

#include <vector>
#include <string>
#include "../graph/graph.hpp"
#include "mpi.h"
// #include "mesh-util.hpp"
#include "../util/util.hpp"

struct TetElementWithFaces {
    uint64_t element_tag;
    uint64_t global_idx;
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[4];
    bool operator==(const TetElementWithFaces& other) const {
        return morton_encoding == other.morton_encoding;
    }
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
    uint64_t global_idx;
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[6];
    bool operator==(const HexElementWithFaces& other) const {
        return morton_encoding == other.morton_encoding;
    }
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



struct ElementWithFace
{
    uint64_t element_tag;
    uint64_t global_idx;
    uint64_t face_tag;
    bool operator==(const ElementWithFace& other) const {
        return face_tag < other.face_tag;
    }
    bool operator<(const ElementWithFace& other) const {
        return face_tag < other.face_tag;
    }
    bool operator<=(const ElementWithFace& other) const {
        return face_tag <= other.face_tag;
    }
    bool operator>=(const ElementWithFace& other) const {
        return face_tag >= other.face_tag;
    }
    bool operator>(const ElementWithFace& other) const {
        return face_tag > other.face_tag;
    }
};

std::ostream& operator<<(std::ostream& os, const ElementWithFace& obj);

struct ElementWithCoord
{
    uint64_t element_tag;
    uint64_t global_idx;
    double x;
    double y;
    double z;
};
std::ostream& operator<<(std::ostream& os, const ElementWithCoord& obj);


struct ElementWithTag
{
    uint64_t element_tag;
    uint64_t global_idx;
    bool operator==(const ElementWithTag& other) const {
        return global_idx < other.global_idx;
    }
    bool operator<(const ElementWithTag& other) const {
        return global_idx < other.global_idx;
    }
    bool operator<=(const ElementWithTag& other) const {
        return global_idx <= other.global_idx;
    }
    bool operator>=(const ElementWithTag& other) const {
        return global_idx >= other.global_idx;
    }
    bool operator>(const ElementWithTag& other) const {
        return global_idx > other.global_idx;
    }
};
std::ostream& operator<<(std::ostream& os, const ElementWithTag& obj);


enum ElementType { TET=4, HEX=5 };


Graph GmshGetElementGraph(const std::string &mesh_file_path , std::vector<double>& elem_coordinates_out, std::vector<size_t>& elem_tags);

ElementType GetElementType(const std::string &mesh_file_path, MPI_Comm comm);

template <class T>
void GetElementsWithFacesCentroids(const std::string &mesh_file_path, std::vector<T> &elements_out,
                          ElementType element_type, MPI_Comm comm);


template <class T>
void ResolveLocalElementConnectivity(const std::vector<T> &elements, ElementType element_type,
                                std::vector<std::pair<ElementWithTag, ElementWithTag>> &connected_element_pairs_out,
                                std::vector<ElementWithFace> &unconnected_elements_faces_out);


void ResolveBoundaryElementConnectivity(std::vector<ElementWithFace> &unpaired_element_faces,
                                        std::vector<uint64_t> &proc_element_counts,
                                        std::vector<std::pair<ElementWithTag, ElementWithTag>> &boundary_connected_element_pairs_out,
                                        MPI_Comm comm);


void ExtractGhostElements(std::vector<std::pair<ElementWithTag, ElementWithTag>>& boundary_connected_element_pairs,
                          std::vector<uint64_t>& proc_element_counts,
                          std::vector<uint64_t>& proc_element_counts_scanned,
                          std::vector<ElementWithTag>& ghost_elements_out, std::vector<int>& ghost_element_counts_out,
                          MPI_Comm comm);
#include "mesh-util.tcc"

#endif