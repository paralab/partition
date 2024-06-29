#ifndef MESH_UTIL_H
#define MESH_UTIL_H

#include <vector>
#include <string>
// #include "../graph/graph.hpp"
#include "mpi.h"
// #include "mesh-util.hpp"
#include "../util/util.hpp"
#include "../usort/dtypes.h"

#include <unordered_map>

struct TetElementWithFacesNodes {
    uint64_t element_tag;
    uint64_t global_idx;
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[4];
    uint64_t node_tags[4];
    bool operator==(const TetElementWithFacesNodes& other) const {
        return morton_encoding == other.morton_encoding;
    }
    bool operator<(const TetElementWithFacesNodes& other) const {
        return morton_encoding < other.morton_encoding;
    }
    bool operator<=(const TetElementWithFacesNodes& other) const {
        return morton_encoding <= other.morton_encoding;
    }
    bool operator>=(const TetElementWithFacesNodes& other) const {
        return morton_encoding >= other.morton_encoding;
    }
    bool operator>(const TetElementWithFacesNodes& other) const {
        return morton_encoding > other.morton_encoding;
    }
};
// Overloading the << operator for TetElementWithFacesNodes
std::ostream& operator<<(std::ostream& os, const TetElementWithFacesNodes& obj);

struct HexElementWithFacesNodes {
    uint64_t element_tag;
    uint64_t global_idx;
    double x;
    double y;
    double z;
    uint64_t morton_encoding;
    uint64_t face_tags[6];
    uint64_t node_tags[8];
    bool operator==(const HexElementWithFacesNodes& other) const {
        return morton_encoding == other.morton_encoding;
    }
    bool operator<(const HexElementWithFacesNodes& other) const {
        return morton_encoding < other.morton_encoding;
    }
    bool operator<=(const HexElementWithFacesNodes& other) const {
        return morton_encoding <= other.morton_encoding;
    }
    bool operator>=(const HexElementWithFacesNodes& other) const {
        return morton_encoding >= other.morton_encoding;
    }
    bool operator>(const HexElementWithFacesNodes& other) const {
        return morton_encoding > other.morton_encoding;
    }
};

// Overloading the << operator for HexElementWithFacesNodes
std::ostream& operator<<(std::ostream& os, const HexElementWithFacesNodes& obj);



struct ElementWithFace
{
    uint64_t element_tag;
    uint64_t global_idx;
    uint64_t face_tag;
    bool operator==(const ElementWithFace& other) const {
        return face_tag == other.face_tag;
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
        return global_idx == other.global_idx;
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


/**
 * node_tag and its current MPI rank
*/
struct NodeLocationPair
{
    uint64_t node_tag;
    int location_rank;

    bool operator==(const NodeLocationPair& other) const {
        return node_tag == other.node_tag;
    }
    bool operator<(const NodeLocationPair& other) const {
        return node_tag < other.node_tag;
    }
    bool operator<=(const NodeLocationPair& other) const {
        return node_tag <= other.node_tag;
    }
    bool operator>=(const NodeLocationPair& other) const {
        return node_tag >= other.node_tag;
    }
    bool operator>(const NodeLocationPair& other) const {
        return node_tag > other.node_tag;
    }

};

std::ostream& operator<<(std::ostream& os, const NodeLocationPair& obj);



struct NodeNewInfo
{
    uint64_t node_tag;
    int location_rank;
    int owner_rank;
    uint64_t global_idx;

    bool operator==(const NodeNewInfo& other) const {
        return node_tag == other.node_tag;
    }
    bool operator<(const NodeNewInfo& other) const {
        return node_tag < other.node_tag;
    }
    bool operator<=(const NodeNewInfo& other) const {
        return node_tag <= other.node_tag;
    }
    bool operator>=(const NodeNewInfo& other) const {
        return node_tag >= other.node_tag;
    }
    bool operator>(const NodeNewInfo& other) const {
        return node_tag > other.node_tag;
    }

};

std::ostream& operator<<(std::ostream& os, const NodeNewInfo& obj);

namespace par
{
    template <>
    class Mpi_datatype<TetElementWithFacesNodes> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "TetElementWithFacesNodes"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[8] = {1, 1, 1, 1, 1, 1, 4, 4};
                MPI_Datatype types[8] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[8];
                offsets[0] = offsetof(TetElementWithFacesNodes, element_tag);
                offsets[1] = offsetof(TetElementWithFacesNodes, global_idx);
                offsets[2] = offsetof(TetElementWithFacesNodes, x);
                offsets[3] = offsetof(TetElementWithFacesNodes, y);
                offsets[4] = offsetof(TetElementWithFacesNodes, z);
                offsets[5] = offsetof(TetElementWithFacesNodes, morton_encoding);
                offsets[6] = offsetof(TetElementWithFacesNodes, face_tags);
                offsets[7] = offsetof(TetElementWithFacesNodes, node_tags);



                MPI_Type_create_struct(8, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<HexElementWithFacesNodes> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "HexElementWithFacesNodes"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[8] = {1, 1, 1, 1, 1, 1, 6, 8};
                MPI_Datatype types[8] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[8];
                offsets[0] = offsetof(HexElementWithFacesNodes, element_tag);
                offsets[1] = offsetof(HexElementWithFacesNodes, global_idx);
                offsets[2] = offsetof(HexElementWithFacesNodes, x);
                offsets[3] = offsetof(HexElementWithFacesNodes, y);
                offsets[4] = offsetof(HexElementWithFacesNodes, z);
                offsets[5] = offsetof(HexElementWithFacesNodes, morton_encoding);
                offsets[6] = offsetof(HexElementWithFacesNodes, face_tags);
                offsets[7] = offsetof(HexElementWithFacesNodes, node_tags);



                MPI_Type_create_struct(8, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<ElementWithFace> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "ElementWithFace"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[3] = {1, 1, 1};
                MPI_Datatype types[3] = {MPI_UINT64_T, MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[3];
                offsets[0] = offsetof(ElementWithFace, element_tag);
                offsets[1] = offsetof(ElementWithFace, global_idx);
                offsets[2] = offsetof(ElementWithFace, face_tag);


                MPI_Type_create_struct(3, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<ElementWithTag> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "ElementWithTag"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[2] = {1, 1};
                MPI_Datatype types[2] = {MPI_UINT64_T, MPI_UINT64_T};
                MPI_Aint offsets[2];
                offsets[0] = offsetof(ElementWithTag, element_tag);
                offsets[1] = offsetof(ElementWithTag, global_idx);



                MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<ElementWithCoord> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "ElementWithCoord"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[5] = {1, 1, 1, 1, 1};
                MPI_Datatype types[5] = {MPI_UINT64_T, MPI_UINT64_T, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE};
                MPI_Aint offsets[5];
                offsets[0] = offsetof(ElementWithCoord, element_tag);
                offsets[1] = offsetof(ElementWithCoord, global_idx);
                offsets[2] = offsetof(ElementWithCoord, x);
                offsets[3] = offsetof(ElementWithCoord, y);
                offsets[4] = offsetof(ElementWithCoord, z);




                MPI_Type_create_struct(5, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_pairtype<ElementWithTag,ElementWithTag> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "std::pair<ElementWithTag,ElementWithTag>"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                
                first = false;
                MPI_Datatype inner_type = Mpi_datatype<ElementWithTag>::value();
                
                int second_value_offset;
                MPI_Type_size(inner_type, &second_value_offset);
                int block_lengths[2] = {1, 1};
                MPI_Datatype types[2] = {inner_type, inner_type};
                MPI_Aint offsets[2];
                offsets[0] = 0;
                offsets[1] = static_cast<MPI_Aint>(second_value_offset);


                MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }       



            return custom_mpi_type;
        }
    };

    template <>
    class Mpi_datatype<NodeLocationPair> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "NodeLocationPair"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[2] = {1, 1};
                MPI_Datatype types[2] = {MPI_UINT64_T, MPI_INT};
                MPI_Aint offsets[2];
                offsets[0] = offsetof(NodeLocationPair, node_tag);
                offsets[1] = offsetof(NodeLocationPair, location_rank);



                MPI_Type_create_struct(2, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }

            return custom_mpi_type;
        }
    };


    template <>
    class Mpi_datatype<NodeNewInfo> {

	  /** 
          @return the MPI_Datatype for the C++ datatype "NodeNewInfo"
         **/
        public:
        static MPI_Datatype value() {
            static bool         first = true;
            static MPI_Datatype custom_mpi_type;

            if (first)
            {
                first = false;
                int block_lengths[4] = {1, 1, 1, 1};
                MPI_Datatype types[4] = {MPI_UINT64_T, MPI_INT, MPI_INT, MPI_UINT64_T};
                MPI_Aint offsets[4];
                offsets[0] = offsetof(NodeNewInfo, node_tag);
                offsets[1] = offsetof(NodeNewInfo, location_rank);
                offsets[2] = offsetof(NodeNewInfo, owner_rank);
                offsets[3] = offsetof(NodeNewInfo, global_idx);

                MPI_Type_create_struct(4, block_lengths, offsets, types, &custom_mpi_type);
                MPI_Type_commit(&custom_mpi_type);
            }

            return custom_mpi_type;
        }
    };

}


enum ElementType { TET=4, HEX=5 };


struct DistributionStatus
{
    int return_code;
    int time_us;
};



// Graph GmshGetElementGraph(const std::string &mesh_file_path , std::vector<double>& elem_coordinates_out, std::vector<size_t>& elem_tags);

ElementType GetElementType(const std::string &part_file_prefix, MPI_Comm comm);

template <class T>
void GetElementsWithFacesNodesCentroids(const std::string &mesh_file_path, std::vector<T> &elements_out,
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



template <class T>
DistributionStatus Redistribute(std::vector<T> &elements_in, std::vector<uint16_t>& labeling, std::vector<T> &elements_out, MPI_Comm comm);

// template <class T>
// void GetOwnNodes(const std::vector<T> &local_elements, ElementType element_type, std::vector<u_int64_t>& own_node_tags_out, MPI_Comm comm);
template <class T>
void GetNodetagToGlobalIdx(const std::vector<T> &local_elements, ElementType element_type, 
        std::unordered_map<uint64_t, uint64_t>& mapping_out, uint64_t& global_count_out, uint64_t& local_start_out, uint64_t& local_end_out, MPI_Comm comm);
#include "mesh-util.tcc"

#endif