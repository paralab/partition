#ifndef UTIL_H
#define UTIL_H


#include <iostream>
#include <sstream>
#include <utility>
#include <algorithm>
#include <unordered_map>
#include <vector>




void AssignPartitionLabelsInOrder(std::vector<uint64_t> &ordering, uint64_t count, uint64_t partition_count, std::vector<uint64_t> &labels_out);
void GetSamplesFromOrdered(std::vector<uint64_t> &order, std::vector<uint64_t> &input_arr, uint64_t sample_count, std::vector<uint64_t> &samples_out);

// template <typename T> std::string VectorToString(std::vector<T> vec);

// // Explicit instantiation of the specialized template
// extern template std::string VectorToString<TetElementWithFacesNodes>(std::vector<TetElementWithFacesNodes> vec);
// extern template std::string VectorToString<HexElementWithFacesNodes>(std::vector<HexElementWithFacesNodes> vec);

// Variadic template function to mimic std::cout with spaces between arguments
template<typename T, typename... Args>
void print_log(const T& first, const Args&... args) {
    std::ostringstream oss;
    oss << first;  // Output the first argument directly
    // Use fold expression to concatenate all arguments with spaces in between
    ((oss << ' ' << args), ...);
    // Output the concatenated string
    std::cout << oss.str() <<std::endl;
}

template<typename T, typename U>
std::ostream& operator<<(std::ostream& os, const std::pair<T, U>& p) {
    os << "(" << p.first << ", " << p.second << ")";
    return os;
}

template <typename T> std::string VectorToString(std::vector<T> vec){
    std::ostringstream output;
    output << "[ ";
    for (auto element : vec)
    {
        output << element << ", ";
    }
    output << "]\n";
    return output.str();
}
template <> std::string VectorToString(std::vector<uint8_t> vec);

template <typename T> std::vector<int> GetDisplacementsFromCounts(std::vector<T>& counts){
    std::vector<int> displacements(counts.size());
    displacements[0]=0;
    for (size_t i = 1; i < counts.size(); i++)
    {
        displacements[i] = (int)(displacements[i-1] + counts[i-1]);
    }
    return displacements;
}

// https://stackoverflow.com/a/34937216
// templated function to get the max element from a std::unordered_map
template<typename KeyType, typename ValueType> 
inline std::pair<KeyType,ValueType> get_max( const std::unordered_map<KeyType,ValueType>& x ) {
  using pairtype=std::pair<KeyType,ValueType>; 
  return *std::max_element(x.begin(), x.end(), [] (const pairtype & p1, const pairtype & p2) {
        return p1.second < p2.second;
  }); 
}

void ExportMetricsToJson(
    std::string mesh_file, int file_idx, int run_idx, int partition_count, uint64_t global_vertex_count,
    int graph_setup_time,
    std::vector<uint32_t>& sfc_partition_sizes, std::vector<uint32_t>& sfc_partition_boundaries, int sfc_partition_time, int sfc_mat_assembly_time, int sfc_matvec_time,
    std::vector<uint32_t>& bfs_partition_sizes, std::vector<uint32_t>& bfs_partition_boundaries, int bfs_labeling_time, int bfs_redistribution_time, int bfs_mat_assembly_time, int bfs_matvec_time,
    std::vector<uint32_t>& parmetis_partition_sizes, std::vector<uint32_t>& parmetis_partition_boundaries, int parmetis_labeling_time, int parmetis_redistribution_time, int parmetis_mat_assembly_time, int parmetis_matvec_time,
    std::vector<uint32_t>& ptscotch_partition_sizes, std::vector<uint32_t>& ptscotch_partition_boundaries, int ptscotch_labeling_time, int ptscotch_redistribution_time, int ptscotch_mat_assembly_time, int ptscotch_matvec_time,
    std::string metrics_out_file_path);


// void ExportMetricsToPandasJson(
//     std::string mesh_file, int file_idx, int run_idx, int partition_count, uint64_t global_vertex_count,
//     int graph_setup_time,
//     std::vector<uint32_t>& sfc_partition_sizes, std::vector<uint32_t>& sfc_partition_boundaries, int sfc_partition_time, int sfc_mat_assembly_time, int sfc_matvec_time,
//     std::vector<uint32_t>& bfs_partition_sizes, std::vector<uint32_t>& bfs_partition_boundaries, int bfs_labeling_time, int bfs_redistribution_time, int bfs_mat_assembly_time, int bfs_matvec_time,
//     std::vector<uint32_t>& parmetis_partition_sizes, std::vector<uint32_t>& parmetis_partition_boundaries, int parmetis_labeling_time, int parmetis_redistribution_time, int parmetis_mat_assembly_time, int parmetis_matvec_time,
//     std::vector<uint32_t>& ptscotch_partition_sizes, std::vector<uint32_t>& ptscotch_partition_boundaries, int ptscotch_labeling_time, int ptscotch_redistribution_time, int ptscotch_mat_assembly_time, int ptscotch_matvec_time,
//     std::string metrics_out_file_path) ;

struct PartitionStatus
{
    int return_code;
    int time_us;
};

#endif
