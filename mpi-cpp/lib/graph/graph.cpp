// #include "graph.hpp"
// #include "../util/util.hpp"
// #include <limits>
// #include <stdexcept>
// #include <cassert>

// #include <omp.h>

// #define NO_LABEL SIZE_MAX

// Graph::Graph(/* args */)
// {
//     this->graph_size = 0;
//     this->infinity = std::numeric_limits<unsigned long>::max();
// }

// Graph::~Graph()
// {
// }

// void Graph::AddVertex(vertex_t v)
// {
//     auto index = this->vertex_to_index.find(v);
//     if (index != this->vertex_to_index.end())
//     {
//         return;
//     }

//     auto new_index = this->graph_size++;
//     this->vertex_to_index[v] = new_index;
//     this->index_to_vetex.push_back(v);
//     this->adj_list.push_back({});
//     // this->bfs_status.push_back(infinity);
//     // this->multi_bfs_distances.push_back(infinity);
//     // this->multi_bfs_labels.push_back(NO_LABEL);

// }

// std::vector<vertex_t>& Graph::GetVertices(){
//     return this->index_to_vetex;
// }

// bool Graph::AddEdge(vertex_t v, vertex_t u)
// {
//     auto index_v = this->vertex_to_index.find(v);
//     auto index_u = this->vertex_to_index.find(u);
//     if (index_v == this->vertex_to_index.end() || index_u == this->vertex_to_index.end())
//     {
//         std::cout << "invalid u v" << u << " " << v;
//         throw std::invalid_argument( "u or v does not exist in the graph" );
//         return false;
//     }
//     this->adj_list[index_u->second].insert(index_v->second);
//     this->adj_list[index_v->second].insert(index_u->second);
//     return true;
// }

// std::unordered_set<vertex_t> Graph::GetNeighbors(vertex_t v){
//     auto index_v = this->vertex_to_index.find(v);
//     if (index_v == this->vertex_to_index.end())
//     {
//         throw std::invalid_argument( "invalid vertex" );
//     }
    
//     return this->adj_list[index_v->second];
// }

// bool Graph::FindVertex(vertex_t v){
//     return this->vertex_to_index.find(v) != this->vertex_to_index.end();
// }

// void Graph::InitSingleBFS(vertex_t BFS_seed)
// {
//     this->bfs_status.resize(this->graph_size);
//     std::fill(this->bfs_status.begin(), this->bfs_status.end(), this->infinity);
//     this->bfs_status[this->vertex_to_index[BFS_seed]] = 0;
//     return;
// }

// void Graph::InitMultiBFS(std::vector<vertex_t>& seeds, uint64_t count)
// {
//     assert(seeds.size() == count);
//     this->multi_bfs_distances.resize(this->graph_size);
//     std::fill(this->multi_bfs_distances.begin(), this->multi_bfs_distances.end(), this->infinity);

//     this->multi_bfs_labels.resize(this->graph_size);
//     std::fill(this->multi_bfs_labels.begin(), this->multi_bfs_labels.end(), NO_LABEL);

//     for (size_t i = 0; i < count; i++)
//     {
//         this->multi_bfs_distances[this->vertex_to_index[seeds[i]]] = 0;
//         this->multi_bfs_labels[this->vertex_to_index[seeds[i]]] = i;
//     }    
//     return;
// }

// void Graph::RunBFSToStable()
// {
//     bool is_not_stable = true;
//     std::vector<unsigned long> bfs_status_new_temp(this->graph_size);

//     while (is_not_stable)
//     {

//         is_not_stable = false;
//         // #pragma omp parallel
//         {
//             // #pragma omp for
//             for (unsigned long v_i = 0; v_i < this->graph_size; v_i++)
//             {
//                 bfs_status_new_temp[v_i] = NULL;
//                 auto best_distance = this->bfs_status[v_i];
//                 for (auto neighbor_i : this->adj_list[v_i])
//                 {
//                     if (this->bfs_status[neighbor_i] == this->infinity)
//                     {
//                         continue;
//                     }

//                     if (best_distance > (this->bfs_status[neighbor_i] + 1))
//                     {
//                         best_distance = this->bfs_status[neighbor_i] + 1;
//                         bfs_status_new_temp[v_i] = best_distance;
//                         // #pragma omp critical
//                         is_not_stable = true;
//                         break;
//                     }
//                 }
//             }

//             // #pragma omp for
//             for (unsigned long v_i = 0; v_i < this->graph_size; v_i++)
//             {
//                 if (bfs_status_new_temp[v_i] != NULL)
//                 {
//                     this->bfs_status[v_i] = bfs_status_new_temp[v_i];
//                 }
//             }
//         }
//     }
// }

// void Graph::RunMultiBFSToStable(){
//     bool is_not_stable = true;
//     std::vector<uint64_t> multi_bfs_distances_new_temp(this->graph_size);
//     std::vector<uint64_t> multi_bfs_labels_new_temp(this->graph_size);
//     std::vector<bool> vertex_is_not_stable(this->graph_size,false);


//     while (is_not_stable)
//     {
//         // print_log(VectorToString(multi_bfs_distances));

//         is_not_stable = false;
//         std::fill(vertex_is_not_stable.begin(), vertex_is_not_stable.end(), false);
//         #pragma omp parallel
//         {
//             // print_log("thread count ", omp_get_num_threads());
//             #pragma omp for
//             for (unsigned long v_i = 0; v_i < this->graph_size; v_i++)
//             {
//                 // bfs_status_new_temp[v_i] = NULL;
//                 auto best_distance = this->multi_bfs_distances[v_i];
//                 auto best_label = this->multi_bfs_labels[v_i];
//                 // print_log(best_distance, best_label);
//                 for (auto neighbor_i : this->adj_list[v_i])
//                 {
//                     if (this->multi_bfs_labels[neighbor_i] == NO_LABEL)
//                     {
//                         continue;
//                     }
//                     // print_log(multi_bfs_labels[neighbor_i]);

//                     if (best_distance > (this->multi_bfs_distances[neighbor_i] + 1))
//                     {
//                         best_distance = this->multi_bfs_distances[neighbor_i] + 1;
//                         best_label = this->multi_bfs_labels[neighbor_i];
//                         // #pragma omp critical
//                         // is_not_stable = true;
//                         vertex_is_not_stable[v_i] = true;

//                     }
//                 }
//                 multi_bfs_distances_new_temp[v_i] = best_distance;
//                 multi_bfs_labels_new_temp[v_i] = best_label;

//             }
//             #pragma omp for
//             for (unsigned long v_i = 0; v_i < this->graph_size; v_i++){
//                 this->multi_bfs_distances[v_i]=multi_bfs_distances_new_temp[v_i];
//                 this->multi_bfs_labels[v_i]=multi_bfs_labels_new_temp[v_i];
//             }
//             #pragma omp for reduction(||:is_not_stable)
//             for (unsigned long v_i = 0; v_i < this->graph_size; v_i++) {
//                 is_not_stable = is_not_stable|| vertex_is_not_stable[v_i]; // Perform bitwise OR operation
//             }

//             // print_log(VectorToString(vertex_is_not_stable));


//         }
//     }
// }

// std::vector<unsigned long> & Graph::GetBFSState(){
//     return this->bfs_status;
// }

// std::vector<uint64_t>& Graph::GetMultiBFSLabels(){
//     return this->multi_bfs_labels;
// }


// unsigned long Graph::GETBFSValue(vertex_t v){
//     auto index_v = this->vertex_to_index.find(v);
//     if (index_v == this->vertex_to_index.end())
//     {
//         throw std::invalid_argument( "invalid vertex" );
//     }
    
//     return this->bfs_status[index_v->second];
// }

// unsigned long Graph::GetSize(){
//     return this->graph_size;
// }

// void Graph::Print()
// {
//     std::ostringstream output;
//     unsigned long index = 0;
//     for (std::unordered_set<vertex_t> neighbors : this->adj_list)
//     {
//         output << this->index_to_vetex[index] << " -> ";
//         for (auto neighbor : neighbors)
//         {
//             output << this->index_to_vetex[neighbor] << " ";
//         }
//         output << "\n";
//         index++;
//     }
//     std::cout << output.str();
// }

// void Graph::PrintSingleBFS()
// {
//     std::ostringstream output;
//     for (unsigned long v_i = 0; v_i < this->graph_size; v_i++)
//     {
//         output << this->index_to_vetex[v_i] << " = ";
//         output << this->bfs_status[v_i] << "\n";
//     }
//     std::cout << output.str();
// }


// std::vector<uint64_t> Graph::GetCSR_xadj(){
//     std::vector<uint64_t> xadj(this->graph_size + 1);
//     xadj[0]=0;
//     for (size_t i = 0; i < this->graph_size; i++)
//     {
//         xadj[i+1] = xadj[i] + this->adj_list[i].size();
//     }
//     return xadj;    
// }

// std::vector<uint64_t> Graph::GetCSR_adjncy(){
//     std::vector<uint64_t> adjncy = {};
//     for (size_t i = 0; i < this->graph_size; i++)
//     {
//         for (auto neighbor_i : this->adj_list[i]){
//             adjncy.push_back(neighbor_i);
//         }
//     }
//     return adjncy;


// }