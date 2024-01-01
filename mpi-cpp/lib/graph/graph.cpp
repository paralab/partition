#include "graph.hpp"
#include <limits>
#include <stdexcept>

Graph::Graph(/* args */)
{
    this->graph_size = 0;
    this->infinity = std::numeric_limits<unsigned long>::max();
}

Graph::~Graph()
{
}

void Graph::AddVertex(vertex_t v)
{
    auto index = this->vertex_to_index.find(v);
    if (index != this->vertex_to_index.end())
    {
        return;
    }

    auto new_index = this->graph_size++;
    this->vertex_to_index[v] = new_index;
    this->index_to_vetex.push_back(v);
    this->adj_list.push_back({});
    this->bfs_status.push_back(infinity);
}

std::vector<vertex_t>& Graph::GetVertices(){
    return this->index_to_vetex;
}

bool Graph::AddEdge(vertex_t v, vertex_t u)
{
    auto index_v = this->vertex_to_index.find(v);
    auto index_u = this->vertex_to_index.find(u);
    if (index_v == this->vertex_to_index.end() || index_u == this->vertex_to_index.end())
    {
        std::cout << "invalid u v" << u << " " << v;
        throw std::invalid_argument( "u or v does not exist in the graph" );
        return false;
    }
    this->adj_list[index_u->second].insert(index_v->second);
    this->adj_list[index_v->second].insert(index_u->second);
    return true;
}

std::unordered_set<vertex_t> Graph::GetNeighbors(vertex_t v){
    auto index_v = this->vertex_to_index.find(v);
    if (index_v == this->vertex_to_index.end())
    {
        throw std::invalid_argument( "invalid vertex" );
    }
    
    return this->adj_list[index_v->second];
}

bool Graph::FindVertex(vertex_t v){
    return this->vertex_to_index.find(v) != this->vertex_to_index.end();
}

void Graph::InitSingleBFS(vertex_t BFS_seed)
{
    this->bfs_status[this->vertex_to_index[BFS_seed]] = 0;
    return;
}

void Graph::RunBFSToStable()
{
    bool is_not_stable = true;
    std::vector<unsigned long> bfs_status_new_temp(this->graph_size);

    while (is_not_stable)
    {

        is_not_stable = false;
        // #pragma omp parallel
        {
            // #pragma omp for
            for (unsigned long v_i = 0; v_i < this->graph_size; v_i++)
            {
                bfs_status_new_temp[v_i] = NULL;
                auto best_distance = this->bfs_status[v_i];
                for (auto neighbor_i : this->adj_list[v_i])
                {
                    if (this->bfs_status[neighbor_i] == this->infinity)
                    {
                        continue;
                    }

                    if (best_distance > (this->bfs_status[neighbor_i] + 1))
                    {
                        best_distance = this->bfs_status[neighbor_i] + 1;
                        bfs_status_new_temp[v_i] = best_distance;
                        // #pragma omp critical
                        is_not_stable = true;
                        break;
                    }
                }
            }

            // #pragma omp for
            for (unsigned long v_i = 0; v_i < this->graph_size; v_i++)
            {
                if (bfs_status_new_temp[v_i] != NULL)
                {
                    this->bfs_status[v_i] = bfs_status_new_temp[v_i];
                }
            }
        }
    }
}

std::vector<unsigned long> & Graph::GetBFSState(){
    return this->bfs_status;
}


unsigned long Graph::GETBFSValue(vertex_t v){
    auto index_v = this->vertex_to_index.find(v);
    if (index_v == this->vertex_to_index.end())
    {
        throw std::invalid_argument( "invalid vertex" );
    }
    
    return this->bfs_status[index_v->second];
}

unsigned long Graph::GetSize(){
    return this->graph_size;
}

void Graph::Print()
{
    std::ostringstream output;
    unsigned long index = 0;
    for (std::unordered_set<vertex_t> neighbors : this->adj_list)
    {
        output << this->index_to_vetex[index] << " -> ";
        for (auto neighbor : neighbors)
        {
            output << this->index_to_vetex[neighbor] << " ";
        }
        output << "\n";
        index++;
    }
    std::cout << output.str();
}

void Graph::PrintSingleBFS()
{
    std::ostringstream output;
    for (unsigned long v_i = 0; v_i < this->graph_size; v_i++)
    {
        output << this->index_to_vetex[v_i] << " = ";
        output << this->bfs_status[v_i] << "\n";
    }
    std::cout << output.str();
}
