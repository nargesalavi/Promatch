/*
 *  author: Narges Alavisamani
 *  date:   6 March 2023
 * */

#ifndef PREDECODER_GRAPH_HELPER_h
#define PREDECODER_GRAPH_HELPER_h

#include <iostream>
#include <map>
#include <vector>
#include <unordered_set>
#include "decoding_graph.h"


namespace qpd{

struct CountResult {
    int countANotInB;
    int countBNotInA;
    int countBothAB;
};

void dfs_(qrc::DecodingGraph::Vertex* v, 
        std::unordered_set<qrc::DecodingGraph::Vertex*>& visited, 
        std::vector<qrc::DecodingGraph::Vertex*>& component);

std::vector<std::vector<qrc::DecodingGraph::Vertex*>> 
find_connected_components(std::map<qrc::DecodingGraph::Vertex*, 
                        std::vector<qrc::DecodingGraph::Vertex*>>& adjacency_matrix);

void print_components(std::map<qrc::DecodingGraph::Vertex*, 
                        std::vector<qrc::DecodingGraph::Vertex*>>& adjacency_matrix, qrc::DecodingGraph& decoding_graph);


bool isDuplicate(const std::vector<qrc::DecodingGraph::Vertex*>& vertices, uint detector);

void addUniqueVertex(std::vector<qrc::DecodingGraph::Vertex*>& vertices, qrc::DecodingGraph::Vertex* newVertex);

CountResult countMembersInAB (const std::vector<qrc::DecodingGraph::Vertex*>& a, 
                                            const std::vector<qrc::DecodingGraph::Vertex*>& b);

CountResult countMembersInAB_edges(const std::vector<qrc::DecodingGraph::Edge*>& a,
                             const std::vector<qrc::DecodingGraph::Edge*>& b);

bool compare_edges_by_error_probability(const qrc::DecodingGraph::Edge* edge1, const qrc::DecodingGraph::Edge* edge2);
bool compare_edges_by_future_error_probability(const qrc::DecodingGraph::Edge* edge1, const qrc::DecodingGraph::Edge* edge2);

void delete_vertex_by_detector(std::vector<qrc::DecodingGraph::Vertex*>& vertices, uint detector);

void addUniqueEdge(std::vector<qrc::DecodingGraph::Edge*>& edges, qrc::DecodingGraph::Edge* newEdge);

std::vector<std::vector<bool>> create_combinations(uint m, uint n);

CountResult compare_matchings(std::map<uint, uint>& A, std::map<uint, uint>& B, bool print=true);
CountResult accuracy_coverage(std::map<uint, uint>& A, std::map<uint, uint>& B);
CountResult accuracy_coverage(std::map<uint, uint>& A, std::vector<std::pair<uint, uint>>& B);

void custom_sort(std::vector<qrc::DecodingGraph::Edge*>& edges, const std::map<uint, uint>& degree_map,
qrc::DecodingGraph decoding_graph,std::map<qrc::DecodingGraph::Vertex*, std::vector<qrc::DecodingGraph::Vertex*>> adjacency_matrix);

std::vector<std::map<uint, uint>> getAllCompleteMatchings(std::vector<uint16_t>& numbers);

void generateCompleteMatchings(std::vector<uint16_t>& numbers, std::vector<std::map<uint, uint>>& result, std::map<uint, uint>& currentMatching, uint index);
    
} // namespace qpd





#endif


 
