/*
 *  author: Narges Alavisamani
 *  date:   6 March 2023
 * */

#include "predecoding_graph_helper.h"

namespace qpd{

void dfs_(qrc::DecodingGraph::Vertex* v, 
        std::unordered_set<qrc::DecodingGraph::Vertex*>& visited, 
        std::vector<qrc::DecodingGraph::Vertex*>& component,
        std::map<qrc::DecodingGraph::Vertex*, std::vector<qrc::DecodingGraph::Vertex*>>& adjacency_matrix) {
    visited.insert(v);
    component.push_back(v);
    for (auto neighbor : adjacency_matrix[v]) {
        if (visited.find(neighbor) == visited.end()) {
            dfs_(neighbor, visited, component, adjacency_matrix);
        }
    }
}

std::vector<std::vector<qrc::DecodingGraph::Vertex*>> 
find_connected_components(std::map<qrc::DecodingGraph::Vertex*, std::vector<qrc::DecodingGraph::Vertex*>>& adjacency_matrix) {
    std::unordered_set<qrc::DecodingGraph::Vertex*> visited;
    std::vector<std::vector<qrc::DecodingGraph::Vertex*>> components;
    
    for (auto const& entry : adjacency_matrix) {
        auto v = entry.first;
        if (visited.find(v) == visited.end()) {
            std::vector<qrc::DecodingGraph::Vertex*> component;
            dfs_(v, visited, component, adjacency_matrix);
            components.push_back(component);
        }
    }
    return components;
}

void print_components(std::map<qrc::DecodingGraph::Vertex*, 
                        std::vector<qrc::DecodingGraph::Vertex*>>& adjacency_matrix, qrc::DecodingGraph& decoding_graph){
        
        auto components = find_connected_components(adjacency_matrix);
        uint i = 0;
        for(const auto& c : components){
            i++;
            std::cout << "_____Component " << i << " _________" << std::endl;
            for(const auto& v : c){
                std::cout << "Vertex: " << v->detector << "(" << v->feature << ") Edges: ";
                fp_t sum = 0;
                auto edges_list = decoding_graph.get_vertex_edges_list(v);
                for(const auto& e: edges_list){
                    std::cout << "(" << e->id << ", " << e->error_probability <<"), ";
                    sum += e->error_probability;
                }
                std::cout << sum/edges_list.size() << " - " << edges_list.size() << std::endl;
            }
            std::cout << "_____________________________________" << std::endl;


        }
        
}

bool isDuplicate(const std::vector<qrc::DecodingGraph::Vertex*>& vertices, uint detector) {
    // Check if a Vertex* with the same detector already exists in the vector
    return std::any_of(vertices.begin(), vertices.end(), [detector](const qrc::DecodingGraph::Vertex* v) {
        return v->detector == detector;
    });
}

#include <tuple>

class Vertex {
public:
    uint detector;
    // Additional attributes and methods...
};

CountResult countMembersInAB(const std::vector<qrc::DecodingGraph::Vertex*>& a, 
                            const std::vector<qrc::DecodingGraph::Vertex*>& b) {
    CountResult result = {0, 0, 0};

    for (qrc::DecodingGraph::Vertex* vertexA : a) {
        uint detectorA = vertexA->detector;

        auto itB = std::find_if(b.begin(), b.end(), [&](qrc::DecodingGraph::Vertex* vertexB) {
            return vertexB->detector == detectorA;
        });

        if (itB != b.end()) {
            result.countBothAB++;
        }
        else{
            result.countANotInB++;
        }
    }

    result.countBNotInA = ((int)b.size()) - result.countBothAB;

    return result;
}

void addUniqueVertex(std::vector<qrc::DecodingGraph::Vertex*>& vertices, qrc::DecodingGraph::Vertex* newVertex) {
    if (!isDuplicate(vertices, newVertex->detector)) {
        vertices.push_back(newVertex);
    }
}

bool compare_edges_by_error_probability(const qrc::DecodingGraph::Edge* edge1, const qrc::DecodingGraph::Edge* edge2) {
    return edge1->error_probability > edge2->error_probability;
}
bool compare_edges_by_future_error_probability(const qrc::DecodingGraph::Edge* edge1, const qrc::DecodingGraph::Edge* edge2) {
    return edge1->probability_impact > edge2->probability_impact;
}
void custom_sort(std::vector<qrc::DecodingGraph::Edge*>& edges, const std::map<uint, uint>& degree_map,
qrc::DecodingGraph decoding_graph, std::map<qrc::DecodingGraph::Vertex*, std::vector<qrc::DecodingGraph::Vertex*>> adjacency_matrix) {
    int n = edges.size();
    bool swapped;
    for (int i = 0; i < n - 1; i++) {
        swapped = false;
        for (int j = 0; j < n - i - 1; j++) {
            const qrc::DecodingGraph::Edge* e1 = edges[j];
            const qrc::DecodingGraph::Edge* e2 = edges[j + 1];            

            bool e1_degree_1 = (degree_map.at(e1->detectors.first) == 1 && degree_map.at(e1->detectors.second) == 1);
            bool e2_degree_1 = (degree_map.at(e2->detectors.first) == 1 && degree_map.at(e2->detectors.second) == 1);

            if (!e1_degree_1 && e2_degree_1) {
                std::swap(edges[j], edges[j + 1]);
                swapped = true;
            }
            else if ((!e1_degree_1 && !e2_degree_1) || (e1_degree_1 && e2_degree_1) ){
            
            
                int e1_min_degree = std::min(degree_map.at(e1->detectors.first), degree_map.at(e1->detectors.second));
                int e2_min_degree = std::min(degree_map.at(e2->detectors.first), degree_map.at(e2->detectors.second));
                
                if (e1_min_degree > e2_min_degree || (e1_min_degree == e2_min_degree && e1->error_probability < e2->error_probability)) {
                    std::swap(edges[j], edges[j + 1]);
                    swapped = true;
                }
            }
        }
        if (!swapped) {
            break;
        }
    }
}

void delete_vertex_by_detector(std::vector<qrc::DecodingGraph::Vertex*>& vertices, uint detector) {
    // Use std::remove_if with a lambda function as the predicate
    auto newEnd = std::remove_if(vertices.begin(), vertices.end(), [detector](const qrc::DecodingGraph::Vertex* vertex) {
        return vertex->detector == detector;
    });

    // Erase the matching elements from the vector
    vertices.erase(newEnd, vertices.end());
}

CountResult countMembersInAB_edges(const std::vector<qrc::DecodingGraph::Edge*>& a,
                             const std::vector<qrc::DecodingGraph::Edge*>& b) {
    CountResult result;
    result.countANotInB = 0;
    result.countBNotInA = 0;
    result.countBothAB = 0;

    for (const auto* edgeA : a) {
        auto itB = std::find_if(b.begin(), b.end(), [&](const auto* edgeB) {
            return edgeB->id == edgeA->id;
        });
        if (itB != b.end()) {
            result.countBothAB++;
        } else {
            result.countANotInB++;
        }
    }

    result.countBNotInA = ((int)b.size()) - result.countBothAB;

    return result;
}

void addUniqueEdge(std::vector<qrc::DecodingGraph::Edge*>& edges, qrc::DecodingGraph::Edge* newEdge) {
    auto it = std::find_if(edges.begin(), edges.end(), [&](const qrc::DecodingGraph::Edge* edge) {
        return edge->id == newEdge->id;
    });

    if (it == edges.end()) {
        edges.push_back(newEdge);
    }
}

std::vector<std::vector<bool>> create_combinations(uint m, uint n) {

    std::vector<std::vector<bool>> chooses;
    if(m == 0){
        return chooses;
    }
    if(m < n){
        n = m;
    }
    std::vector<bool> bits(m, false);   // Create a vector of size m, initialized with false
    for (uint i = 0; i < n; ++i) {
        bits[i] = true;  // Set the first n elements to true
    }


    do {
        chooses.push_back(bits);
    } while (std::next_permutation(bits.begin(), bits.end()));
    std::reverse(bits.begin(), bits.end());  // Reverse the vector

    while (std::prev_permutation(bits.begin(), bits.end())){
        chooses.push_back(bits);
    }

    return chooses; 
}

CountResult compare_matchings(std::map<uint, uint>& A, std::map<uint, uint>& B, bool print){
    std::map<uint,uint> printed_pairs;
    CountResult res = {0,0,0};

    for(auto a_map: A){
        uint first = a_map.first;
        uint second = a_map.second;
        if(printed_pairs.count(first)==0 && printed_pairs.count(second)==0 ){
            printed_pairs[first] = second;
            if(B[first] != second || B[second] != first ){
                res.countANotInB++;
                res.countBNotInA++;
                if(print){
                    std::cout << "(" << first << ", " << second << ")-("<< first << ", " << B[first] << "), ";   
                }
            }
            else{
                res.countBothAB++;
            }

        }
    }
    std::cout << std::endl;
    return res;
}

CountResult accuracy_coverage(std::map<uint, uint>& A, std::map<uint, uint>& B){
    CountResult res = {0,0,0};

    for(auto a_map: A){
        uint first = a_map.first;
        uint second = a_map.second;
        {
            if(B.count(first) == 0 && B.count(second) == 0){
                
                res.countANotInB++; //skipped
            }
        }
    }
    for(auto b_map: B){
        uint first = b_map.first;
        uint second = b_map.second;
        {
            if(A[first] != second || A[second] != first ){
                
                res.countBNotInA++; // incorrect
            }
            else{
                res.countBothAB++; //  accuracy
            }

        }
    }
    // if(res.countBNotInA !=0 & res.countANotInB!=0){
    //     res.countANotInB+=2;
    // }
    return res;
}

CountResult accuracy_coverage(std::map<uint, uint>& A, std::vector<std::pair<uint, uint>>& B){
    CountResult res = {0,0,0};

    for (const auto& a_map : A) {
        uint32_t first = a_map.first;
        uint32_t second = a_map.second;

        bool foundInB = false;
        for (const auto& b_pair : B) {
            if ((first == b_pair.first || second == b_pair.second || first == b_pair.second || second == b_pair.first)) {
                foundInB = true;
                break;
            }
        }

        if (!foundInB) {
            res.countANotInB++; // Skipped
        }
    }

    for (const auto& b_pair : B) {
        uint32_t first = b_pair.first;
        uint32_t second = b_pair.second;

        if ((A[first] != second || A[second] != first)) {
            res.countBNotInA++; // Incorrect
        }
        if( (A[first] == second)){
            res.countBothAB++; // Accuracy
        }
    }
    res.countANotInB /= 2;
    res.countBNotInA /= 2;
    res.countBothAB /= 2;
    

    return res;
}

std::vector<qrc::DecodingGraph::Edge*> prioritizeEdgesByDegree(qrc::DecodingGraph* graph_, const std::map<qrc::DecodingGraph::Vertex*, 
std::vector<qrc::DecodingGraph::Vertex*>>& adjacency_matrix, 
const std::vector<qrc::DecodingGraph::Edge*>& sorted_predecoding_edges) {
    std::vector<qrc::DecodingGraph::Edge*> degree_priority;
    
    for (qrc::DecodingGraph::Edge* edge : sorted_predecoding_edges) {
        qrc::DecodingGraph::Vertex* vertex1 = graph_->get_vertex(edge->detectors.first);
        qrc::DecodingGraph::Vertex* vertex2 = graph_->get_vertex(edge->detectors.second);
        
        size_t degree1 = adjacency_matrix.at(vertex1).size();
        size_t degree2 = adjacency_matrix.at(vertex2).size();
        
        if (degree1 < degree2 || (degree1 == degree2 && edge->error_probability > degree_priority[0]->error_probability)) {
            degree_priority.insert(degree_priority.begin(), edge);
        } else {
            degree_priority.push_back(edge);
        }
    }
    
    return degree_priority;
}

void generateCompleteMatchings(std::vector<uint16_t>& numbers, std::vector<std::map<uint, uint>>& result, std::map<uint, uint>& currentMatching, uint index) {
    std::cout << result.size() << std::endl;
    if (index == numbers.size()) {
        result.push_back(currentMatching);
        currentMatching.clear();
        return;
    }

    uint currentNumber = numbers[index];
    for (uint i = index + 1; i < numbers.size(); i++) {
        uint nextNumber = numbers[i];
        if (currentMatching.count(currentNumber) == 0 && currentMatching.count(nextNumber) == 0) {
            currentMatching[currentNumber] = nextNumber;
            currentMatching[nextNumber] = currentNumber;
            
            // currentMatching.erase(currentNumber);
            // currentMatching.erase(nextNumber);
            
            

            generateCompleteMatchings(numbers, result, currentMatching, index);
            
        }
    }
}

std::vector<std::map<uint, uint>> getAllCompleteMatchings(std::vector<uint16_t>& numbers) {
    std::vector<std::map<uint, uint>> result;
    std::map<uint, uint> currentMatching;
    generateCompleteMatchings(numbers, result, currentMatching, 0);
    return result;
}

} // qpd namespace