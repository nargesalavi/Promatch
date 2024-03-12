/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DECODING_GRAPH_h
#define DECODING_GRAPH_h

#include "defs.h"
#include "graph/dijkstra.h"

#include <stim.h>

#include <array>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

#include <math.h>

#define N_COORD 100

namespace qrc {

#define BOUNDARY_INDEX ((uint)-1)

class DecodingGraph {
public:
    DecodingGraph();
    ~DecodingGraph();

    struct Vertex {
        int32_t id;
        std::array<fp_t, N_COORD> coord;
        uint detector;
        fp_t feature = -1;
        uint n_neighbors = -1;
        uint n_dependents = 0;
        uint n_cluster_members = 1;

        Vertex()
            :id(-1), coord(), detector(0)
        {}

        Vertex(uint id, std::array<fp_t, N_COORD> coord, uint detector)
            :id(id), coord(coord), detector(detector)
        {}

        Vertex(const Vertex& other)
            :id(other.id), coord(other.coord), detector(other.detector), feature(other.feature),
            n_neighbors(other.n_neighbors), n_dependents(other.n_dependents), n_cluster_members(other.n_cluster_members)
        {}

        bool operator==(const Vertex& other) const {
            return id == other.id; 
        }

        bool operator <(const Vertex& other) const {
            return id < other.id; 
        }

        bool operator!=(const Vertex& other) const {
            return !(*this == other); 
        }
    };

    struct Edge {
        int32_t id;
        std::pair<uint, uint> detectors;
        fp_t edge_weight;
        fp_t error_probability;
        std::set<uint> frames;
        fp_t probability_impact;

        Edge()
            :id(-1), detectors(std::make_pair(0,0)), edge_weight(0),
            error_probability(0), frames(std::set<uint>())
        {}

        Edge(int32_t id, uint di, uint dj, fp_t w, fp_t p, std::set<uint> frames)
            :id(id), detectors(std::make_pair(di, dj)), edge_weight(w),
            error_probability(p), frames(frames)
        {}

        Edge(const Edge& other)
            :id(other.id), detectors(other.detectors), 
            edge_weight(other.edge_weight), 
            error_probability(other.error_probability),
            frames(other.frames)
        {}

        bool operator==(const Edge& other) const {
            return id == other.id; 
        }

        bool operator <(const Edge& other) const {
            return id < other.id; 
        }

        bool operator!=(const Edge& other) const {
            return !(*this == other); 
        }
    };

    uint count_detectors(void);

    void add_detector(uint, std::array<fp_t, N_COORD>& coord);
    void add_edge(uint det1, uint det2, fp_t weight,
            fp_t e_prob, std::set<uint>& frames);

    void remove_vertex(Vertex*);
    void remove_edge(Edge*);

    Vertex * get_vertex(uint det_id);
    Vertex * get_next_round(uint det_id);
    Vertex * get_next_round(Vertex*);
    Vertex * get_prev_round(uint det_id);
    Vertex * get_prev_round(Vertex*);
    Edge * get_edge(uint, uint);
    Edge * get_edge(Vertex*, Vertex*);

    uint32_t get_chain_length(uint det1, uint det2);
   
    std::vector<Vertex*> vertices(void);
    std::vector<Vertex*> adjacency_list(Vertex*);
    std::vector<Edge*> get_vertex_edges_list(Vertex*);

    void set_vertices_features();

    //Gathering some statistics about the graph
    std::map<fp_t, uint64_t> feature_group;
    std::map<uint, uint64_t> node_degree_group;
    std::map<std::pair<fp_t, uint>, uint64_t>  feature_degree_group;


private:
    std::map<uint, Vertex*> detector_to_vertex;
    std::map<std::pair<Vertex*, Vertex*>, Edge*> vertices_to_edge;
    std::map<Vertex*, Vertex*> vertex_to_next_round;
    std::map<Vertex*, Vertex*> vertex_to_prev_round;
    std::array<fp_t, N_COORD> boundary_coord;

    std::vector<Vertex*> vertex_list;
    std::map<Vertex*, std::vector<Vertex*>> adjacency_matrix;
};

DecodingGraph
to_decoding_graph(const stim::Circuit&);

graph::PathTable<DecodingGraph::Vertex>
compute_path_table(DecodingGraph&);

std::vector<DecodingGraph::Edge*> vertex_edges_list(DecodingGraph::Vertex*);



typedef std::function<void(fp_t, std::vector<uint>, std::set<uint>)>
    error_callback_f;
typedef std::function<void(uint, std::array<fp_t, N_COORD>)>
    detector_callback_f;

void
_read_detector_error_model(const stim::DetectorErrorModel&, 
        uint n_iter, uint& det_offset, 
        std::array<fp_t, N_COORD>& coord_offset,
        error_callback_f, detector_callback_f);

}  // qrc

#endif
