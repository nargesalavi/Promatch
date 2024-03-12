/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#include "decoding_graph.h"

namespace qrc {

static int32_t V_INDEX = 0;
static int32_t E_INDEX = 0;

DecodingGraph::DecodingGraph() 
    :detector_to_vertex(),
    vertices_to_edge(),
    vertex_to_next_round(),
    vertex_to_prev_round(),
    boundary_coord(),
    vertex_list(),
    adjacency_matrix()
{
    boundary_coord.fill((uint)-1);
    add_detector(BOUNDARY_INDEX, boundary_coord);
}

DecodingGraph::~DecodingGraph() {
}

uint
DecodingGraph::count_detectors() {
    return vertex_list.size();
}

void
DecodingGraph::add_detector(uint det, std::array<fp_t, N_COORD>& coord) {
    if (detector_to_vertex.count(det) && detector_to_vertex[det] != nullptr) {
        // Simple update the coord data.
        detector_to_vertex[det]->coord = coord;
    } else {
        Vertex * v = new Vertex(V_INDEX++, coord, det);
        detector_to_vertex[det] = v;
        vertex_list.push_back(v);
        adjacency_matrix[v] = std::vector<Vertex*>();
    }
}

void
DecodingGraph::add_edge(uint det1, uint det2,
        fp_t weight, fp_t e_prob, std::set<uint>& frames) 
{
    Vertex * v1 = detector_to_vertex[det1];
    Vertex * v2 = detector_to_vertex[det2];
    
    Edge * e = new Edge(E_INDEX++, det1, det2, weight, e_prob, frames);
    vertices_to_edge[std::make_pair(v1,v2)] = e;
    vertices_to_edge[std::make_pair(v2,v1)] = e;
    adjacency_matrix[v1].push_back(v2);
    adjacency_matrix[v2].push_back(v1);
}

void
DecodingGraph::remove_vertex(Vertex * v) {
    // Delete v from the vertex list
    for (auto it = vertex_list.begin(); it != vertex_list.end(); ) {
        if (*it == v) {
            it = vertex_list.erase(it); 
        } else {
            it++; 
        }
    }
    // And adjacency lists of its neighbors. Delete any
    // edges containing v.
    for (Vertex * w : adjacency_matrix[v]) {
        auto adj_list = adjacency_matrix[w];
        for (auto it = adj_list.begin(); it != adj_list.end(); ) {
            if (*it == v) {
                it = adj_list.erase(it); 
                vertices_to_edge[std::make_pair(v,w)] = nullptr;
                vertices_to_edge[std::make_pair(w,v)] = nullptr;
            } else {
                it++; 
            }
        } 
    }
    adjacency_matrix[v].clear();
    detector_to_vertex[v->detector] = nullptr;
    delete v;
}

std::vector<DecodingGraph::Edge*> 
DecodingGraph::get_vertex_edges_list(DecodingGraph::Vertex* v){
    std::vector<DecodingGraph::Edge*>  edges_list;
    for(Vertex* v_ : adjacency_matrix[v]){
        Edge* e = get_edge(v,v_);
        if( e == nullptr){
            std::cout << "THERE IS AN ERROR IN THE CODE. This message should not be printed! Adjacent node with null edge pointer"; 
        }
        else{
            edges_list.push_back(e);
        }
    }
    return edges_list;
}

void 
DecodingGraph::set_vertices_features(){
    for(Vertex* v : vertex_list){
        v->feature = 0;
        for(Vertex* v_ : adjacency_matrix[v]){
            Edge* e = get_edge(v,v_);
            if( e == nullptr){
                std::cout << "THERE IS AN ERROR IN THE CODE. This message should not be printed! Adjacent node with null edge pointer"; 
            }
            else{
                v->feature+= e->error_probability;
            }
        }

        v->feature /= adjacency_matrix[v].size();
        v->n_neighbors = adjacency_matrix[v].size();
        feature_group[v->feature]++;
        node_degree_group[v->n_neighbors]++;
        std::pair<fp_t, uint> key(v->feature, v->n_neighbors);

        // Increment the count for the pair
        feature_degree_group[key]++;
    }
}

void
DecodingGraph::remove_edge(Edge * e) {
    Vertex * v = get_vertex(e->detectors.first);
    Vertex * w = get_vertex(e->detectors.second);
    vertices_to_edge[std::make_pair(v,w)] = nullptr;
    vertices_to_edge[std::make_pair(w,v)] = nullptr;
    // Remove v from adjacency list of w and vice versa.
    auto& adj_list_v = adjacency_matrix[v];
    auto& adj_list_w = adjacency_matrix[w];

    for (auto it = adj_list_v.begin(); it != adj_list_v.end(); ) {
        if (*it == w) {
            adj_list_v.erase(it); 
        } else {
            it++; 
        }
    }

    for (auto it = adj_list_w.begin(); it != adj_list_w.end(); ) {
        if (*it == v) {
            adj_list_w.erase(it); 
        } else {
            it++; 
        }
    }

    delete e;
}


DecodingGraph::Vertex*
DecodingGraph::get_vertex(uint det_id) {
    if (!detector_to_vertex.count(det_id) 
        || detector_to_vertex[det_id] == nullptr) 
    {
        // Add it to the graph.
        add_detector(det_id, boundary_coord);
    }
    return detector_to_vertex[det_id];
}

DecodingGraph::Vertex*
DecodingGraph::get_next_round(uint det_id) {
    return get_next_round(get_vertex(det_id));
}

DecodingGraph::Vertex*
DecodingGraph::get_next_round(Vertex * v) {
    if (vertex_to_next_round.count(v)) {
        return vertex_to_next_round[v];
    }
    if (v->detector == BOUNDARY_INDEX) {
        vertex_to_next_round[v] = v;
        return v;
    }

    Vertex * next_round = nullptr;
    for (Vertex * w : vertex_list) {
        if (v->coord[0] == w->coord[0]
            && v->coord[1] == w->coord[1]
            && v->coord[2] + 1 == w->coord[2])
        {
            next_round = w;
            break;
        }
    }
    vertex_to_next_round[v] = next_round;
    return next_round;
}

DecodingGraph::Vertex*
DecodingGraph::get_prev_round(uint det_id) {
    return get_prev_round(get_vertex(det_id));
}

DecodingGraph::Vertex*
DecodingGraph::get_prev_round(Vertex * v) {
    if (vertex_to_prev_round.count(v)) {
        return vertex_to_prev_round[v];
    }
    if (v->detector == BOUNDARY_INDEX) {
        vertex_to_prev_round[v] = v;
        return v;
    }

    Vertex * prev_round = nullptr;
    for (Vertex * w : vertex_list) {
        if (v->coord[0] == w->coord[0]
            && v->coord[1] == w->coord[1]
            && v->coord[2] == w->coord[2] + 1)
        {
            prev_round = w;
            break;
        }
    }
    vertex_to_prev_round[v] = prev_round;
    return prev_round;
}

DecodingGraph::Edge*
DecodingGraph::get_edge(uint det1, uint det2) {
    return get_edge(get_vertex(det1), get_vertex(det2));
}

DecodingGraph::Edge*
DecodingGraph::get_edge(Vertex * v1, Vertex * v2) {
    auto v1_v2 = std::make_pair(v1, v2);
    auto v2_v1 = std::make_pair(v2, v1);
    if (vertices_to_edge.count(v1_v2)) {
        return vertices_to_edge[v1_v2]; 
    } else if (vertices_to_edge.count(v2_v1)) {
        return vertices_to_edge[v2_v1];
    } else {
        return nullptr; 
    }
}

uint32_t
DecodingGraph::get_chain_length(uint det1, uint det2) {
    Vertex * src = get_vertex(det1);
    Vertex * dst = get_vertex(det2);
    
    std::deque<Vertex*> bfs_queue{src};
    std::set<Vertex*> visited;
    std::map<Vertex*, uint32_t> distance;
    distance[src] = 0;

    while (!bfs_queue.empty()) {
        Vertex * v = bfs_queue.front();
        bfs_queue.pop_front();
        if (visited.count(v)) {
            continue;
        }

        for (Vertex * w : adjacency_list(v)) {
            if (!distance.count(w) || distance[v] + 1 < distance[w]) {
                distance[w] = distance[v] + 1;
            }
            bfs_queue.push_back(w);
        }
        visited.insert(v);
    }
    return distance[dst];
}

std::vector<DecodingGraph::Vertex*>
DecodingGraph::vertices() {
    return vertex_list;
}

std::vector<DecodingGraph::Vertex*>
DecodingGraph::adjacency_list(Vertex * v) {
    return adjacency_matrix[v];
}

DecodingGraph
to_decoding_graph(const stim::Circuit& qec_circ) {
    DecodingGraph graph;

    stim::DetectorErrorModel dem = 
        stim::ErrorAnalyzer::circuit_to_detector_error_model(
            qec_circ,
            true,  // decompose_errors
            true,  // fold loops
            false, // allow gauge detectors
            1.0,   // approx disjoint errors threshold
            false, // ignore decomposition failures
            false
        );
    // Create callbacks.
    error_callback_f err_f = 
        [&graph](fp_t e_prob, std::vector<uint> dets,
                std::set<uint> frames) 
        {
            if (e_prob == 0 || dets.size() == 0 || dets.size() > 2) {
                return;  // Zero error probability -- not an edge.
            }
            
            if (dets.size() == 1) {
                // We are connecting to the boundary here.
                dets.push_back(BOUNDARY_INDEX);
            }
            // Now, we should only have two entries in det.
            uint det1 = dets[0];
            uint det2 = dets[1];
            auto graph_edge = graph.get_edge(det1, det2);
            if (graph_edge != nullptr) {
                // Get old edge data.
                fp_t old_e_prob = graph_edge->error_probability;
                std::set<uint> old_frames = graph_edge->frames;
                if (frames == old_frames) {
                    e_prob = 
                        e_prob * (1-old_e_prob) + old_e_prob * (1-e_prob);
                    // We will introduce a new edge index, so just
                    // delete this.
                    graph.remove_edge(graph_edge);
                }
            }
            fp_t edge_weight = (fp_t)log10((1-e_prob)/e_prob);
            graph.add_edge(det1, det2, edge_weight, e_prob, frames);
        };
    detector_callback_f det_f =
        [&graph](uint det, std::array<fp_t, N_COORD> coords) 
        {
            graph.add_detector(det, coords);
        };
    // Declare coord offset array.
    uint det_offset = 0;
    std::array<fp_t, N_COORD> coord_offset;
    coord_offset.fill(0);  // Zero initialize.
    // Use callbacks to build graph.
    _read_detector_error_model(dem, 1, det_offset, coord_offset,
                                err_f, det_f);
    return graph;
}

graph::PathTable<DecodingGraph::Vertex>
compute_path_table(DecodingGraph& graph) {
    typedef DecodingGraph G;
    typedef DecodingGraph::Vertex V;

    graph::ewf_t<G, V> w = [] (G& g, V * v1, V * v2)
    {
        auto e = g.get_edge(v1, v2);
        if (e == nullptr) {
            return 10000000.0;
        } else {
            return e->edge_weight;
        }
    };

    return graph::compute_path_table(graph, w);
}

void 
_read_detector_error_model(
        const stim::DetectorErrorModel& dem, uint n_iter,
        uint& det_offset, std::array<fp_t, N_COORD>& coord_offset,
        error_callback_f err_f, detector_callback_f det_f) 
{
    while (n_iter--) {  // Need this to handle repeats.
        for (stim::DemInstruction inst : dem.instructions) {
            stim::DemInstructionType type = inst.type;
            if (type == stim::DemInstructionType::DEM_REPEAT_BLOCK) {
                // The targets for this instruction are
                // (1) number of repeats, and
                // (2) block number.
                uint n_repeats = (uint)inst.target_data[0].data;
                stim::DetectorErrorModel subblock = 
                    dem.blocks[inst.target_data[1].data];
                _read_detector_error_model(subblock, n_repeats,
                           det_offset, coord_offset, err_f, det_f);
            } else if (type == stim::DemInstructionType::DEM_ERROR) {
                std::vector<uint> detectors;
                std::set<uint> frames;
                 
                fp_t e_prob = (fp_t)inst.arg_data[0];
                for (stim::DemTarget target : inst.target_data) {
                    if (target.is_relative_detector_id()) {
                        // This is a detector, add it to the list.
                        detectors.push_back(
                                (uint)target.data + det_offset);
                    } else if (target.is_observable_id()) {
                        frames.insert(target.data); 
                    } else if (target.is_separator()) {
                        // This is just due to decomposition.
                        // Handle each part of the decomposition
                        // separately.
                        err_f(e_prob, detectors, frames);
                        // Clear detectors and frames.
                        // We have already done the callback.
                        detectors.clear();
                        frames.clear();
                    }
                }
                // Handle last error.
                err_f(e_prob, detectors, frames);
            } else if (type == stim::DemInstructionType::DEM_SHIFT_DETECTORS) {
                det_offset += inst.target_data[0].data;
                uint k = 0;
                for (double a : inst.arg_data) {
                    coord_offset[k++] += (fp_t)a;
                }
            } else if (type == stim::DemInstructionType::DEM_DETECTOR) {
                // Compute coordinates.
                std::array<fp_t, N_COORD> coords(coord_offset);
                uint k = 0;
                for (double a : inst.arg_data) {
                    coords[k++] += (fp_t)a; 
                }
                // Now go through all declared detectors.
                for (stim::DemTarget target : inst.target_data) {
                    det_f(target.data + det_offset, coords);
                }
            }
        }
    }
}


}  // qrc
