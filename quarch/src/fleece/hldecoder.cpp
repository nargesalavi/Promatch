/*
 *  author: Suhas Vittal
 *  date:   15 March 2023
 *
 *  A more hardware amenable version of the Heralded Leakage Decoder proposed by
 *  Suchara et al. 
 * */

#include "fleece/hldecoder.h"

namespace qrc {
namespace fleece {

HLDecoder::HLDecoder(const stim::Circuit& circ)
    :MWPMDecoder(circ),
    base_path_table(),
    lattice_graph(to_lattice_graph(circ))
{
    base_path_table = compute_path_table(graph);
}

DecoderShotResult
HLDecoder::decode_error(const std::vector<uint8_t>& M_syndrome, const std::vector<uint8_t>& L_syndrome) {
    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();

    std::vector<uint8_t> adjusted_syndrome(M_syndrome);

    std::map<DecodingGraph::Edge*, fp_t> edge_to_prob_old;

    bool any_leaked = false;
    for (uint i = 0; i < n_detectors; i++) {
        if (L_syndrome[i]) {
            uint32_t severity = lattice_graph.get_leakage_severity_by_detector(i);
            adjusted_syndrome[i] ^= L_syndrome[i];
            if (severity == ((uint32_t)-1)) {
                continue;
            }
            any_leaked = true;
            auto dec_v = graph.get_vertex(i);
            for (auto dec_w : graph.adjacency_list(dec_v)) {
                auto e = graph.get_edge(dec_v, dec_w);
                edge_to_prob_old[e] = e->error_probability;
                fp_t s = 1.0/((fp_t)severity);
                e->error_probability = 0.5*s + e->error_probability - e->error_probability * s;
                e->edge_weight = log10((1-e->error_probability)/e->error_probability);
            }
        }
    }

    if (any_leaked) {
        path_table = compute_path_table(graph);
    } else {
        path_table = base_path_table;
    }
    auto res = MWPMDecoder::decode_error(adjusted_syndrome);
    // Undo graph modifications
    for (auto pair : edge_to_prob_old) {
        auto e = pair.first;
        auto e_prob = pair.second;
        e->error_probability = e_prob;
        e->edge_weight = log10((1.0-e_prob)/e_prob);
    }
    return res;
}

}   // fleece
}   // qrc
