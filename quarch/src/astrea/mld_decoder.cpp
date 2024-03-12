/*
 *  author: Suhas Vittal
 *  date:   30 March 2023
 * */

#include "astrea/mld_decoder.h"

namespace qrc {
namespace astrea {

MLDDecoder::MLDDecoder(const stim::Circuit& circ)
    :MWPMDecoder(circ),
    curr_detector_list(),
    representatives(),
    obs_prob()
{}

std::string
MLDDecoder::name() {
    return "Brute Force MLD Decoder";
}

DecoderShotResult
MLDDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    obs_prob.fill(0);
    curr_detector_list.clear();

    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();
    
    for (uint i = 0; i < n_detectors; i++) {
        if (syndrome[i]) {
            curr_detector_list.push_back(i);
        }
    }

    if (curr_detector_list.size() & 0x1) {
        curr_detector_list.push_back(BOUNDARY_INDEX);
    }

    std::map<uint, uint> matching;
    explore(matching, 0);

    if (obs_prob[0] > obs_prob[1]) {
        matching = representatives[0];
    } else {
        matching = representatives[1];
    }
    std::vector<uint8_t> correction = get_correction_from_matching(matching);
    bool is_error = is_logical_error(correction, syndrome, n_detectors, n_observables);

    DecoderShotResult res = {
        0.0,
        0.0,
        is_error,
        correction,
        matching
    };
    return res;
}

fp_t
MLDDecoder::get_obs_prob(uint obs) {
    return obs_prob[obs];
}

std::map<uint, uint>
MLDDecoder::get_representative(uint obs) {
    return representatives[obs];
}

void
MLDDecoder::explore(const std::map<uint, uint>& base, fp_t weight) {
    if (base.size() == curr_detector_list.size()) {
        fp_t p = pow(10, -weight);
        std::vector<uint8_t> corr = get_correction_from_matching(base);
        if (obs_prob[corr[0]] == 0.0) {
            representatives[corr[0]] = base;
        }
        obs_prob[corr[0]] += p;
        return;
    }
    for (uint i = 0; i < curr_detector_list.size(); i++) {
        uint di = curr_detector_list[i];
        if (base.count(di)) {
            continue;
        }
        for (uint j = i+1; j < curr_detector_list.size(); j++) {
            uint dj = curr_detector_list[j];
            if (base.count(dj)) {
                continue;
            }
            std::map<uint, uint> m(base);
            m[di] = dj;
            m[dj] = di;

            auto vi = graph.get_vertex(di);
            auto vj = graph.get_vertex(dj);
            fp_t nw = weight + path_table[std::make_pair(vi, vj)].distance;

            explore(m, nw);
        }
    }
}

}   // astrea
}   // qrc
