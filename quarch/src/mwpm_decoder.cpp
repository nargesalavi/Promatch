/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#include "mwpm_decoder.h"

namespace qrc {

//#define TRY_FILTER
#ifdef TRY_FILTER
static bool did_get_search_space_size = false;
static std::set<std::map<uint, uint>> search_set;
uint64_t _gsss_helper(
        const std::map<uint, std::vector<uint>>& mate_list,
        const std::map<uint, uint>& curr_matching)
{
    if (curr_matching.size() == mate_list.size()) {
        if (search_set.count(curr_matching)) {
            return 0;
        }
        search_set.insert(curr_matching);
        return 1;
    }
    uint64_t n_matchings = 0;
    for (auto pair : mate_list) {
        uint vi = pair.first;
        auto vj_list = pair.second;
        if (curr_matching.count(vi)) {
            continue;
        }
        for (uint vj : vj_list) {
            if (curr_matching.count(vj) || vj <= vi) {
                continue;
            }
            std::map<uint, uint> new_matching(curr_matching);
            new_matching[vi] = vj;
            new_matching[vj] = vi;
            n_matchings += _gsss_helper(mate_list, new_matching);
        }
    }
    return n_matchings;
}

uint64_t get_search_space_size(const std::map<uint, std::vector<uint>>& mate_list) {
    std::map<uint, uint> init;
    return _gsss_helper(mate_list, init); 
}
#endif

MWPMDecoder::MWPMDecoder(const stim::Circuit& circ, uint max_detector) 
    :Decoder(circ),
    longest_error_chain(0),
    path_table(),
    max_detector(max_detector)
{
    path_table = compute_path_table(graph);
}

std::string
MWPMDecoder::name() {
    return std::string(MWPM_DECODER_NAME);
}

bool
MWPMDecoder::is_software() {
    return true;
}

uint64_t
MWPMDecoder::sram_cost() {
    // We examine the size of the path_table.
    uint64_t n_bytes_sram = 0;   
    for (auto kv_pair : path_table) {
        auto res = kv_pair.second;
        // We assume that the entry is stored
        // in a hardware matrix. We don't count the 
        // SRAM required to store keys.
        // 
        // Instead, we will only consider the SRAM
        // required to store the entry.
        uint64_t bytes_for_path = sizeof(res.path[0]) * res.path.size();
        uint64_t bytes_for_distance = sizeof(res.distance);
        n_bytes_sram += bytes_for_path + bytes_for_distance;
    }
    return n_bytes_sram;
}

DecoderShotResult
MWPMDecoder::decode_error(const std::vector<uint8_t>& syndrome) {
    // Build Boost graph for MWPM.
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();

    // Note to self: fault ids in pymatching are the frames in DecodingGraph.
    // Log start time.
#ifdef __APPLE__
    auto start_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
    struct timespec start_time_data;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time_data);
    auto start_time = start_time_data.tv_nsec;
#endif
    // Count number of detectors.
    uint8_t syndrome_is_even = 0x1;
    std::vector<uint> detector_list;
    for (uint di = 0; di < n_detectors; di++) {
        auto syndrome_bit = syndrome[di];
        if (di > max_detector) {
            syndrome_bit = 0;
        }
        if (syndrome_bit) {
            syndrome_is_even ^= 0x1;
            detector_list.push_back(di);
        }
    }

    if (!syndrome_is_even) {
        // Add boundary to matching graph.
        detector_list.push_back(BOUNDARY_INDEX);
    }
    // Build Blossom V instance.
    uint n_vertices = detector_list.size();
    uint n_edges = n_vertices * (n_vertices + 1) / 2;  // Graph is complete.
    PerfectMatching pm(n_vertices, n_edges);
    pm.options.verbose = false;
#ifdef TRY_FILTER
    PerfectMatching pmfilt(n_vertices, n_edges);
    pmfilt.options.verbose = false;
    uint32_t n_filtedges = 0;
    std::map<uint, std::vector<uint>> possible_mates;
#endif
    // Add edges.
    for (uint vi = 0; vi < n_vertices; vi++) {
        uint di = detector_list[vi];
        DecodingGraph::Vertex * vdi = graph.get_vertex(di);
        for (uint vj = vi + 1; vj < n_vertices; vj++) {
            uint dj = detector_list[vj];
            DecodingGraph::Vertex * vdj = graph.get_vertex(dj);
            auto vdi_vdj = std::make_pair(vdi, vdj);
            if (path_table[vdi_vdj].distance >= 1000.0) 
            {
                continue;  // There is no path.
            }
            fp_t raw_weight = path_table[vdi_vdj].distance;
            // Note that we have already typedef'd qfp_t as wgt_t.
            wgt_t edge_weight = (wgt_t) (MWPM_INTEGER_SCALE * raw_weight);
            pm.AddEdge(vi, vj, edge_weight);
#ifdef TRY_FILTER
            if (edge_weight <= 9000) {
                pmfilt.AddEdge(vi, vj, edge_weight); 
                n_filtedges++;
                if (!possible_mates.count(di)) {
                    possible_mates[di] = std::vector<uint>();
                }
                possible_mates[di].push_back(dj);
                if (!possible_mates.count(dj)) {
                    possible_mates[dj] = std::vector<uint>();
                }
                possible_mates[dj].push_back(di);

            }
#endif
        }
    }
    // Solve instance.
    pm.Solve();
#ifdef TRY_FILTER
    fp_t weight = 0.0;
#endif
    std::map<uint, uint> matching;
    for (uint vi = 0; vi < n_vertices; vi++) {
        uint vj = pm.GetMatch(vi);
        uint di = detector_list[vi];
        uint dj = detector_list[vj];
        // Update matching data structure.
        matching[di] = dj;
        matching[dj] = di;
#ifdef TRY_FILTER
        weight += path_table[std::make_pair(di, dj)].distance;
#endif
    }
    fp_t pm_weight = 0;
    // Compute logical correction.
    std::vector<uint8_t> correction = get_correction_from_matching(matching);
    bool is_error = 
        is_logical_error(correction, syndrome, n_detectors, n_observables);
#ifdef TRY_FILTER
    uint hw = std::accumulate(syndrome.begin(),
                                syndrome.begin()+n_detectors,
                                0);
    if (hw > 10) {
        pmfilt.Solve();
        std::map<uint, uint> filtmatching;
        std::cout << "Syndrome (HW = " << hw << "): ";
        for (uint i = 0; i < n_detectors; i++) {
            if (syndrome[i]) {
                std::cout << possible_mates[i].size();
            } else {
                std::cout << ".";
            }
        }
        std::cout << "\n";
        fp_t filt_weight = 0.0;
        for (uint vi = 0; vi < n_vertices; vi++) {
            uint vj = pmfilt.GetMatch(vi);
            uint di = detector_list[vi];
            uint dj = detector_list[vj];
            // Update matching data structure.
            filtmatching[di] = dj;
            filtmatching[dj] = di;
            fp_t w = path_table[std::make_pair(di, dj)].distance;
            filt_weight += w;
            std::cout << "\tUsed weight: " << w << "\n";
        }
        auto filt_correction = get_correction_from_matching(filtmatching);
        bool filt_is_error =
            is_logical_error(filt_correction, syndrome, n_detectors, n_observables);
        std::cout << "\tOriginal error: " << is_error 
            << ", Filtered error: " << filt_is_error << "\n";
        std::cout << "\tFiltered edges: " << n_filtedges << " of " << n_edges << "\n";
        std::cout << "\tWeights: " << weight << ", " << filt_weight << "\n";
        std::cout << "\tWeights are equal: " << (weight == filt_weight) << "\n";
        if (!did_get_search_space_size && hw >= 16) {
            did_get_search_space_size = true;
            uint64_t n_matchings = get_search_space_size(possible_mates);
            for (auto pair : possible_mates) {
                std::cout << pair.first << ":";
                uint di = detector_list[pair.first];
                for (uint vj : pair.second) {
                    uint dj = detector_list[vj];
                    fp_t w = path_table[std::make_pair(di, dj)].distance;
                    std::cout << " " << vj << "(";
                    if (w <= 7) {
                        std::cout << "L";
                    } else {
                        std::cout << "M";
                    }
                    std::cout << ")";
                }
                std::cout << "\n";
            }
            std::cout << "\tSEARCH SPACE SIZE = " << n_matchings << "\n";
        }
    }
#endif
    // Stop time here.
#ifdef __APPLE__
    auto end_time = clock_gettime_nsec_np(CLOCK_MONOTONIC_RAW);
#else
    struct timespec end_time_data;
    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time_data);
    auto end_time = end_time_data.tv_nsec;
#endif
    auto time_taken = end_time - start_time;
    // Build result.
   DecoderShotResult res = {
        time_taken,
        0.0, // TODO
        is_error,
        correction,
        matching,
        pm_weight
    };
    return res;
}

std::vector<uint8_t>
MWPMDecoder::get_correction_from_matching(const std::map<uint, uint>& matching) {
    std::set<uint> visited;
    std::vector<uint8_t> correction(circuit.count_observables(), 0);
    for (auto di_dj : matching) {
        uint di = di_dj.first;
        uint dj = di_dj.second;
        if (visited.count(di) || visited.count(dj)) {
            continue;
        }
        // Check path between the two detectors.
        // This is examining the error chain.
        DecodingGraph::Vertex * vdi = graph.get_vertex(di);
        DecodingGraph::Vertex * vdj = graph.get_vertex(dj);
        auto vdi_vdj = std::make_pair(vdi, vdj);
        std::vector<DecodingGraph::Vertex*> detector_path(path_table[vdi_vdj].path);
        for (uint i = 1; i < detector_path.size(); i++) {
            // Get edge from decoding graph.
            auto wi = detector_path[i-1];
            auto wj = detector_path[i];
            auto edge = graph.get_edge(wi, wj);
            // The edge should exist.
            for (uint obs : edge->frames) {
                // Flip the bit.
                if (obs >= 0) {
                    correction[obs] = !correction[obs];
                }
            }
        }
        if (detector_path.size()-1 > longest_error_chain) {
            longest_error_chain = detector_path.size() - 1;
        }
        visited.insert(di);
        visited.insert(dj);
    }
    return correction;
}

}  // qrc
