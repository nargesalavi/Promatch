/*
 *  author: Suhas Vittal
 *  date:   13 August 2022
 * */

#include "decoder.h"

namespace qrc {

Decoder::Decoder(const stim::Circuit& circ)
    :circuit(circ), 
    graph(to_decoding_graph(circ)),
    syndromes(),
    execution_times(),
    memory_overheads(),
    n_logical_errors(0),
    mean_execution_time(0),
    max_execution_time(0),
    max_execution_time_for_correctable(0),
    hamming_weight_dist()
{}

void
Decoder::clear_stats() {
    syndromes.clear();
    execution_times.clear();
    memory_overheads.clear();
    n_logical_errors = 0;
    mean_execution_time = 0;
    max_execution_time = 0;
    max_execution_time_for_correctable = 0;
}

std::vector<Decoder::SyndromeEvent>
Decoder::nonzero_syndromes_and_time_taken() {
    std::vector<Decoder::SyndromeEvent> syndrome_list;
    uint32_t n_syndromes = syndromes.size();

    for (uint32_t i = 0; i < n_syndromes; i++) {
        uint hw = std::accumulate(syndromes[i].begin(), syndromes[i].end(), 0);
        if (hw) {
            auto pair = std::make_pair(syndromes[i], execution_times[i]);
            syndrome_list.push_back(pair);
        }
    }

    return syndrome_list;
}

std::vector<Decoder::SyndromeEvent>
Decoder::nonzero_syndromes_completed_within(fp_t max_time, 
        uint32_t& n_nonzero_syndromes) 
{
    std::vector<Decoder::SyndromeEvent> syndrome_list;
    uint32_t n_syndromes = syndromes.size();

    n_nonzero_syndromes = 0;
#ifdef USE_OMP
#pragma omp parallel for reduction (+: n_nonzero_syndromes)
#endif
    for (uint32_t i = 0; i < n_syndromes; i++) {
        uint hw = std::accumulate(syndromes[i].begin(), syndromes[i].end(), 0);
        n_nonzero_syndromes += hw ? 1 : 0;
#ifdef USE_OMP
#pragma omp critical
#endif
        {
            if (hw && execution_times[i] <= max_time) {
                auto pair = std::make_pair(syndromes[i], execution_times[i]);
                syndrome_list.push_back(pair);
            }
        }
    }

    return syndrome_list;
}

uint
Decoder::max_hamming_wgt_completed_within(fp_t max_time) {
    std::vector<uint8_t> can_complete_all(circuit.count_detectors(), 0x1);
    uint32_t n_syndromes = syndromes.size();
#ifdef USE_OMP
#pragma omp parallel for
#endif
    for (uint32_t i = 0; i < n_syndromes; i++) {
        auto syndrome = syndromes[i];
        // Compute Hamming weight.
        uint hw = std::accumulate(syndrome.begin(), syndrome.end(), 0);
#ifdef USE_OMP
#pragma omp atomic update
#endif
        can_complete_all[hw] &= execution_times[i] <= max_time;
    }
    // Go compute the max possible hamming weight.
    uint max_weight = 0;
    for (uint w = 1; w < can_complete_all.size(); w++) {
        if (can_complete_all[w]) {
            max_weight = w;
        } else {
            break;
        }
    }
    return max_weight;
}

uint64_t
Decoder::sram_cost() {
    return 0;
}

uint32_t
syndrome_to_int(const std::vector<uint8_t>& syndrome, uint size) {
    if (size > 32) {
        std::cout << "syndrome is too large to be converted to 32 bit int.\n";
        std::cout << "size is " << size << ".\n";
        return (uint32_t)-1;
    }
    uint32_t x = 0;
    for (uint i = 0; i < size; i++) {
        x |= (syndrome[i] & 0x1) << i;
    }
    return x;
}

bool
is_logical_error(const std::vector<uint8_t>& correction,
        const std::vector<uint8_t>& syndrome,
        uint n_detectors, uint n_observables) 
{
    bool corr_matches_obs = true;
    for (uint i = 0; i < n_observables; i++) {
        uint8_t corr = correction[i];
        auto actual = syndrome[n_detectors + i];
        corr_matches_obs = corr_matches_obs && (corr == actual);
    }
    return !corr_matches_obs;
}

std::vector<uint8_t> _to_vector(const stim::simd_bits_range_ref& array,
        uint n_detectors, uint n_observables) 
{
    std::vector<uint8_t> syndrome(n_detectors+n_observables);
    for (uint i = 0; i < n_detectors + n_observables; i++) {
        syndrome[i] = array[i];
    }
    return syndrome;
}

}  // qrc
