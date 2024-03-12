/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#include "astrea.h"

namespace qrc {

Astrea::Astrea(const stim::Circuit circuit,
        uint n_detectors_per_round,
        uint32_t weight_filter_cutoff,
        const AstreaParams& params,
        fp_t time_limit)
    :astrea::MLDDecoder(circuit), 
    // Statistics
    n_nonzero_syndromes(0),
    n_hhw_syndromes(0), total_bfu_cycles(0),
    total_prefetch_cycles(0),
    total_cycles_to_converge(0),
    total_logfilter_savings(0),
    max_bfu_cycles(0),
    max_prefetch_cycles(0),
    max_cycles_to_converge(0),
    min_filter_savings(std::numeric_limits<fp_t>::max()),
    max_hamming_weight(0),
    // Memory system
    simulator(nullptr),
    // Properties
    n_rounds(circuit.count_detectors()/n_detectors_per_round),
    main_clock_frequency(params.main_clock_frequency),
    use_mld(params.use_mld),
    time_limit(time_limit),
    baseline(circuit)
{
    // Initialize Memory System.
    // Define memory event table and associated callbacks.

    astrea::AstreaSimulatorParams sim_params = {
        circuit.count_detectors()+1,
        n_detectors_per_round,
        params.n_registers,
        params.bfu_fetch_width,
        params.bfu_compute_stages,
        params.bfu_priority_queue_size,
        weight_filter_cutoff,
        params.use_mld
    };
    simulator = new astrea::AstreaSimulator(graph, circuit, sim_params);
}

Astrea::~Astrea() {
    delete simulator;
}

std::string
Astrea::name() {
    return "Astrea";
}

bool
Astrea::is_software() {
    return false;
}

uint64_t
Astrea::sram_cost() {
    return 0;   // TODO
}

DecoderShotResult
Astrea::decode_error(const std::vector<uint8_t>& syndrome) {
    uint n_detectors = circuit.count_detectors();
    uint n_observables = circuit.count_observables();

    // Compute Hamming weight.
    // Don't count the observables.
    uint hw = std::accumulate(syndrome.begin(), syndrome.end()-n_observables, 0);
    if (hw > max_hamming_weight) {
        max_hamming_weight = hw;
    }
    if (hw > 0) {
        n_nonzero_syndromes++;
        if (hw > 10) {
            n_hhw_syndromes++;
        }
    }
    if (hw <= 10) {
        uint64_t cycles;
        if (hw & 0x1) {
            hw++;
        }
        if (hw <= 2) {
            cycles = 0;
        } else if (hw <= 6) {
            cycles = 3*hw*(hw-1) + 1;
        } else if (hw == 8) {
            cycles = 3*hw*(hw-1) + 15;
        } else if (hw == 10) {
            cycles = 3*hw*(hw-1) + 73;
        }
        fp_t time_taken = cycles / main_clock_frequency * 1e9;
        DecoderShotResult res;
        res = MWPMDecoder::decode_error(syndrome);
        res.execution_time = time_taken;
        return res;
    }
    // Invoke BFUs.
    // We use a simulator to count the number of cycles.
    std::vector<uint> detector_array;
    for (uint di = 0; di < n_detectors; di++) {
        if (syndrome[di]) {
            detector_array.push_back(di);
        }
    }
    if (detector_array.size() & 0x1) {
        // Add boundary.
        detector_array.push_back(BOUNDARY_INDEX);
    }
    // Run simulator.
#ifdef ASTREA_DEBUG
    std::cout << "==============================================\n";
#endif
    simulator->reset_stats();
    simulator->load_detectors(detector_array);

    uint64_t n_cycles = 0;
    uint round = 1;

    uint32_t total_cycles = 0;
    while (!simulator->is_idle()) {
        fp_t t = n_cycles / main_clock_frequency * 1e9;
        if (t >= time_limit) {
            if (round < n_rounds) {
                simulator->sig_end_round();
                n_cycles = 0;
                round++;
            } else {
                break;
            }
        }
        simulator->tick();
        n_cycles++;
        total_cycles++;
    }
    // Get matching from simulator.
    auto matching = simulator->get_matching();
    std::vector<uint8_t> correction = get_correction_from_matching(matching);

    // Update stats.
    total_bfu_cycles += simulator->bfu_cycles;
    total_prefetch_cycles += simulator->prefetch_cycles;
    total_cycles_to_converge += simulator->cycles_to_converge;
    if (simulator->bfu_cycles > max_bfu_cycles) {
        max_bfu_cycles = simulator->bfu_cycles;
    }
    if (simulator->prefetch_cycles > max_prefetch_cycles) {
        max_prefetch_cycles = simulator->prefetch_cycles;
    }
    if (simulator->cycles_to_converge > max_cycles_to_converge) {
        max_cycles_to_converge = simulator->cycles_to_converge;
    }
    if (hw > 10) {
        fp_t filter_savings;
        if (hw & 0x1) {
            fp_t max_pairs = hw * (hw+1) * 0.5;
            filter_savings = ((fp_t)simulator->valid_weights_after_filter) 
                                / max_pairs;
        } else {
            fp_t max_pairs = hw * (hw-1) * 0.5;
            filter_savings = ((fp_t)simulator->valid_weights_after_filter) 
                                / max_pairs;
        }
        total_logfilter_savings += log(filter_savings);
        if (filter_savings < min_filter_savings) {
            min_filter_savings = filter_savings;
        }
    }

    fp_t time_taken = n_cycles / main_clock_frequency * 1e9;

    bool is_error = 
        is_logical_error(correction, syndrome, n_detectors, n_observables);
#ifdef ASTREA_DEBUG
        std::cout << "Is error: " << is_error << "\n";
        std::cout << "Time taken: " << time_taken << "ns.\n";
        std::cout << "Prefetch cycles: " << simulator->prefetch_cycles << "\n";
        std::cout << "BFU Cycles: " << simulator->bfu_cycles << "\n";
        if (is_error) {
            auto mwpm_res = baseline.decode_error(syndrome);
            auto mwpm_matching = mwpm_res.matching;
            std::cout << "MWPM Matching:\n";
            fp_t mwpm_weight = 0.0;
            for (auto pair : mwpm_matching) {
                auto vdi = graph.get_vertex(pair.first);
                auto vdj = graph.get_vertex(pair.second);
                fp_t w = path_table[std::make_pair(vdi, vdj)].distance;
                std::cout << "\t" << pair.first << " --> " << pair.second 
                    << " ( w = " << w << " )\n";
                mwpm_weight += w;
            }
            fp_t weight = 0.0;
            std::cout << "Astrea Matching:\n";
            for (auto pair : matching) {
                auto vdi = graph.get_vertex(pair.first);
                auto vdj = graph.get_vertex(pair.second);
                fp_t w = path_table[std::make_pair(vdi, vdj)].distance;
                std::cout << "\t" << pair.first << " --> " << pair.second 
                    <<  " ( w = " << w << " )\n";
                weight += w;
            }
            std::cout << "mwpm = " << mwpm_weight 
                << " , astrea = " << weight << "\n";
        }
#endif
    DecoderShotResult res = {
        time_taken,
        0.0, // TODO
        is_error,
        correction,
        matching
    };
    return res;
}

}  // namespace qrc
