/*
 *  author: Suhas Vittal
 *  date:   9 August 2022
 * */

#include "benchmark.h"

namespace qrc {

/* Helper Functions */
inline fp_t __CHS(fp_t x, fp_t y, fp_t z) {
    return z >= 0 ? z : (y >=0 ? y : x);
}

/* Main Functions */

fp_t
min(const std::vector<fp_t>& data) {
    return *std::min_element(data.begin(), data.end());
}

fp_t
max(const std::vector<fp_t>& data) {
    return *std::max_element(data.begin(), data.end());
}

fp_t
mean(const std::vector<fp_t>& data) {
    fp_t sum = std::accumulate(data.begin(), data.end(), 0);
    return sum / data.size();
}
fp_t
stdev(const std::vector<fp_t>& data) {
    fp_t sum = 0.0;
    fp_t m = mean(data);
    for (uint32_t i = 0; i < data.size(); i++) {
        sum += (data[i] - m)*(data[i] - m);
    }
    return sqrt(sum/data.size());
}

void
b_decoder_ler(Decoder * decoder_p, uint64_t shots, std::mt19937_64& rng,
        bool save_per_shot_data) 
{
    const uint64_t _shots = shots;
    // Clear stats.
    decoder_p->clear_stats();
    // Declare statistics
    uint32_t array_size = save_per_shot_data ? shots : 1;

    uint64_t total_shots = shots;
    std::vector<std::vector<uint8_t>> syndromes(array_size);
    std::vector<fp_t> execution_times(array_size);
    std::vector<fp_t> memory_overheads(array_size);
    uint32_t n_logical_errors = 0;
    fp_t mean_execution_time = 0.0,
         max_execution_time = 0.0,
         max_execution_time_for_correctable = 0.0;

    std::map<uint, uint64_t> hamming_weight_freq;

    uint32_t sn = 0;
    uint32_t bn = 0;
    while (shots > 0) {
        uint32_t shots_this_round = shots > MAX_SHOTS ? MAX_SHOTS : shots;
        // Get samples from Stim.
        stim::simd_bit_table leakage_buffer(1, 1);
        stim::simd_bit_table sample_buffer = 
            stim::detector_samples(decoder_p->circuit, shots_this_round,
                    false, true, rng, true, leakage_buffer);
        // Last part of samples is the actual observable.
        // We are trying to match that.
        sample_buffer = sample_buffer.transposed();
        leakage_buffer = leakage_buffer.transposed();
        uint n_detectors = decoder_p->circuit.count_detectors();
        uint n_observables = decoder_p->circuit.count_observables();

        for (uint32_t i = 0; i < shots_this_round; i++) {
            auto syndrome = 
                _to_vector(sample_buffer[i], n_detectors, n_observables);
            uint hw = 
                std::accumulate(syndrome.begin(), syndrome.end() - n_observables, 0);
            if (hw & 0x1) {
                hw++;
            }
            if (!hamming_weight_freq.count(hw)) {
                hamming_weight_freq[hw] = 0;
            }
            hamming_weight_freq[hw]++;
            // Update some stats before decoding.
            if (save_per_shot_data) {
                syndromes[sn] = syndrome;
            }
            if (hw > 0) {
                DecoderShotResult res = decoder_p->decode_error(syndrome);
                // Update statistics.
//                n_logical_errors += (res.is_logical_error || leakage_buffer[i][n_detectors]);
                n_logical_errors += res.is_logical_error;
                mean_execution_time += res.execution_time / ((fp_t)total_shots);
                if (res.execution_time > max_execution_time) {
                    max_execution_time = res.execution_time;
                    if (!res.is_logical_error) {
                        max_execution_time_for_correctable = res.execution_time;    
                    }
                }
                if (save_per_shot_data) {
                    execution_times[sn] = res.execution_time;
                    memory_overheads[sn] = res.memory_overhead;
                }
            } else {
                if (save_per_shot_data) {
                    execution_times[sn] = 0;
                    memory_overheads[sn] = 0;
                }
            }
            sn++;
        }
        shots -= shots_this_round;
        bn++;
    }
    // Update stats in decoder.
    decoder_p->syndromes = syndromes;
    decoder_p->execution_times = execution_times;
    decoder_p->memory_overheads = memory_overheads;
    decoder_p->n_logical_errors = n_logical_errors;
    decoder_p->mean_execution_time = mean_execution_time;
    decoder_p->max_execution_time = max_execution_time;
    decoder_p->max_execution_time_for_correctable = 
        max_execution_time_for_correctable;
    decoder_p->hamming_weight_dist.clear();
    for (auto pair : hamming_weight_freq) {
        decoder_p->hamming_weight_dist[pair.first] = ((fp_t) pair.second) / ((fp_t) _shots);
    }
}

benchmark::StatisticalResult
b_statistical_ler(Decoder * decoder, uint64_t shots_per_batch, std::mt19937_64& rng, bool use_mpi, uint n_faults, uint distance) {
    int world_rank = 0, world_size = 1;
    if (use_mpi) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }
    if(world_rank == 0)
        std::cout << "Distance = " << distance << " - Shots per k = " << shots_per_batch << " - Decoder: " << decoder->name() << std::endl;


    stim::Circuit circuit(decoder->circuit);
    DecodingGraph decoding_graph = to_decoding_graph(circuit);

    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();
    const uint n_results = n_detectors + n_observables;

    uint64_t local_shots = shots_per_batch / world_size;
    if (world_rank == 0) {
        local_shots += shots_per_batch % world_size;
    }
    stim::simd_bit_table result_table(local_shots, n_results);

#define DELTA(a, b) ((a)-(b))/(a)
#define B_STAT_LER_USE_POLY

#ifdef B_STAT_LER_USE_POLY
    std::array<fp_t, 100'000> log_poly;
    log_poly.fill(1);

    uint eventno = 0;
#else
    fp_t mean_flip_prob = 0.0;
#endif
    std::vector<fp_t> edge_probs;
    std::vector<DecodingGraph::Edge*> edge_list;
    for (uint d = 0; d < n_detectors; d++) {
        auto v = decoding_graph.get_vertex(d);
        for (auto w : decoding_graph.adjacency_list(v)) {
            if(w->detector < d){
                continue;
            }
            auto edge = decoding_graph.get_edge(v, w);
            edge_list.push_back(edge);
            edge_probs.push_back(edge->error_probability);
#ifdef B_STAT_LER_USE_POLY
            fp_t ep = edge->error_probability;
            if (eventno == 0) {
                log_poly[0] = log(1-ep);
                log_poly[1] = log(ep);
            } else {
                std::array<fp_t, 100'000> log_pa, log_pb;
                log_pa.fill(1);
                log_pb.fill(1);
                for (uint i = 0; i <= eventno; i++) {
                    log_pa[i] = log_poly[i] + log(1-ep);
                    log_pb[i+1] = log_poly[i] + log(ep);
                }
                for (uint i = 0; i < log_poly.size(); i++) {
                    if (log_pa[i] == 1 && log_pb[i] == 1) {
                        log_poly[i] = 1;
                    } else if (log_pa[i] == 1) {
                        log_poly[i] = log_pb[i];
                    } else if (log_pb[i] == 1) {
                        log_poly[i] = log_pa[i];
                    } else {
                        log_poly[i] = log(pow(M_E, log_pa[i]) + pow(M_E, log_pb[i]));
                    }
                }
            }
            eventno++;
#else
            mean_flip_prob += 1.0/edge->error_probability;
#endif
        }
    }

    benchmark::StatisticalResult statres;

    // std::cout << "Size of edges: " << edge_probs.size() << std::endl;
    // for(int i = 0; i < 28; i++){
    //     std::cout << "K = " << i << " : "<<  pow(M_E, log_poly[i]) << std::endl;
    // }
    // uint64_t s = 0;
    // uint64_t s_ex = 0;
    // for(auto i : decoding_graph.vertices()){
    //     for(auto j : decoding_graph.vertices()){

    //         auto ed = decoding_graph.get_edge(i,j);
    //         if(ed != nullptr){
    //             s++;
    //         }
    //     }
    // }
    // std::cout << "All EdgesX2: " << s << " Number of detectors: "<< n_detectors << " v = "<< decoding_graph.vertices().size()<< std::endl;
    // std::cout << "BOUNDRY NEIGHBORS = " << decoding_graph.adjacency_list(decoding_graph.get_vertex(BOUNDARY_INDEX)).size() << std::endl;


    // return statres;

#ifndef B_STAT_LER_USE_POLY
    mean_flip_prob = edge_list.size() / mean_flip_prob;
#endif
    std::discrete_distribution<> edge_dist(edge_probs.begin(), edge_probs.end());

    for (uint64_t s = 0; s < local_shots; s++) {
        for (uint i = 0; i < n_faults-1; i++) {
            // Add a random error.
            auto edge = edge_list[edge_dist(rng)];
            uint d1 = edge->detectors.first;
            uint d2 = edge->detectors.second;
            if (d1 != BOUNDARY_INDEX) {
                result_table[s][d1] ^= 1;
            }
            if (d2 != BOUNDARY_INDEX) {
                result_table[s][d2] ^= 1;
            }
            for (uint obs : edge->frames) {
                if (obs >= 0) {
                    result_table[s][n_detectors+obs] ^= 1;
                }
            }
        }
    }

    fp_t prev_logical_error_rate = 0.0;
    
    while (statres.n_logical_errors < 10 || DELTA(statres.logical_error_rate, prev_logical_error_rate) > 0.1) {
        uint64_t local_errors = 0;

        for (uint64_t s = 0; s < local_shots; s++) {
            // Add a random error.
            auto edge = edge_list[edge_dist(rng)];
            uint d1 = edge->detectors.first;
            uint d2 = edge->detectors.second;
            if (d1 != BOUNDARY_INDEX) {
                result_table[s][d1] ^= 1;
            }
            if (d2 != BOUNDARY_INDEX) {
                result_table[s][d2] ^= 1;
            }
            for (uint obs : edge->frames) {
                if (obs >= 0) {
                    result_table[s][n_detectors+obs] ^= 1;
                }
            }

            auto syndrome = _to_vector(result_table[s], n_detectors, n_observables);
            auto res = decoder->decode_error(syndrome);
            local_errors += res.is_logical_error;
        }

        uint64_t fault_errors;
        MPI_Allreduce(&local_errors, &fault_errors, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

#ifdef B_STAT_LER_USE_POLY
        fp_t log_fault_rate = log_poly[n_faults];
#else
        fp_t log_possible_faults = lgamma(edge_list.size()+1) 
                                    - lgamma(n_faults+1) - lgamma(edge_list.size()-n_faults+1);
        fp_t log_fault_rate = log_possible_faults
                                + n_faults*log(mean_flip_prob)
                                + (edge_list.size() - n_faults)*log(1-mean_flip_prob);
#endif
        if (world_rank == 0) {
            std::cout << "n_faults = " << n_faults << "\n"
                    << "\tlog fault rate = " << log_fault_rate << " (" << pow(M_E, log_fault_rate) << ")\n";
        }

        if (fault_errors > 0) {
            fp_t log_failure_rate = log(fault_errors) - log(shots_per_batch);
        
            prev_logical_error_rate = statres.logical_error_rate;
            statres.logical_error_rate += pow(M_E, log_failure_rate+log_fault_rate);
            statres.n_logical_errors += fault_errors;

            if (world_rank == 0) {
                std::cout << "\tlog failure rate = " << log_failure_rate << " (" << pow(M_E, log_failure_rate) << ")\n"
                        << "\tnumber of errors = " << statres.n_logical_errors << "\n";
                std::cout << "\tlogical error rate = " << statres.logical_error_rate << "\n";
            }
        }

        n_faults++;
    }
    if (world_rank == 0) {
    std::cout << "DONE!" << std::endl;
    }

    return statres;
}

void
generate_traces(std::string output_folder, const stim::Circuit& circuit, uint64_t shots, uint64_t shots_per_batch,
        uint64_t hw_cutoff, uint64_t base, uint64_t offset, std::mt19937_64& rng) 
{
    uint64_t fileno = base;
    const uint64_t n_results = circuit.count_detectors() + circuit.count_observables();
    const stim::simd_bits ref(n_results);

    uint64_t t = 0;
    stim::simd_bit_table filtered_samples(shots_per_batch, n_results);
    while (shots) {
        const uint64_t shots_this_round = shots < shots_per_batch ? shots : shots_per_batch;

        stim::simd_bit_table samples = stim::detector_samples(circuit, shots_this_round, false, true, rng);

        samples = samples.transposed();
        for (uint64_t s = 0; s < shots_this_round; s++) {
            if (samples[s].popcnt() > hw_cutoff) {
                filtered_samples[t++] |= samples[s];
            }
            if (t == shots_per_batch) {
                filtered_samples = filtered_samples.transposed();
                std::string filename = output_folder + "/shots_" + std::to_string(fileno) + ".dets";
                FILE * shots_out = fopen(filename.c_str(), "w");
                stim::write_table_data(shots_out, t, n_results, ref, filtered_samples, 
                                        stim::SampleFormat::SAMPLE_FORMAT_DETS, 'D', 'L', circuit.count_detectors());

                fclose(shots_out);
                fileno += offset;

                filtered_samples.clear();
                t = 0;
            }
        }

        shots -= shots_this_round;
    }
    filtered_samples = filtered_samples.transposed();
    std::string filename = output_folder + "/shots_" + std::to_string(fileno) + ".dets";
    FILE * shots_out = fopen(filename.c_str(), "w");
    stim::write_table_data(shots_out, t, n_results, ref, filtered_samples, 
                            stim::SampleFormat::SAMPLE_FORMAT_DETS, 'D', 'L', circuit.count_detectors());

    fclose(shots_out);
}

void
read_traces(std::string input_folder, Decoder * decoder, uint64_t max_shots_per_file, uint64_t base, uint64_t offset) {
    uint64_t fileno = base;

    const uint64_t n_detectors = decoder->circuit.count_detectors();
    const uint64_t n_observables = decoder->circuit.count_observables();

    while (true) {
        std::string filename = input_folder + "/shots_" + std::to_string(fileno) + ".dets";
        if (base == 0) {
            std::cout << "reading " << filename << "\n";
        }
        FILE * shots_in = fopen(filename.c_str(), "r");
        if (!shots_in) {
            break;
        }

        stim::simd_bit_table samples(max_shots_per_file, n_detectors+n_observables);
        uint64_t true_shots = stim::read_file_data_into_shot_table(shots_in, max_shots_per_file, n_detectors,
                                    stim::SampleFormat::SAMPLE_FORMAT_DETS, 'D', samples, true, 0, n_detectors,
                                    n_observables);
        if (base == 0) {
            std::cout << "\tshots in file = " << true_shots << "\n";
        }
        for (uint64_t s = 0; s < true_shots; s++) {
            if (base == 0 && s % 100000 == 0) {
                std::cout << "\tshot " << s << "\n";
            }
            auto syndrome = _to_vector(samples[s], n_detectors, n_observables);
            auto res = decoder->decode_error(syndrome);
            decoder->n_logical_errors += res.is_logical_error;
        }

        fileno += offset;
        fclose(shots_in);
    }
}

stim::Circuit
build_circuit(
    uint code_dist,
    fp_t error_mean,
    fp_t error_stddev,
    bool is_memory_z,
    bool is_rotated,
    bool both_stabilizers,
    uint8_t other_flags,
    uint rounds,
    fp_t clevel_error_mean,
    fp_t clevel_error_stddev,
    fp_t pauliplus_error_mean,
    fp_t pauliplus_error_stddev,
    fp_t round_dp_mean,
    fp_t sq_dp_mean,
    fp_t cx_dp_mean,
    fp_t reset_flip_mean,
    fp_t meas_flip_mean,
    fp_t round_dp_stddev,
    fp_t sq_dp_stddev,
    fp_t cx_dp_stddev,
    fp_t reset_flip_stddev,
    fp_t meas_flip_stddev,
    fp_t round_leak_mean,
    fp_t clifford_leak_mean,
    fp_t leak_transport_mean,
    fp_t round_leak_stddev,
    fp_t clifford_leak_stddev,
    fp_t leak_transport_stddev)
{
    if (rounds == 0) {
        rounds = code_dist;
    }
    std::string circ_type = (is_rotated ? "" : "un");
    circ_type += "rotated_memory_";
    circ_type += (is_memory_z ? "z" : "x");

    stim::CircuitGenParameters params(rounds, code_dist, circ_type);
    // Declare error rates.
    params.before_round_data_depolarization = __CHS(error_mean, clevel_error_mean, round_dp_mean);
    params.after_clifford_depolarization = __CHS(error_mean, clevel_error_mean, cx_dp_mean);
    params.after_clifford_sq_depolarization = __CHS(error_mean, clevel_error_mean, sq_dp_mean);
    params.after_reset_flip_probability = __CHS(error_mean, clevel_error_mean, reset_flip_mean);
    params.before_measure_flip_probability = __CHS(error_mean, clevel_error_mean, meas_flip_mean);

    params.before_round_data_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, round_dp_stddev);
    params.after_clifford_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, cx_dp_stddev);
    params.after_clifford_sq_depolarization_stddev = __CHS(error_stddev, clevel_error_stddev, sq_dp_stddev);
    params.after_reset_flip_probability_stddev = __CHS(error_stddev, clevel_error_stddev, reset_flip_stddev);
    params.before_measure_flip_probability_stddev = __CHS(error_stddev, clevel_error_stddev, meas_flip_stddev);

    params.before_round_leakage_probability = __CHS(error_mean, pauliplus_error_mean, round_leak_mean);
    params.after_clifford_leakage_probability = __CHS(error_mean, pauliplus_error_mean, clifford_leak_mean);
    params.after_clifford_leakage_transport = __CHS(error_mean, pauliplus_error_mean, leak_transport_mean);
    
    params.before_round_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, round_leak_stddev);
    params.after_clifford_leakage_probability_stddev = __CHS(error_stddev, pauliplus_error_stddev, clifford_leak_stddev);
    params.after_clifford_leakage_transport_stddev = __CHS(error_stddev, pauliplus_error_stddev, leak_transport_stddev);

    params.both_stabilizers = both_stabilizers;
    params.swap_lru = other_flags & BC_FLAG_SWAP_LRU_V1;
    params.swap_lru_with_no_swap = other_flags & BC_FLAG_SWAP_LRU_V2;
    params.initial_state_is_basis_1 = other_flags & BC_FLAG_INVERT_STATE;

    stim::Circuit circ = generate_surface_code_circuit(params).circuit;
    return circ;
}

}  // qrc
