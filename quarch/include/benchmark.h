/*
 *  author: Suhas Vittal
 *  date:   21 August 2022
 * */

#ifndef BENCHMARK_h
#define BENCHMARK_h

#include <stim.h>

#include "benchmark/statbench.h"
#include "defs.h"
#include "decoder.h"

#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <utility>

#include <math.h>
#include <mpi.h>

#define MAX_SHOTS 100000

namespace qrc {

fp_t 
min(const std::vector<fp_t>&);
fp_t
max(const std::vector<fp_t>&);
fp_t
mean(const std::vector<fp_t>&);
fp_t
stddev(const std::vector<fp_t>&);

void
b_decoder_ler(Decoder*, uint64_t shots, std::mt19937_64&, bool save_per_shot_data=false);
/*
 *  Pre-condition: MPI is initialized before call and exited after call.
 * */
benchmark::StatisticalResult
b_statistical_ler(Decoder*, uint64_t shots, std::mt19937_64&, bool use_mpi=true, uint n_faults=1, uint distance = 5);
/*
 *  Save syndromes above some given Hamming weight.
 * */
void
generate_traces(std::string output_folder, const stim::Circuit&, uint64_t shots, uint64_t shots_per_batch,
                uint64_t hw_cutoff, uint64_t base, uint64_t offset, std::mt19937_64&);
void
read_traces(std::string input_folder, Decoder*, uint64_t max_shots_per_file, uint64_t base, uint64_t offset);

#define BC_FLAG_SWAP_LRU_V1     0x1
#define BC_FLAG_SWAP_LRU_V2     0x2
#define BC_FLAG_INVERT_STATE    0x4

stim::Circuit
build_circuit(
    // Required
    uint code_dist, 
    fp_t error_mean,
    fp_t error_stddev,
    // Optionals
    bool is_memory_z=true,
    bool is_rotated=true,
    bool both_stabilizers=false,
    uint8_t other_flags=0,
    // Level 1 Specificity
    uint rounds=0,
    fp_t clevel_error_mean=-1,
    fp_t clevel_error_stddev=-1,
    fp_t pauliplus_error_mean=-1,
    fp_t pauliplus_error_stddev=-1,
    // Level 2 Specificity
    fp_t round_dp_mean=-1,
    fp_t sq_dp_mean=-1,
    fp_t cx_dp_mean=-1,
    fp_t reset_flip_mean=-1,
    fp_t meas_flip_mean=-1,
    fp_t round_dp_stddev=-1,
    fp_t sq_dp_stddev=-1,
    fp_t cx_dp_stddev=-1,
    fp_t reset_flip_stddev=-1,
    fp_t meas_flip_stddev=-1,
    fp_t round_leak_mean=-1,
    fp_t clifford_leak_mean=-1,
    fp_t leak_transport_mean=-1,
    fp_t round_leak_stddev=-1,
    fp_t clifford_leak_stddev=-1,
    fp_t leak_transport_stddev=-1);

}  // qrc

#endif
