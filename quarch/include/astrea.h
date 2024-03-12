/*
 *  author: Suhas Vittal
 *  date:   22 August 2022
 * */

#ifndef ASTREA_h
#define ASTREA_h

#include "astrea/mld_decoder.h"
#include "benchmark.h"
#include "astrea/simulator.h"

#include <algorithm>
#include <array>
#include <deque>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stack>
#include <string>
#include <tuple>
#include <vector>

#include <math.h>

namespace qrc {

struct AstreaParams {
    // Fetch width for brute force unit
    uint bfu_fetch_width;             
    uint bfu_compute_stages;
    uint bfu_priority_queue_size;
    // Memory Parameters
    uint n_registers;

    bool use_mld;
    
    fp_t main_clock_frequency;   // in Hz
};

class Astrea : public astrea::MLDDecoder {
public:
    Astrea(const stim::Circuit, 
            uint n_detectors_per_round, 
            uint32_t weight_filter_cutoff,
            const AstreaParams&,
            fp_t time_limit=1000);
    ~Astrea();

    DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
    std::string name(void) override;
    bool is_software(void) override;

    uint64_t sram_cost(void) override;

    // Statistics
    uint64_t n_nonzero_syndromes;
    uint64_t n_hhw_syndromes;
    uint64_t total_bfu_cycles;
    uint64_t total_prefetch_cycles;
    uint64_t total_cycles_to_converge;
    fp_t total_logfilter_savings;
    uint64_t max_bfu_cycles;
    uint64_t max_prefetch_cycles;
    uint64_t max_cycles_to_converge;
    fp_t min_filter_savings;
    uint64_t max_hamming_weight;
    // More statistics are in the simulator.
    astrea::AstreaSimulator * simulator;
private:
    uint n_rounds;
    fp_t main_clock_frequency;
    fp_t time_limit;
    bool use_mld;
    
    MWPMDecoder baseline;
};

}  // qrc

#endif
