/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#ifndef ASTREA_SIMULATOR_h
#define ASTREA_SIMULATOR_h

#include "astrea/mld_decoder.h"
#include "defs.h"
#include "decoding_graph.h"
#include "graph/dijkstra.h"
#include "mwpm_decoder.h"

#include <algorithm>
#include <array>
#include <deque>
#include <map>
#include <vector>
#include <utility>

//#define ASTREA_DEBUG
#define HW_CUTOFF 10

namespace qrc {
namespace astrea {

struct AstreaSimulatorParams {
    uint n_detectors;
    uint n_detectors_per_round;

    uint n_registers;
    uint bfu_fetch_width;
    uint bfu_compute_stages;
    uint bfu_priority_queue_size;

    uint32_t weight_filter_cutoff;

    bool use_mld;
};

typedef std::pair<uint, uint> addr_t;

class AstreaDeque;

class AstreaSimulator {
public:
    AstreaSimulator(DecodingGraph&, const stim::Circuit&, const AstreaSimulatorParams&);
    ~AstreaSimulator();

    void load_detectors(const std::vector<uint>&);
    void load_graph(DecodingGraph&, 
            const graph::PathTable<DecodingGraph::Vertex>&);
    void load_qubit_number(uint);
    void load_base_address(uint8_t bankgroup, uint8_t bank, uint32_t row_offset);

    void tick(void);
    void sig_end_round(uint=1);  // Number of rounds to jump. 
    
    bool is_idle(void);
    void force_idle(void);
    std::map<uint, uint> get_matching(void);

    void reset_stats(void);

    // Statistics
    uint64_t prefetch_cycles;
    uint64_t bfu_cycles;
    uint64_t cycles_to_converge;
    uint64_t valid_weights_after_filter;
protected:
    void tick_prefetch(void);
    void tick_bfu(void);

    void tick_bfu_compute(uint stage);
    void tick_bfu_sort(uint stage);
    void tick_bfu_fetch(void);

    bool access(addr_t, bool set_evictable_on_hit);
    void update_state(void);

    void clear(void);

    enum class State { prefetch, bfu, idle };

    struct Register {
        addr_t address;
        uint64_t last_use;
        bool valid;
    };

    struct DequeEntry {
        std::map<uint, uint> running_matching;
        uint32_t matching_weight;
        uint next_unmatched_index;
        
        bool operator<(const DequeEntry& other) const {
            return matching_weight/running_matching.size() 
                < other.matching_weight/other.running_matching.size();
        }
    };

    struct BFUPipelineLatch {
        std::deque<std::pair<uint, uint32_t>> proposed_matches;
        DequeEntry base_entry;
        bool stalled;
        bool valid;
    };

/* Microarchitectural components.*/
    // Global memory
    std::vector<Register> register_file;    
    std::vector<uint> detector_vector_register; // Holds the detectors from the
                                                // current syndrome.
    uint32_t mean_weight_register;
    uint access_counter;

    // Prefetch
    uint min_detector_register;
    uint major_detector_register;
    std::map<uint, uint> minor_detector_table;

    // There are fetch_width hardware deques of size 2*fetch_width
    std::vector<AstreaDeque> hardware_deques;

    // Size of latches is fetch_width by (1 + 4 + 1)
    // There is one FETCH stage,
    //          four SORT stages,
    //          and 1 compute stage.
    std::vector<std::vector<BFUPipelineLatch>> bfu_pipeline_latches;   
    
    std::deque<addr_t> replacement_queue;

    // MWPM variant
    DequeEntry best_matching_register;
    // MLD variant
    std::array<DequeEntry, 2> mld_bm_register;
    std::array<fp_t, 2> mld_prob_register;
    // Global state machine
    State state; 
    bool bfu_idle;

    /* Data */
    DecodingGraph graph;
    graph::PathTable<DecodingGraph::Vertex> path_table;
    /* Configuation parameters. */
private:
    uint curr_max_detector;
    uint32_t weight_filter_cutoff;

    uint n_detectors;
    uint n_detectors_per_round;
    uint bfu_fetch_width;
    uint bfu_compute_stages; 
    uint curr_qubit;
    bool has_boundary;

    bool use_mld;

    uint32_t cycles_after_last_converge;

    MWPMDecoder * hw6decoder;

    friend class AstreaDeque;
};

// Fixed-size priority queue.
class AstreaDeque {
public:
    AstreaDeque(uint);

    void push(AstreaSimulator::DequeEntry);
    AstreaSimulator::DequeEntry top(void);
    void pop(void);
    uint size(void);
    bool empty(void);
    void clear(void);
private:
    void downheap(uint);
    void upheap(uint);

    uint parent(uint);
    uint left(uint);
    uint right(uint);

    std::vector<AstreaSimulator::DequeEntry> backing_array;

    uint max_size;

    friend class AstreaSimulator;
};

uint bound_detector(uint, uint n_detectors);
uint unbound_detector(uint, uint n_detectors);

} // astrea
} // qrc

#endif
