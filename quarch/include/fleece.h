/*
 *  author: Suhas Vittal
 *  date:   6 December 2022
 * */

#ifndef FLEECE_h
#define FLEECE_h

#include <stim.h>

#include <PerfectMatching.h>

#include "decoder.h"
#include "defs.h"
#include "fleece/lattice_graph.h"
#include "fleece/rtanalysis.h"
#include "graph/dijkstra.h"
#include "mwpm_decoder.h"

#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set> 
#include <string>
#include <vector>
#include <utility>

/*
 *  FLEECE = Fast LEakagE CorrEction
 *
 *  TODO: Implement basis-aware FLEECE. Only supports
 *  Z basis operations at the moment.
 * */

#define EN_STATE_DRAIN      0x01
#define EN_TMP_STAB_EXT     0x02 // enables EN_STATE_DRAIN by default.
#define EN_SWAP_LRU         0x04
#define NO_MITIGATION       0x08
#define NO_RESTART          0x10

namespace qrc {

class Fleece {
public:
    Fleece(const stim::CircuitGenParameters&,
            uint16_t flags,
            std::mt19937_64& rng,
            char reset_basis='Z',
            char output_basis='Z');
    ~Fleece();

    stim::simd_bit_table create_syndromes(uint64_t shots, 
            bool maintain_failure_log=false,
            bool record_in_rtanalyzer=false);

    std::string failure_log;

    fleece::RealTimeAnalyzer * rtanalyzer;

    uint64_t n_restarts;
    std::vector<uint64_t> parity_leakage_population;
    std::vector<uint64_t> data_leakage_population;
private:
    void compute_optimal_swap_set(void);

    void write_leakage_condition_to_log(std::string&);
    void write_aliases_to_log(std::string&, 
            const std::map<fleece::LatticeGraph::Vertex*, fleece::LatticeGraph::Vertex*>&);
    void write_syndrome_to_log(std::string&, const std::vector<uint8_t>&, const std::vector<uint8_t>&);

    void apply_reset(uint32_t qubit, bool add_error=true);
    void apply_round_start_error(uint32_t qubit, fp_t dp_error_mult=1.0);
    void apply_H(uint32_t qubit, bool add_error=true);
    void apply_CX(uint32_t, uint32_t);
    void apply_measure(uint32_t qubit, bool add_error=true);
    void apply_SWAP(uint32_t, uint32_t, bool add_error=true);

    const stim::CircuitGenParameters circuit_params;
    const uint16_t flags;
    const char reset_basis;
    const char output_basis;

    bool leakage_enabled;

    stim::Circuit base_circuit;
    stim::FrameSimulator * sim;

    fleece::LatticeGraph lattice_graph;
    graph::PathTable<fleece::LatticeGraph::Vertex> lattice_path_table;
    std::vector<fleece::LatticeGraph::Vertex*> data_qubits;
    std::vector<fleece::LatticeGraph::Vertex*> parity_qubits;

    std::map<fleece::LatticeGraph::Vertex*, fleece::LatticeGraph::Vertex*> swap_set;
    fleece::LatticeGraph::Vertex * unlucky_data_qubit;

    std::map<uint32_t, fp_t> rounddp_table;
    std::map<uint32_t, fp_t> clifforddp1_table;
    std::map<std::vector<uint32_t>, fp_t> clifforddp2_table;
    std::map<uint32_t, fp_t> premeasflip_table;
    std::map<uint32_t, fp_t> postresflip_table;

    std::map<uint32_t, fp_t> roundleak_table;
    std::map<std::vector<uint32_t>, fp_t> cliffordleak_table;
    std::map<std::vector<uint32_t>, fp_t> leaktransport_table;

    std::mt19937_64 rng;
};

} // qrc

#endif
