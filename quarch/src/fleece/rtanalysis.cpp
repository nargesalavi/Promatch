/*
 *  author: Suhas Vittal
 *  date:   23 February 2023
 * */

#include "fleece/rtanalysis.h"

namespace qrc {
namespace fleece {

RealTimeAnalyzer::RealTimeAnalyzer(
        const stim::Circuit& circuit, std::mt19937_64& rng)
    :circuit(circuit),
    leakage_persist_table(),
    lattice_graph(to_lattice_graph(circuit)),
    sim(circuit.count_qubits(), MAX_SHOTS, SIZE_MAX, rng),
    leakage_occurrence_table()
{
    leakage_persist_table.fill(std::make_tuple(0, 0, 0));
}

void
RealTimeAnalyzer::setup() {
    leakage_occurrence_table.clear();
    sim.reset_all();
}

bool
RealTimeAnalyzer::step(uint64_t shots) {
    if (sim.cycle_level_simulation(circuit)) {
        flush_table();
        return false;
    }

    for (uint64_t s = 0; s < shots; s++) {
        for (auto v : lattice_graph.vertices()) {
            if (v->qubit < 0) {
                continue;
            }

            if (sim.leakage_table[v->qubit][s]) {
                increment(v, s);
            } else {
                read_and_flush(v, s);
            }
        } 
    }
    
    return true;
}

void
RealTimeAnalyzer::increment(LatticeGraph::Vertex * v, uint64_t s) {
    if (!leakage_occurrence_table.count(s)) {
        leakage_occurrence_table[s] = std::map<LatticeGraph::Vertex*, uint>();
    }

    if (!leakage_occurrence_table[s].count(v)) {
        leakage_occurrence_table[s][v] = 0;
    }
    leakage_occurrence_table[s][v]++;
}

void
RealTimeAnalyzer::read_and_flush(LatticeGraph::Vertex * v, uint64_t s) {
    uint length = leakage_occurrence_table[s][v];
    persist_t& entry = leakage_persist_table[length];

    std::get<0>(entry)++;
    if (v->is_data) {
        std::get<1>(entry)++;
    } else {
        std::get<2>(entry)++;
    }

    leakage_occurrence_table[s][v] = 0;
}

void
RealTimeAnalyzer::flush_table() {
    for (auto pair1 : leakage_occurrence_table) {
        uint64_t s = pair1.first;
        for (auto pair2 : pair1.second) {
            auto v = pair2.first;
            read_and_flush(v, s);
        }
    }
}

}   // fleece
}   // qrc
