/*
 *  author: Suhas Vittal
 *  date:   23 February 2023
 * */

#ifndef FLEECE_RTANALYSIS_h
#define FLEECE_RTANALYSIS_h

#include <stim.h>

#include <array>
#include <map>
#include <tuple>

#include "benchmark.h"
#include "defs.h"
#include "lattice_graph.h"

namespace qrc {
namespace fleece {

/*
 *  Want to identify:
 *      (1) Average number of rounds that occur until leakage is eliminated for
 *          --> data qubit
 *          --> parity qubit.
 * */

typedef std::tuple<uint64_t, uint64_t, uint64_t> persist_t;

class RealTimeAnalyzer {
public:
    RealTimeAnalyzer(const stim::Circuit&, std::mt19937_64&);

    void setup(void);
    bool step(uint64_t shots);  // shots < MAX_SHOTS (100'000).

    const stim::Circuit circuit;

    std::array<persist_t, 100> leakage_persist_table;

    void increment(LatticeGraph::Vertex*, uint64_t s);
    void read_and_flush(LatticeGraph::Vertex*, uint64_t s);
    void flush_table(void);
protected:
    LatticeGraph lattice_graph;
    stim::FrameSimulator sim;
private:
    std::map<uint64_t, std::map<LatticeGraph::Vertex*, uint>> leakage_occurrence_table;
};

} // fleece
} // qrc

#endif
