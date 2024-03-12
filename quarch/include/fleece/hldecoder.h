/*
 *  author: Suhas Vittal
 *  date:   15 March 2023
 * */

#ifndef FLEECE_HLDECODER_h
#define FLEECE_HLDECODER_h

#include "decoder.h"
#include "fleece/lattice_graph.h"
#include "mwpm_decoder.h"

namespace qrc {
namespace fleece {

class HLDecoder : public MWPMDecoder {
public:
    HLDecoder(const stim::Circuit&);

    DecoderShotResult decode_error(const std::vector<uint8_t>&, const std::vector<uint8_t>&);
protected:
    graph::PathTable<DecodingGraph::Vertex> base_path_table;
    fleece::LatticeGraph lattice_graph;
};

} // fleece
} // qrc

#endif
