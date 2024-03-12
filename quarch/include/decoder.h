/*
 *  author: Suhas Vittal
 *  date:   5 August 2022
 * */

#ifndef DECODER_h
#define DECODER_h

#include <stim.h>

#include "defs.h"
#include "decoding_graph.h"

#include <string>
#include <utility>

namespace qrc {

/* Benchmark functions and structures*/
struct DecoderShotResult {
    fp_t execution_time;
    fp_t memory_overhead;
    bool is_logical_error;
    std::vector<uint8_t> correction;
    std::map<uint, uint> matching;
    fp_t weight;
};

class Decoder {
public: 
    Decoder(const stim::Circuit&);
    virtual ~Decoder() {}

    virtual DecoderShotResult decode_error(const std::vector<uint8_t>&) =0;

    virtual std::string name(void) =0;
    virtual bool is_software(void) =0;

    virtual void clear_stats();

    // General decoder statistics
    // Functions and member variables.

    typedef std::pair<std::vector<uint8_t>, fp_t> SyndromeEvent;
    /*
     * Gets all nonzero syndromes and the time required to complete
     * each of them. Each output pair has format (syndrome, time_taken).
     * */
    std::vector<SyndromeEvent> nonzero_syndromes_and_time_taken();
    /*
     *  Gets all nonzero syndromes completd within the given time.
     *  Important for real time decoding. Each output pair has 
     *  the format: (syndrome, time taken).
     * */
    std::vector<SyndromeEvent>
        nonzero_syndromes_completed_within(fp_t max_time,
                uint32_t& n_nonzero_syndromes);
    /*
     *  Gets the maximum hamming weight that the decoder will
     *  always complete within max_time.
     * */
    uint max_hamming_wgt_completed_within(fp_t max_time);
    /*
     *  Gets the potential SRAM cost if all data structures were
     *  implemented in hardware. Returns amount in bytes.
     * */
    virtual uint64_t sram_cost(void);
    
    std::vector<std::vector<uint8_t>> syndromes;
    std::vector<fp_t> execution_times;  // in nanoseconds
    std::vector<fp_t> memory_overheads; // in bytes
    uint64_t n_logical_errors;
    fp_t mean_execution_time;
    fp_t max_execution_time;
    fp_t max_execution_time_for_correctable;
    std::map<uint, fp_t> hamming_weight_dist;
    // Benchmarking circuit.
    stim::Circuit circuit;
protected:
    // Member variables.
    DecodingGraph graph;
    // Give all benchmarking functions
    // friendship.
    friend void 
        b_decoder_ler(Decoder*, uint64_t, std::mt19937_64&, bool);
};

uint32_t
syndrome_to_int(const std::vector<uint8_t>&, uint size);
bool
is_logical_error(const std::vector<uint8_t>& correction,
        const std::vector<uint8_t>& syndrome,
        uint n_detectors, uint n_observables);

std::vector<uint8_t> 
_to_vector(const stim::simd_bits_range_ref&,
        uint n_detectors, uint n_observables);

} // qrc

#endif
