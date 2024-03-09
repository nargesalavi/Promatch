/*
    author: Narges Alavisamani
*/

#ifndef PREDECODER_h
#define PREDECODER_h

#include <stim.h>
#include <mwpm_decoder.h>
#include <benchmark.h>
#include <defs.h>
#include <mpi.h>
#include <math.h>
#include <random>
#include <cfloat>

#include <algorithm>
#include <fstream>
#include <iostream>
#include "decoder.h"
#include "decoding_graph.h"
#include "predecoding_graph_helper.h"
#include "syndrome_surgery.h"

#include <stdexcept>


// For 6, 8 cycles (limit = 250 - 10 - 8 = 232); for 8, 20(limit = 250 - 10 - 20 = 220); for 10, 114 cycles (limit = 250 - 10 - 114 = 126);
#define MAX_BF6_HW 6
#define MAX_BF6_CYCLE 232

#define MAX_BF8_HW 8
#define MAX_BF8_CYCLE 220

#define MAX_BF10_HW 10
#define MAX_BF10_CYCLE 126

//typedef struct qrc::DecoderShotResult DecoderShotResult_t;


namespace qpd {

struct NoneIsolatedInf{
    
    std::vector<qrc::DecodingGraph::Vertex*> v1_adjacent_flipped_parity;
    std::vector<qrc::DecodingGraph::Vertex*> v2_adjacent_flipped_parity;

    std::vector<fp_t> v1_adjacent_prob;
    std::vector<fp_t> v2_adjacent_prob;
};

struct  SingleMatchingInfo{
    uint first;
    uint second;
    fp_t chain_probability;
    uint length;
    int8_t isolated = -1;
    NoneIsolatedInf* non_isolated_inf = nullptr;
    
};

struct PredecoderShotResult{
    std::vector<uint8_t> post_syndrome;
    std::vector<uint8_t> prematch_syndrome;
    std::map<uint, uint> pre_matches;
};

struct PredecoderShotResult_v2{
    std::vector<uint8_t> post_syndrome;
    std::vector<std::pair<uint, uint>> pre_matches;
    std::map<std::pair<uint, uint>, uint> edge_counting;

};

struct PostPredecoderResult{
    std::vector<uint8_t> correction;
    std::vector<uint8_t> post_syndrome;
};

struct PredecoderSubgraph{
    std::map<qrc::DecodingGraph::Vertex*, std::vector<qrc::DecodingGraph::Vertex*>> adjacency_matrix;
    std::vector<uint> detector_list;
    uint distance_with_neighbor;
    std::map<qrc::DecodingGraph::Vertex*, std::vector<qrc::DecodingGraph::Edge*>> edge_lists;
};

struct EnsembleEntry {
    wgt_t predecoding_weight;
    std::vector<uint8_t> syndrome;
    std::map<uint, uint> pre_matching;
};

struct PredecoderDecoderShotResult{
    bool is_logical_error;
    std::vector<uint8_t> correction;
    std::map<std::pair<uint, uint>, uint> correcting_edge;
    std::vector<std::pair<uint,uint>> pre_matching;
    std::vector<std::pair<uint,uint>> dec_matching;
    fp_t weight;
};

class Predecoder: public qrc::Decoder{
    public:
        Predecoder(const stim::Circuit& circ, bool general_ = false);
    
        ~Predecoder(){}
        //std::vector<uint8_t> syndrome;
        qrc::DecodingGraph decoding_graph;
        uint n_detectors;
        uint n_observables;
        qrc::Decoder* decoder;

        //When ready is true it means that hw is low enough(6 or less) for the maind decoder
        bool ready;
        bool exist_degree_one_no_poor;
        bool exist_multi_degree_no_poor;
        bool general;

        // prioritized fast pairs contains edges that are
        // in the matching and classified into 6 groups
        std::vector<std::vector<qrc::DecodingGraph::Edge*>> pf_pairs;
        std::vector<qrc::DecodingGraph::Edge*> predecoding_edges;
        std::map<uint,uint> adaptive_predecoding_map;

        // prioritized fast groups contains vertices that are
        // predecoded and classified into 6 groups. Vertices are associated 
        // vertices of the edges in pf_pair
        std::vector<std::vector<qrc::DecodingGraph::Vertex*>> pf_groups;
        std::vector<qrc::DecodingGraph::Vertex*> predecoding_vertices;

        std::vector<qrc::DecodingGraph::Vertex*> degree_zero_vertices;
        
        uint number_of_edges;
        uint starting_number_of_edges;
        std::array<bool,6> reached_stage;


        virtual qrc::DecoderShotResult decode_error(const std::vector<uint8_t>&) override;
        std::string name(void) override;
        bool is_software(void) override;
        void set_decoder(qrc::Decoder* decoder_);
        
        // This function counts stand alone degree one error chains
        uint count_degree_one();

        //void send_sydrome_to_predecoder(std::vector<uint8_t> syndrome_);
        PredecoderShotResult matching_ec_length_one(const std::vector<uint8_t>& syndrome_, bool add_boundry = true);
        //MWPMDecoder::get_correction_from_matching(const std::map<uint, uint>& matching) 
        //Bunch of printing stuffs for validation
        // For checking the correctness
        void print_graph();
        void print_subgraph(const std::vector<uint8_t>& syndrome_);

        std::vector<uint8_t> get_correction_from_matching(const std::map<uint, uint>& matching);

        PostPredecoderResult predecode_error(PredecoderShotResult prematching_result);

        void print_shot_info(const std::vector<uint8_t>& syndrome,
                    qpd::PredecoderShotResult& predecoder_results,  qrc::DecoderShotResult &predecoder_decoder_results, 
                    qrc::DecoderShotResult &decoder_results, int n_1s, bool only_predecoder_matches=true, bool full_info=false, bool print_syndromes = false);

        PredecoderShotResult matching_ec(uint ec_length, const std::vector<uint8_t>& syndrome_, bool add_boundry = true);

        const PredecoderSubgraph get_adjacency_matrix();

        void update_adjacency_matrix_and_detector_list(const std::vector<uint8_t>& syndrome_, uint ec_length, bool add_boundry = true);
        void update_adjacency_matrix_and_detector_list_from_flipped_vec(const std::vector<uint16_t>& data, uint ec_length, bool add_boundry = true);
        void is_matched_parity_isolated(SingleMatchingInfo& m_info);
        std::vector<qrc::DecodingGraph::Vertex*> get_fast_matching_group(const std::vector<uint8_t>& syndrome_);
        std::vector<qrc::DecodingGraph::Vertex*> get_fast_matching_group_v2(const std::vector<uint8_t>& syndrome_);
        std::vector<qrc::DecodingGraph::Vertex*> get_fast_matching_group_all_neighbors(const std::vector<uint8_t>& syndrome_);

        std::vector<uint8_t> set_the_priorities(const std::vector<uint8_t>& syndrome_, bool all_in = true);
        std::vector<uint8_t> prioritize_and_set_potential_matching(const std::vector<uint8_t>& syndrome_);

        std::vector<EnsembleEntry> create_syndromes_ensembles(const std::vector<uint8_t>& syndrome_, bool all_in = true);
        std::vector<EnsembleEntry> create_syndromes_ensembles_upair_3(const std::vector<uint8_t>& syndrome_, bool all_in = true);
        std::vector<EnsembleEntry> create_syndromes_ensembles_unpair_1(const std::vector<uint8_t>& syndrome_, bool all_in = true);
        std::vector<EnsembleEntry> create_syndromes_ensembles_unpair_2(const std::vector<uint8_t>& syndrome_, bool all_in = true);
        std::vector<EnsembleEntry> create_syndromes_no_unpairing(const std::vector<uint8_t>& syndrome_, bool all_in = true);

        std::vector<EnsembleEntry> create_syndromes_selected_unpairing(const std::vector<uint8_t>& syndrome_);
        PredecoderDecoderShotResult smith_decode_error(const std::vector<uint8_t>& syndrome_);
    
        
        
        void stage1_decoding(std::vector<uint8_t>& syndrome_, uint round);
        void stage2_decoding(std::vector<uint8_t>& syndrome_, uint& round, bool no_poor);
        void stage3_decoding(std::vector<uint8_t>& syndrome_, uint& round, bool no_poor);
        void stage0_decoding(std::vector<uint8_t>& syndrome_, uint& round, bool no_poor);
        
    

        qrc::DecoderShotResult ensemble_decode_error(const std::vector<uint8_t>& syndrome_, uint max_unpair_number = 3, bool all_in = true);
        fp_t calc_matching_weight(std::map<uint,uint> );

        void prioritize_edges_by_degree(const std::vector<uint8_t>& syndrome_ , bool all_in = true);

        qrc::DecoderShotResult adaptively_decode_error(const std::vector<uint8_t>& syndrome_);

        PredecoderShotResult_v2 smith_predecoder(const std::vector<uint8_t>& syndrome_, bool add_boundry );
        std::vector<std::pair<uint,uint>> degree_zero_sorted_queue_paths();

        PredecoderShotResult_v2 clique_predecoder(const std::vector<uint8_t>& syndrome_, bool add_boundry );
        PredecoderDecoderShotResult clique_decode_error(const std::vector<uint8_t>& syndrome_);
        uint print_paths_matching(std::map<uint, uint> matching);
        fp_t graph_sparcity_calculator(std::vector<uint16_t> compressed_synd, fp_t& sparcity_unweighted, uint distance);


        void reset_clocks();
        uint total_cycles;
        uint total_cycles_parallel_updating;
        uint updating_clk_parallel;
        uint updating_clk;
        uint number_of_rounds;
        uint check_condition;
        uint number_paths_singleton_matching_candidate;

        std::vector<uint8_t> syndrome_register;
        void update_degrees();
        std::vector<std::vector<uint8_t>> simulator_adjacenct_matrix;
        void set_simulator_adjacenct_matrix();

        uint MAX_BF_HW;

        // in Hardware it is a register that has a valid bit that is connected to 
        // syndrome register. So, it is valid if the bit is flipped.
        std::map<uint, uint> detector_degree_map;
        


        


        
    private:
        
        PredecoderSubgraph prematching_subgraph;
        std::vector<qrc::DecodingGraph::Edge*> sorted_predecoding_edges;





    protected:
        qrc::graph::PathTable<qrc::DecodingGraph::Vertex> path_table;


};


} //namesapce qpd

#endif