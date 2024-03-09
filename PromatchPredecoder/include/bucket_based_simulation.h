/*
    author: Narges Alavisamani
*/

#ifndef BUCKET_BASED_BENCHMARK_h
#define BUCKET_BASED_BENCHMARK_h
#include <algorithm>
#include "decoding_graph.h"
#include "error_event.h"
#include "syndrome_surgery.h"
#include <typeinfo>


namespace bbsim{


class BucketBasedSim {
    public:
        BucketBasedSim(const stim::Circuit& circ, uint k_max_);
        ~BucketBasedSim();

        void print_decoding_graph();
        void print_sorted_event_vector();
        void print_path_matrix();
        void print_weight_matrix();
        void print_vertices();
        void print_path_list();

        std::vector<std::vector<uint8_t>> create_syndromes(uint k, uint64_t shots_per_k, std::mt19937_64& rng, bool print_log =  false);
        stim::simd_bit_table create_syndromes_simd(uint k, uint64_t shots_per_k, std::mt19937_64& rng, bool print_log);


        fp_t return_probab_sum();

        uint64_t number_of_events;
        std::vector<fp_t> list_of_error_probabs;

        fp_t event_probability(int k);
        uint k_max;

        std::vector<std::vector<int>> combinations(const std::vector<fp_t>& nums, int k);
        void  combinationsHelper(const std::vector<fp_t>& nums, int k, int start, std::vector<int>& temp, std::vector<std::vector<int>>& result);

        qrc::DecodingGraph decoding_graph;
        std::vector<fp_t> prob_k;
        std::vector<uint64_t> shots_per_k;

        void calculate_prob_k();
        void set_shots_per_k(fp_t expected_ler, uint64_t max_n_iteration, bool set_shots_to_all = true, uint64_t filling_n_iterations = 20'000'000, uint64_t min_n_iteration = 10'000'000);

        uint64_t n_detectors;
        uint64_t n_observables;


    private:
        
        fp_t probab_sum;

        std::map<fp_t, uint64_t> probability_group;
        
        std::vector<ErrorEvent*> sorted_events_vector;

        void create_the_sorted_event_vector_calc_prob_k();

        

};

}

#endif