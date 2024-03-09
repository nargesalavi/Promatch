
/*
    author: Narges Alavisamani
*/

#ifndef BB_BENCHMARK_h
#define BB_BENCHMARK_h


#include <benchmark.h>
#include "bucket_based_simulation.h"
#include "experiments.h"
#include <predecoder.h>
#include <mwpm_decoder.h>
#include "syndrome_surgery.h"
#include <astrea.h>
#include <iomanip>
#include "predecoding_graph_helper.h"

#include <functional>

#include <stdexcept>


// const uint64_t arr_N = 300;

fp_t binomial_probability(uint k, uint n, fp_t p);

fp_t poisson_probability(fp_t lambda, int k);

void HW_distributions(std::ofstream& out, uint64_t max_shot, uint distance, 
    fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint64_t hshots_replc, uint64_t lshots_rplac);

void edge_distributions(std::ofstream& out, uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint64_t hshots_replc, uint64_t lshots_rplac);

void test_sim(std::ofstream& out, std::string syndromes_directory, uint64_t total_shots, uint distance, 
fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry);

void bb_ler_calculation_and_stats(std::ofstream& out, std::string syndromes_directory, uint64_t shots_per_k, uint distance, 
fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry);

void printing_samples_of_syndromes_and_matchings(uint64_t shots_per_k, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool print_ec_length = false);

void bb_ler_calculation_AstreaG(std::ofstream& out, std::string syndromes_directory, uint distance, 
fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry);

void test_b_statistical_ler(uint64_t shots_per_k, uint distance, fp_t physcial_error, fp_t meas_er);

void bb_ler_calculation(std::string decoder_name, std::ofstream& out, uint64_t shots_per_k, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k = 1, uint64_t max_k = 20, uint64_t SHOTS_PER_BATCH = 1'000'000);

void test_b_statistical_ler_AstreaG(uint64_t shots_per_k, uint distance, fp_t physcial_error, fp_t meas_er);

void error_chains_distribution(uint64_t shots_per_batch, bool use_mpi, uint n_faults, uint distance,  fp_t physcial_error, fp_t meas_er
                        , uint64_t min_k, uint64_t max_k);
void test_groups( uint distance, fp_t physcial_error, fp_t meas_er);

void test_finding_fast_group(uint64_t shots_per_k, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, std::string decoder_name="MWPM");

void printing_incorrect_predictions(uint64_t shots_per_k, uint distance, 
                            fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k,
                            std::string decoder_name="MWPM");

void test_priority_groups(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name="MWPM");

void test_size_of_predecoding(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, std::string decoder_name="MWPM");

void bb_ler_calc_promatch(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name="MWPM", bool generalized = false,
uint round_n = 0, bool save_syndromes = false, std::string syndrome_folder_name="dummyfolder", std::string syndrome_clock_challenging_file="dummyfolder_clock",  uint64_t hshots_replc = 20'000'000, uint64_t lshots_rplac = 1'000'000);

void degree_distribution(std::ofstream& out, uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint64_t hshots_replc = 20'000'000, uint64_t lshots_rplac = 1'000'000);

void check_the_incorrect_syndrome(uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, std::ofstream& out, bool generalized, std::string decoder_name="MWPM", bool write_file = false);

void stats_of_incorrect_syndromes(uint distance,
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool generalized = false, std::string decoder_name="MWPM");

void bb_decoder_ler_calc(uint64_t max_shot, uint distance,fp_t physcial_error, fp_t meas_er, uint64_t min_k, 
    uint64_t max_k, bool& print_time, uint round_n = 0, bool save_syndromes = false, std::string syndrome_folder_name="dummyfolder",  
    std::string decoder_name="MWPM", uint64_t hshots_replc = 20'000'000, uint64_t lshots_rplac = 1'000'000);
    
void bb_predecoder_ler_calc(uint64_t max_shot, uint distance,fp_t physcial_error, fp_t meas_er, uint64_t min_k, 
    uint64_t max_k, bool& print_time, uint round_n = 0, bool save_syndromes = false, std::string syndrome_folder_name="dummyfolder",  
    std::string decoder_name="MWPM", std::string predecoder_name="CLIQUE", uint64_t hshots_replc = 20'000'000, uint64_t lshots_rplac = 1'000'000);

void bb_avg_ec1(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name="MWPM");

void sample_test_degree_priority(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name="MWPM");

void ler_predecoders(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, fp_t threshold_scale, bool& print_time, std::string decoder_name, bool generalized, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string syndrome_clock_challenging_file, uint64_t hshots_replc, uint64_t lshots_rplac);    

void PAG_synergies(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name, bool generalized, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string syndrome_clock_challenging_file, uint64_t hshots_replc, uint64_t lshots_rplac);

void print_mwpm_mistakes(uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, std::ofstream& out, bool generalized, std::string decoder_name="MWPM", bool write_file = false);

 void testing_all(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, fp_t threshold_scale, bool& print_time, std::string decoder_name, bool generalized, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string syndrome_clock_challenging_file, uint64_t hshots_replc, uint64_t lshots_rplac, bool only_important_bucket = false);    


#endif