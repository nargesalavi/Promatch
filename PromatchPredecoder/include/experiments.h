/*
    author: Narges Alavisamani
*/

#ifndef EXPERIMENTS_h
#define EXPERIMENTS_h
/*
 *  author: Narges Alavisamani
 *  date:   6 Feburary 2023
 * */
 
#include <stim.h>
#include <mwpm_decoder.h>
#include <astrea.h>
#include <benchmark.h>
#include <defs.h>
#include <math.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <predecoder.h>
#include <iostream>
#include <string>
#include <sstream>
#include "predecoding_graph_helper.h"
#include "syndrome_surgery.h"
#include <sys/stat.h>
#include <dirent.h>
#include <sys/types.h>
// #include <dirent.h>

using namespace std;

const uint64_t arr_N = 300;

void error_chain_distr(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er);

void error_chain_distr_v2(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, FILE *file = NULL);

void degree_one_distr(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er);

// This function calculates the average of remained HWs for each original HW
void HW_remained_distr(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, uint ec_length = 1);

void HW_distr(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er);

// This function calculates the frequnecy of HWs after prematching regardless of their original HW
void HW_remained_distr_v2(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, uint ec_length = 1);

void decoding_cost_experiments(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er); 

std::string generate_experiment_name(std::string exp_type, uint distance, fp_t physical_er, fp_t m_error, uint64_t total_shots = 5'000'000, bool mention_shots = true);

std::string generate_matching_file_name(std::string exp_type, uint distance, 
fp_t physical_er, fp_t m_error, uint64_t total_shots, bool mention_shots);

std::string double_to_scientific_notation_string(double number);

std::string number_of_shots_to_string(uint64_t shots_number);

void print_info(const std::vector<uint8_t>& syndrome,
                qpd::PredecoderShotResult predecoder_result, uint n_detectors, 
                bool print_syndrome = false, bool check_syndrome_correctness = false);

void check_subgraphs(uint distance, fp_t physcial_error, fp_t meas_er, uint ec_length = 1, bool print_l1_prematches = false);

void predecoder_plus_MWPM_decoder(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er);

//void measurements_stacked(uint distance, fp_t physcial_error, fp_t meas_er, uint ec_length = 1, bool print_l1_prematches = false);

//void high_distance_experiments(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, FILE *file = NULL);

void print_decoding_results(const std::vector<uint8_t>& syndrome, qrc::DecoderShotResult results, uint n_detectors);

void high_distance_experiments(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er);

void print_decoding_results_v2(const std::vector<uint8_t>& syndrome, qrc::DecoderShotResult results,
 uint n_detectors, std::map<uint, uint> pre_matches, qrc::DecoderShotResult ppost_results, std::ofstream& out);

void high_distance_predecoder_plus_MWPM_decoder(std::ofstream& out,uint64_t total_shots, 
uint distance, fp_t physcial_error, fp_t meas_er, std::ofstream& out1);

void saving_hhw_syndromes(std::ofstream& out, uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, std::string syndrome_folder_name, bool print_logs = false);

std::string syndrome_file_name(std::string experiment_name, uint distance, fp_t physical_er, fp_t m_error, uint64_t total_shots, uint8_t core_number ,uint64_t round, bool partially = false);

void predecoding_hhw_syndromes(std::ofstream& out, std::string syndromes_directory, uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, bool print_logs = false, bool adding_boundry = true);

void predecoding_hhw_syndromes_astreaG(std::ofstream& out, std::string syndromes_directory, uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, bool print_logs = false, bool adding_boundry = true);

uint64_t parseMagnitude(const std::string& input);


#endif