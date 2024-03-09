/*
 *  author: Narges Alavisamani
 *  date:   28 March 2023
 * */


#include <stim.h>
#include <mwpm_decoder.h>
#include <benchmark.h>
#include <defs.h>
#include <mpi.h>
#include <math.h>
#include <random>
#include <cstdint>

#ifndef SYNDROME_SURGERY_H
#define SYNDROME_SURGERY_H

#include <algorithm>
#include <fstream>
#include <iostream>
#include "decoder.h"
#include "decoding_graph.h"
#include "predecoding_graph_helper.h"
#include "predecoder.h"

namespace qpd{

void print_graph(qrc::DecodingGraph decoding_graph);

void print_graph_and_syndrome(qrc::DecodingGraph graph, const std::vector<uint8_t>& syndrome);

void testing_syndroms_and_graphs(uint64_t total_shots, uint distance, 
                fp_t physcial_error, fp_t meas_er);

std::vector<uint16_t> syndrome_compressed(const std::vector<uint8_t>& syndrome, bool print_syndrome = true);

void write_vector(const std::vector<uint16_t>& vec, const std::string& filename);

std::vector<uint16_t> read_vector(const std::string& filename);

void write_vector_of_vectors(const std::vector<std::vector<uint16_t>>& data, const std::string& filename);

std::vector<std::vector<uint16_t>> read_vector_of_vectors(const std::string& filename);

std::vector<uint8_t> syndrome_decompressed(std::vector<uint16_t> synd_vec, uint d);

void print_syndrome(const std::vector<uint8_t>& syndrome);

}//namespace qpd

#endif