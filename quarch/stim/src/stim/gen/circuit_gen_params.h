/*
 *  Modified by Suhas Vittal on 23 August 2022
 * */

#ifndef _STIM_GEN_CIRCUIT_GEN_PARAMS_H
#define _STIM_GEN_CIRCUIT_GEN_PARAMS_H

#include <chrono>
#include <map>
#include <random>

#include <math.h>
#include <stddef.h>
#include <stdint.h>

#include "stim/circuit/circuit.h"

namespace stim {
class CircuitGenParameters {
public:
    uint64_t rounds;
    uint32_t distance;
    std::string task;

    enum Distribution {normal, gamma, lognormal};

    // If the stddev are nonzero, then
    // these values are the mean.
    // Note: these haven't been renamed to avoid
    // breaking code.
    double after_clifford_depolarization = 0;
    double before_round_data_depolarization = 0;
    double before_measure_flip_probability = 0;
    double after_reset_flip_probability = 0;

    double after_clifford_sq_depolarization = -1;  // If negative, then will use after_clifford_depolarization.
    
    double before_round_leakage_probability = 0;
    double after_clifford_leakage_probability = 0;

    double after_clifford_leakage_transport = 0;

    // These will be the standard deviation.
    double after_clifford_depolarization_stddev = 0;
    double before_round_data_depolarization_stddev = 0;
    double before_measure_flip_probability_stddev = 0;
    double after_reset_flip_probability_stddev = 0;

    double after_clifford_sq_depolarization_stddev = 0;

    double before_round_leakage_probability_stddev = 0;
    double after_clifford_leakage_probability_stddev = 0;

    double after_clifford_leakage_transport_stddev = 0;

    Distribution dist = Distribution::lognormal;
    
    bool both_stabilizers = false;
    bool swap_lru = false;
    bool swap_lru_with_no_swap = false;

    bool initial_state_is_basis_1 = false;

    void validate_params() const;

    CircuitGenParameters(uint64_t rounds, uint32_t distance, std::string task);
    void append_begin_round_tick(
            Circuit &circuit, const std::vector<uint32_t> &data_qubits) const;
    void append_unitary_1(
            Circuit &circuit, const std::string &name, 
            const std::vector<uint32_t> targets) const;
    void append_unitary_2(
            Circuit &circuit, const std::string &name, 
            const std::vector<uint32_t> targets) const;
    void append_reset(
            Circuit &circuit, const std::vector<uint32_t> targets,
            char basis = 'Z') const;
    void append_measure(
            Circuit &circuit, const std::vector<uint32_t> targets,
            char basis = 'Z') const;
    void append_measure_reset(
            Circuit &circuit, const std::vector<uint32_t> targets,
            char basis = 'Z') const;
    void reset_data(void) const;

    double get_after_clifford_depolarization(bool single_qubit_gate=false) const;
    double get_before_round_data_depolarization(void) const;
    double get_before_measure_flip_probability(void) const;
    double get_after_reset_flip_probability(void) const;
    double get_before_round_leakage_probability(void) const;
    double get_after_clifford_leakage_probability(void) const;
    double get_after_clifford_leakage_transport(void) const;
private:
    double get_error(double mean, double stddev) const;

    std::mt19937_64 rng; 
};

struct GeneratedCircuit {
    Circuit circuit;
    std::map<std::pair<uint32_t, uint32_t>, std::pair<std::string, uint32_t>> layout;
    std::string hint_str;
    std::string layout_str() const;
};
}  // namespace stim

#endif
