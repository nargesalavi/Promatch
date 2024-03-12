/*
 *  Modified by Suhas Vittal on 23 August 2022
 * */

#include "stim/gen/circuit_gen_params.h"
#include <random>
#include <map>

#include "stim/arg_parse.h"

using namespace stim;

const auto seed = std::chrono::system_clock::now().time_since_epoch().count();
static std::mt19937_64 CIRCGEN_RNG(seed);

#define K(m,s)  ( ((m)*(m)) / ((s)*(s)) )
#define T(m,s)  ( ((s)*(s)) / (m) )

static std::map<uint32_t, double> rounddp_table;
static std::map<uint32_t, double> unitary1dp_table;
static std::map<std::vector<uint32_t>, double> unitary2dp_table;
static std::map<uint32_t, double> premeasflip_table;
static std::map<uint32_t, double> postresflip_table;

static std::map<uint32_t, double> roundleak_table;
static std::map<std::vector<uint32_t>, double> cliffordleak_table;

static std::map<std::vector<uint32_t>, double> cliffordtransport_table;

void
CircuitGenParameters::reset_data() const {
    rounddp_table.clear();
    unitary1dp_table.clear();
    unitary2dp_table.clear();
    premeasflip_table.clear();
    postresflip_table.clear();
    roundleak_table.clear();
    cliffordleak_table.clear();
    cliffordtransport_table.clear();
}

template <typename K> double
get_from(std::map<K, double>& table, K k, double store_on_fail) {
    double p;
    if (table.count(k)) {
        p = table[k];
    } else {
        p = store_on_fail;
        table[k] = p;
    }
    return p;
}


void append_anti_basis_error(Circuit &circuit, 
        const std::vector<uint32_t> &targets, double p, char basis) 
{
    if (p > 0) {
        if (basis == 'X') {
            circuit.append_op("Z_ERROR", targets, p);
        } else {
            circuit.append_op("X_ERROR", targets, p);
        }
    }
}

void 
CircuitGenParameters::validate_params() const {
    // TODO: Update for standard deviation.
    if (before_measure_flip_probability < 0 || before_measure_flip_probability > 1) {
        throw std::invalid_argument("not 0 <= before_measure_flip_probability <= 1");
    }
    if (before_round_data_depolarization < 0 || before_round_data_depolarization > 1) {
        throw std::invalid_argument("not 0 <= before_round_data_depolarization <= 1");
    }
    if (after_clifford_depolarization < 0 || after_clifford_depolarization > 1) {
        throw std::invalid_argument("not 0 <= after_clifford_depolarization <= 1");
    }
    if (after_reset_flip_probability < 0 || after_reset_flip_probability > 1) {
        throw std::invalid_argument("not 0 <= after_reset_flip_probability <= 1");
    }
}

CircuitGenParameters::CircuitGenParameters(
        uint64_t rounds, uint32_t distance, std::string task)
    : rounds(rounds), distance(distance), task(task) 
{}

void
CircuitGenParameters::append_begin_round_tick(
        Circuit &circuit, const std::vector<uint32_t> &data_qubits) 
const {
    circuit.append_op("TICK", {});
    for (uint32_t d : data_qubits) {
        std::vector<uint32_t> singleton{d};
        double p = get_from(rounddp_table, d, 
                get_before_round_data_depolarization());
        if (p > 0) {
            circuit.append_op("DEPOLARIZE1", singleton, p);
        }
    }
    for (uint32_t d : data_qubits) {
        std::vector<uint32_t> singleton{d};
        double p = get_from(roundleak_table, d, 
                    get_before_round_leakage_probability());
        if (p > 0) {
            circuit.append_op("L_ERROR", singleton, p);
        }
    }
}

void 
CircuitGenParameters::append_unitary_1(
    Circuit &circuit, const std::string &name, const std::vector<uint32_t> targets)
const {
    circuit.append_op(name, targets);
    for (uint32_t t : targets) {
        std::vector<uint32_t> singleton{t};
        double p = get_from(unitary1dp_table, t, get_after_clifford_depolarization(true));
        if (p > 0) {
            circuit.append_op("DEPOLARIZE1", singleton, p);
        }
    }
}

void
CircuitGenParameters::append_unitary_2(
    Circuit &circuit, const std::string &name, const std::vector<uint32_t> targets)
const {
    circuit.append_op(name, targets);
    for (uint32_t i = 0; i < targets.size(); i += 2) {
        std::vector<uint32_t> cx{targets[i], targets[i+1]};
        double p = get_from(unitary2dp_table, cx, 
                            get_after_clifford_depolarization());
        if (p > 0) {
            circuit.append_op("DEPOLARIZE2", cx, p);
        }
    }
    // I know this is "less optimal", but it is easier to read.
    for (uint32_t i = 0; i < targets.size(); i += 2) {
        std::vector<uint32_t> pair{targets[i], targets[i+1]};
        double p = get_from(cliffordleak_table, pair, get_after_clifford_leakage_probability());
        if (p > 0) {
            circuit.append_op("L_ERROR", pair, p);
        }
    }

    for (uint32_t i = 0; i < targets.size(); i += 2) {
        std::vector<uint32_t> pair{targets[i], targets[i+1]};
        double p = get_from(cliffordtransport_table, pair, get_after_clifford_leakage_transport());
        if (p > 0) {
            circuit.append_op("L_TRANSPORT", pair, p);
        }
    }
}

void
CircuitGenParameters::append_reset(
        Circuit &circuit, const std::vector<uint32_t> targets, char basis)
const {
    circuit.append_op(std::string("R") + basis, targets);
    for (uint32_t t : targets) {
        std::vector<uint32_t> singleton{t};
        double p = get_from(postresflip_table, t, 
                            get_after_reset_flip_probability());
        append_anti_basis_error(circuit, singleton, p, basis);
    }
}

void 
CircuitGenParameters::append_measure(
        Circuit &circuit, const std::vector<uint32_t> targets, char basis)
const {
    for (uint32_t t : targets) {
        std::vector<uint32_t> singleton{t};
        double p = get_from(premeasflip_table, t, 
                            get_before_measure_flip_probability());
        append_anti_basis_error(circuit, singleton, p, basis);
    }
    circuit.append_op(std::string("M") + basis, targets);
}

void 
CircuitGenParameters::append_measure_reset(
    Circuit &circuit, const std::vector<uint32_t> targets, char basis)
const {
    for (uint32_t t : targets) {
        std::vector<uint32_t> singleton{t};
        double p = get_from(premeasflip_table, t,
                            get_before_measure_flip_probability());
        append_anti_basis_error(circuit, singleton, p, basis);
    }
    circuit.append_op(std::string("MR") + basis, targets);
    for (uint32_t t : targets) {
        std::vector<uint32_t> singleton{t};
        double p = get_from(postresflip_table, t, 
                            get_after_reset_flip_probability());
        append_anti_basis_error(circuit, singleton, p, basis);
    }
}

double
CircuitGenParameters::get_after_clifford_depolarization(bool single_qubit_gate) const {
    if (single_qubit_gate && after_clifford_sq_depolarization >= 0) {
        return get_error(after_clifford_sq_depolarization, after_clifford_sq_depolarization_stddev);
    } else {
        return get_error(after_clifford_depolarization, after_clifford_depolarization_stddev);
    }
}

double
CircuitGenParameters::get_before_round_data_depolarization() const {
    return get_error(before_round_data_depolarization, before_round_data_depolarization_stddev);
}

double
CircuitGenParameters::get_before_measure_flip_probability() const {
    return get_error(before_measure_flip_probability, before_measure_flip_probability_stddev);
}

double
CircuitGenParameters::get_after_reset_flip_probability() const {
    return get_error(after_reset_flip_probability, after_reset_flip_probability_stddev);
}

double
CircuitGenParameters::get_before_round_leakage_probability() const {
    return get_error(before_round_leakage_probability, before_round_leakage_probability_stddev);
}

double
CircuitGenParameters::get_after_clifford_leakage_probability() const {
    return get_error(after_clifford_leakage_probability, after_clifford_leakage_probability_stddev);
}

double
CircuitGenParameters::get_after_clifford_leakage_transport() const {
    return get_error(after_clifford_leakage_transport, after_clifford_leakage_transport_stddev);
}

double
CircuitGenParameters::get_error(double mean, double stddev) const {
    if (stddev == 0.0) {
        return mean;
    }
    double e;
    if (dist == Distribution::normal) {
        std::normal_distribution<double> dist{mean, stddev};
        e = dist(CIRCGEN_RNG);
        if (e < after_reset_flip_probability) {
            e = 2*after_reset_flip_probability - e;
        }
    } else if (dist == Distribution::gamma) {
        double k = K(mean, stddev);
        double t = T(mean, stddev);
        std::gamma_distribution<double> dist{k, t};
        e = dist(CIRCGEN_RNG) + mean;
    } else {
        double m = log(mean);
        // Assume stddev is for the logarithm.
        std::lognormal_distribution<double> dist{m, stddev};
        e = dist(CIRCGEN_RNG);
    }
    return e;
}

std::string GeneratedCircuit::layout_str() const {
    std::stringstream ss;
    std::vector<std::vector<std::string>> lines;
    for (const auto &kv : layout) {
        auto x = kv.first.first;
        auto y = kv.first.second;
        while (lines.size() <= y) {
            lines.push_back({});
        }
        while (lines[y].size() <= x) {
            lines[y].push_back("");
        }
        lines[y][x] = kv.second.first + std::to_string(kv.second.second);
    }
    size_t max_len = 0;
    for (const auto &line : lines) {
        for (const auto &entry : line) {
            max_len = std::max(max_len, entry.size());
        }
    }
    for (auto p = lines.crbegin(); p != lines.crend(); p++) {
        const auto &line = *p;
        ss << "#";
        for (const auto &entry : line) {
            ss << ' ' << entry;
            ss << std::string(max_len - entry.size(), ' ');
        }
        ss << "\n";
    }
    return ss.str();
}
