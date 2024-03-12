// Copyright 2021 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "stim/simulators/detection_simulator.h"
#include <random>

#include "stim/simulators/frame_simulator.h"

using namespace stim;

static std::mt19937_64 RNG(0);

template <typename T>
void xor_measurement_set_into_result(
    const T &measurement_set, simd_bit_table &frame_samples, simd_bit_table &output, size_t output_index_ticker,
    simd_bit_table &leak_samples, bool use_leak_samples=false) 
{
    simd_bits_range_ref dst = output[output_index_ticker];
    for (auto i : measurement_set) {
        dst.for_each_word(
                frame_samples[i], 
                leak_samples[i], 
                [&](simd_word& d, simd_word& f, simd_word& l) {
                    const simd_word rand(RNG(), RNG());
                    if (use_leak_samples) {
                        d ^= l.andnot(f) | (rand & l);  // Leakage is always random.
                    } else {
                        d ^= l | f;
                    }
                });
    }
}

template <typename T>
void or_measurement_set_into_result(
    const T &measurement_set, 
    simd_bit_table &leak_samples,
    simd_bit_table &output,
    size_t output_index_ticker)
{
    simd_bits_range_ref dst = output[output_index_ticker];
    for (auto i : measurement_set) {
        dst |= leak_samples[i];
    }
}

void stim::read_from_sim(
    FrameSimulator &sim,
    const DetectorsAndObservables &det_obs,
    bool prepend_observables,
    bool append_observables,
    bool get_leakage_data,
    simd_bit_table &result_table,
    simd_bit_table &leakage_table,
    uint64_t max_detector)
{
    auto frame_samples = sim.m_record.storage;
    auto leakage_samples = sim.leak_record.storage;

    // Xor together measurement samples to form detector samples.
    size_t offset = 0;
    if (prepend_observables) {
        for (const auto &obs : det_obs.observables) {
            xor_measurement_set_into_result(obs, frame_samples, result_table, offset, leakage_samples, true);
            if (get_leakage_data) {
                or_measurement_set_into_result(obs, leakage_samples,
                        leakage_table, offset);
            }
            offset++;
        }
    }
    for (const auto &det : det_obs.detectors) {
        xor_measurement_set_into_result(det, frame_samples, result_table, offset, leakage_samples, false);
        if (get_leakage_data) {
            or_measurement_set_into_result(det, leakage_samples,
                    leakage_table, offset);
        }
        offset++;
        if (offset >= max_detector) {
            break;
        }
    }
    if (append_observables) {
        for (const auto &obs : det_obs.observables) {
            xor_measurement_set_into_result(obs, frame_samples, result_table, offset, leakage_samples, true);
            if (get_leakage_data) {
                or_measurement_set_into_result(obs, leakage_samples,
                        leakage_table, offset);
            }
            offset++;
        }
    }
}


simd_bit_table stim::detector_samples(
    const Circuit &circuit,
    const DetectorsAndObservables &det_obs,
    size_t num_shots,
    bool prepend_observables,
    bool append_observables,
    std::mt19937_64 &rng)
{
    simd_bit_table false_table(1, 1);
    return detector_samples(circuit, det_obs, num_shots, prepend_observables, append_observables, rng, false, false_table);
}


simd_bit_table stim::detector_samples(
    const Circuit &circuit,
    const DetectorsAndObservables &det_obs,
    size_t num_shots,
    bool prepend_observables,
    bool append_observables,
    std::mt19937_64 &rng,
    bool get_leakage_data,
    simd_bit_table& leakage_table) 
{
    // Start from measurement samples.
    FrameSimulator sim(circuit.count_qubits(), num_shots, SIZE_MAX, rng);
    sim.reset_all_and_run(circuit);

    auto num_detectors = det_obs.detectors.size();
    auto num_obs = det_obs.observables.size();
    size_t num_results = num_detectors + num_obs * (prepend_observables + append_observables);
    simd_bit_table result_table(num_results, num_shots);
    if (get_leakage_data) {
        leakage_table = simd_bit_table(num_results, num_shots);
    }

    read_from_sim(sim, det_obs, prepend_observables, append_observables,
            get_leakage_data, result_table, leakage_table);

    return result_table;
}

simd_bit_table stim::detector_samples(
    const Circuit &circuit,
    size_t num_shots,
    bool prepend_observables,
    bool append_observables,
    std::mt19937_64 &rng)
{
    simd_bit_table false_table(1, 1);
    return detector_samples(circuit, num_shots, prepend_observables, 
            append_observables, rng, false, false_table);
}

simd_bit_table stim::detector_samples(
    const Circuit &circuit,
    size_t num_shots,
    bool prepend_observables,
    bool append_observables,
    std::mt19937_64 &rng,
    bool get_leakage_data,
    simd_bit_table& leakage_table) 
{
    return detector_samples(
        circuit, DetectorsAndObservables(circuit), num_shots, prepend_observables, append_observables, rng, get_leakage_data, leakage_table);
}

void detector_sample_out_helper_stream(
    const Circuit &circuit,
    FrameSimulator &sim,
    size_t num_samples,
    bool append_observables,
    FILE *out,
    SampleFormat format) {
    MeasureRecordBatchWriter writer(out, num_samples, format);
    std::vector<simd_bits> observables;
    sim.reset_all();
    writer.begin_result_type('D');
    simd_bit_table detector_buffer(1024, num_samples);
    size_t buffered_detectors = 0;
    circuit.for_each_operation([&](const Operation &op) {
        if (op.gate->id == gate_name_to_id("DETECTOR")) {
            simd_bits_range_ref result = detector_buffer[buffered_detectors];
            result.clear();
            for (auto t : op.target_data.targets) {
                assert(t.data & TARGET_RECORD_BIT);
                result ^= sim.m_record.lookback(t.data ^ TARGET_RECORD_BIT);
            }
            buffered_detectors++;
            if (buffered_detectors == 1024) {
                writer.batch_write_bytes(detector_buffer, 1024 >> 6);
                buffered_detectors = 0;
            }
        } else if (op.gate->id == gate_name_to_id("OBSERVABLE_INCLUDE")) {
            if (append_observables) {
                size_t id = (size_t)op.target_data.args[0];
                while (observables.size() <= id) {
                    observables.emplace_back(num_samples);
                }
                simd_bits_range_ref result = observables[id];

                for (auto t : op.target_data.targets) {
                    assert(t.data & TARGET_RECORD_BIT);
                    result ^= sim.m_record.lookback(t.data ^ TARGET_RECORD_BIT);
                }
            }
        } else {
            (sim.*op.gate->frame_simulator_function)(op.target_data);
            sim.m_record.mark_all_as_written();
        }
    });
    for (size_t k = 0; k < buffered_detectors; k++) {
        writer.batch_write_bit(detector_buffer[k]);
    }
    writer.begin_result_type('L');
    for (const auto &result : observables) {
        writer.batch_write_bit(result);
    }
    writer.write_end();
}

void detector_samples_out_in_memory(
    const Circuit &circuit,
    size_t num_shots,
    bool prepend_observables,
    bool append_observables,
    FILE *out,
    SampleFormat format,
    std::mt19937_64 &rng,
    FILE *obs_out,
    SampleFormat obs_out_format) {
    bool obs_in_det_output = prepend_observables || append_observables;
    if (prepend_observables + append_observables + (obs_out != nullptr) > 1) {
        throw std::out_of_range("Can't combine --prepend_observables, --append_observables, or --obs_out");
    }

    DetectorsAndObservables det_obs(circuit);
    size_t no = det_obs.observables.size();
    size_t nd = det_obs.detectors.size();

    char c1, c2;
    size_t ct;
    if (prepend_observables) {
        c1 = 'L';
        c2 = 'D';
        ct = no;
    } else if (append_observables) {
        c1 = 'D';
        c2 = 'L';
        ct = nd;
    } else {
        c1 = 'D';
        c2 = 'D';
        ct = 0;
    }

    auto table = detector_samples(
        circuit, det_obs, num_shots, prepend_observables, append_observables || obs_out != nullptr, rng);

    if (obs_out != nullptr) {
        simd_bit_table obs_data(no, num_shots);
        for (size_t k = 0; k < no; k++) {
            obs_data[k] = table[nd + k];
            table[nd + k].clear();
        }
        write_table_data(obs_out, num_shots, no, simd_bits(0), obs_data, obs_out_format, 'L', 'L', no);
    }

    write_table_data(out, num_shots, nd + no * obs_in_det_output, simd_bits(0), table, format, c1, c2, ct);
}

void detector_sample_out_helper(
    const Circuit &circuit,
    FrameSimulator &sim,
    size_t num_shots,
    bool prepend_observables,
    bool append_observables,
    FILE *out,
    SampleFormat format,
    std::mt19937_64 &rng,
    FILE *obs_out,
    SampleFormat obs_out_format) {
    uint64_t d = circuit.count_detectors() + circuit.count_observables();
    uint64_t approx_mem_usage = std::max(num_shots, size_t{256}) * std::max(circuit.count_measurements(), d);
    if (!prepend_observables && obs_out == nullptr && should_use_streaming_instead_of_memory(approx_mem_usage)) {
        detector_sample_out_helper_stream(circuit, sim, num_shots, append_observables, out, format);
    } else {
        detector_samples_out_in_memory(
            circuit, num_shots, prepend_observables, append_observables, out, format, rng, obs_out, obs_out_format);
    }
}

void stim::detector_samples_out(
    const Circuit &circuit,
    size_t num_shots,
    bool prepend_observables,
    bool append_observables,
    FILE *out,
    SampleFormat format,
    std::mt19937_64 &rng,
    FILE *obs_out,
    SampleFormat obs_out_format) {
    constexpr size_t GOOD_BLOCK_SIZE = 768;
    size_t num_qubits = circuit.count_qubits();
    size_t max_lookback = circuit.max_lookback();
    if (num_shots >= GOOD_BLOCK_SIZE) {
        auto sim = FrameSimulator(num_qubits, GOOD_BLOCK_SIZE, max_lookback, rng);
        while (num_shots > GOOD_BLOCK_SIZE) {
            detector_sample_out_helper(
                circuit,
                sim,
                GOOD_BLOCK_SIZE,
                prepend_observables,
                append_observables,
                out,
                format,
                rng,
                obs_out,
                obs_out_format);
            num_shots -= GOOD_BLOCK_SIZE;
        }
    }
    if (num_shots) {
        auto sim = FrameSimulator(num_qubits, num_shots, max_lookback, rng);
        detector_sample_out_helper(
            circuit,
            sim,
            num_shots,
            prepend_observables,
            append_observables,
            out,
            format,
            rng,
            obs_out,
            obs_out_format);
    }
}
