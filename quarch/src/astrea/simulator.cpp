/*
 *  author: Suhas Vittal
 *  date:   7 September 2022
 * */

#include "astrea/simulator.h"
#include <limits>

namespace qrc {
namespace astrea {

static const uint BFU_SORT_STAGES = 4;
static const uint RADIX_WIDTH = 16 / (BFU_SORT_STAGES);

AstreaSimulator::AstreaSimulator(DecodingGraph& graph, const stim::Circuit& circ, const AstreaSimulatorParams& params)
    :
    /* Statistics */
    prefetch_cycles(0),
    bfu_cycles(0),
    cycles_to_converge(0),
    valid_weights_after_filter(0),
    /* Microarchitecture */
    register_file(params.n_registers, (Register){std::make_pair(0,0), 0, false}),
    detector_vector_register(),
    mean_weight_register(0),
    access_counter(0),
    major_detector_register(0),
    minor_detector_table(),
    hardware_deques(params.bfu_fetch_width, 
            AstreaDeque(params.bfu_priority_queue_size)),
    bfu_pipeline_latches(params.bfu_fetch_width, 
            std::vector<BFUPipelineLatch>(
                BFU_SORT_STAGES+params.bfu_compute_stages)),
    replacement_queue(),
    best_matching_register(),
    mld_bm_register(),
    mld_prob_register(),
    state(AstreaSimulator::State::prefetch),
    bfu_idle(false),
    /* Data */
    graph(graph),
    path_table(compute_path_table(graph)),
    /* Config parameters */
    curr_max_detector(0),
    weight_filter_cutoff(params.weight_filter_cutoff),
    n_detectors(params.n_detectors),
    n_detectors_per_round(params.n_detectors_per_round),
    bfu_fetch_width(params.bfu_fetch_width),
    bfu_compute_stages(params.bfu_compute_stages),
    curr_qubit(0),
    has_boundary(false),
    use_mld(params.use_mld),
    cycles_after_last_converge(0),
    hw6decoder(nullptr)
{
    clear();
    if (use_mld) {
        hw6decoder = new MLDDecoder(circ);
    } else {
        hw6decoder = new MWPMDecoder(circ);
    }
}

AstreaSimulator::~AstreaSimulator() {
    delete hw6decoder;
}

void
AstreaSimulator::load_detectors(const std::vector<uint>& detector_array) {
    if (detector_array.size() == 0) {
        state = State::idle;
    }
    clear();
    detector_vector_register = detector_array;
#ifdef ASTREA_DEBUG
    std::cout << "LOAD:";
    for (uint d : detector_array) {
        std::cout << " " << d; 
    }
    std::cout << "\n";
#endif
    has_boundary = detector_vector_register.back() == BOUNDARY_INDEX;
    // Compute the critical index.
    curr_max_detector = n_detectors_per_round;
}

void
AstreaSimulator::tick() {
    // Update last use for every register.
    for (Register& r : register_file) {
        r.last_use++;
    }
    // Go perform any replacements if necessary.
    if (replacement_queue.size()) {
        addr_t address = replacement_queue.front();
        // Check if there are any evictable registers.
        Register * victim = &register_file[0];
        for (uint i = 1; i < register_file.size(); i++) {
            Register& r = register_file[i];
            if (victim->valid &&
                    (!r.valid || r.last_use > victim->last_use)) 
            {
                victim = &register_file[i];
            }
        }
        victim->address = address;
        victim->valid = true;
        victim->last_use = 0;

        replacement_queue.pop_front();
    }

    switch (state) {
    case State::prefetch:
        tick_prefetch();
        break;
    case State::bfu:
#ifdef ASTREA_DEBUG
        if (bfu_cycles % 100 == 0) {
            for (uint i = 0; i < hardware_deques.size(); i++) {
                std::cout << "[PQ] Priority Queue " << i << ":\n";
                for (auto entry : hardware_deques[i].backing_array) {
                    std::cout << "\t";
                    for (auto pair : entry.running_matching) {
                        std::cout << " " << pair.first << " --> " << pair.second;
                    }
                    std::cout << "\n";
                }
            }
        }
#endif
        tick_bfu();
        break;
    default:
        break;
    }
}

void
AstreaSimulator::sig_end_round(uint rounds_ended) {
    // We subtract by 1 because we don't count the boundary.
    curr_max_detector += n_detectors_per_round * rounds_ended;
    if (curr_max_detector > n_detectors-1) {
        curr_max_detector = n_detectors - 1;
    }
    major_detector_register = 0;
#ifdef ASTREA_DEBUG
    std::cout << "SIGNAL -- ROUND END | Max detector = " 
        << curr_max_detector << "\n";
#endif
}

bool
AstreaSimulator::is_idle() {
    return state == State::idle;
}

void
AstreaSimulator::force_idle() {
    state = State::idle;
}

std::map<uint, uint>
AstreaSimulator::get_matching() {
    if (use_mld) {
        if (mld_prob_register[0] > mld_prob_register[1]) {
            return mld_bm_register[0].running_matching;
        } else {
            return mld_bm_register[1].running_matching;
        }
    } else {
        return best_matching_register.running_matching;
    }
}

void
AstreaSimulator::reset_stats() {
    prefetch_cycles = 0;
    bfu_cycles = 0;
    cycles_to_converge = 0;
    valid_weights_after_filter = 0;
    cycles_after_last_converge = 0;
}

void
AstreaSimulator::tick_prefetch()  {
    uint& ai = major_detector_register;
    if (minor_detector_table.count(ai) == 0) {
        // Start with boundary first if it exists.
        if (has_boundary) {
            minor_detector_table[ai] = detector_vector_register.size() - 1;
        } else {
            minor_detector_table[ai] = ai + 1;
        }
    }
    uint& aj = minor_detector_table[ai];
    
    while (aj >= detector_vector_register.size()
        && ai < detector_vector_register.size()-1)
    {
        ai++;
        aj = minor_detector_table[ai];
    }

    if (ai >= detector_vector_register.size()-1) {
        // We only transition state if all DRAM requests
        // have been processed.
        update_state();
        return;
    }
    uint di = detector_vector_register[ai];
    uint dj = detector_vector_register[aj];

    if (di != BOUNDARY_INDEX && di > curr_max_detector) {
        return;
    }

    if (dj != BOUNDARY_INDEX && dj > curr_max_detector) {
        uint next_di = detector_vector_register[ai+1];
        if (next_di == BOUNDARY_INDEX || next_di < curr_max_detector) {
            ai++;
        }
        return;
    }

    addr_t address = std::make_pair(di, dj);
    
    access(address, false);  // We don't care if there is a hit or miss.
    
    auto vdi = graph.get_vertex(di);
    auto vdj = graph.get_vertex(dj);

    fp_t raw_weight = path_table[std::make_pair(vdi, vdj)].distance;
    uint32_t w = (uint32_t) (raw_weight * MWPM_INTEGER_SCALE);
    mean_weight_register += w;
    access_counter++;
    if (w <= weight_filter_cutoff) {
        valid_weights_after_filter++;
    }
    // Update ai, aj.
    // Update varies on whether or not boundary is in the DVR.
    if (has_boundary && dj == BOUNDARY_INDEX) {
        aj = ai + 1;  // Reset, we have finished the boundary.
    } else {
        aj++;
    }

    if ((has_boundary && aj >= detector_vector_register.size()-1) 
            || (!has_boundary && aj >= detector_vector_register.size())) 
    {
        aj = detector_vector_register.size();
        ai++;
    }

    prefetch_cycles++;
}

void
AstreaSimulator::tick_bfu() {
    // We simulate the stages of the pipeline in reverse order.
    // There is one FETCH stage, four SORT stages, and 
    // a variable number of COMPUTE stages.
    bfu_idle = true;
    for (uint stage = bfu_compute_stages; stage > 0; stage--) {
        tick_bfu_compute(stage-1);
    }
    for (uint stage = BFU_SORT_STAGES; stage > 0; stage--) {
        tick_bfu_sort(stage-1);
    }
    tick_bfu_fetch();
    if (bfu_idle) {
        update_state();
    }

    bfu_cycles++;
    cycles_after_last_converge++;
}

void
AstreaSimulator::tick_bfu_compute(uint stage) {
    const uint STAGE_OFFSET = BFU_SORT_STAGES;
    // Get data from the pipeline latch.
    for (uint f = 0; f < bfu_fetch_width; f++) {
        BFUPipelineLatch& latch = bfu_pipeline_latches[f][stage+STAGE_OFFSET];
        if (!latch.valid) {
            continue;
        }
#ifdef ASTREA_DEBUG2
        if (stage == 0 && f == 0) {
            // Check if proposed_matches are sorted.
            std::cout << "\tProposed matches (COMPUTE 0):";
            for (uint i = 0; i < latch.proposed_matches.size(); i++) {
                auto pair = latch.proposed_matches[i];
                std::cout << " " << pair.first << "( w = " << pair.second << " )";
            }
            std::cout << "\n";
        }
#endif
        bfu_idle = false;

        // Add entry from proposed matches to the running matching.
        // Push a matching onto each hardware deque.
        uint offset = f;
        for (uint i = 0; i < hardware_deques.size(); i++) {
            uint deque_index = (i + f) % hardware_deques.size();
            auto& pq = hardware_deques[deque_index];
            if (latch.proposed_matches.empty()) {
                break;
            }
            std::map<uint, uint> matching(latch.base_entry.running_matching);
            uint32_t matching_weight = latch.base_entry.matching_weight;
            auto match = latch.proposed_matches.front();
            uint di = detector_vector_register[
                        latch.base_entry.next_unmatched_index];
            uint dj = match.first;
            matching[di] = dj;
            matching[dj] = di;
            // Update weight
            matching_weight += match.second;
            // Compute next unmatched detector.
            uint next_unmatched_index = detector_vector_register.size();
            for (uint ai = latch.base_entry.next_unmatched_index + 1; 
                    ai < detector_vector_register.size(); ai++) 
            {
                if (matching.count(detector_vector_register[ai]) == 0) {
                    next_unmatched_index = ai;
                    break;
                } 
            }
            // Push this result onto the stack if there is a next unmatched detector
            // and hamming weight > 6.
            DequeEntry ei = {matching, matching_weight, next_unmatched_index};
            uint hw = detector_vector_register.size() - matching.size();
            if (next_unmatched_index < detector_vector_register.size() && hw > 6) {
                pq.push(ei);
            } else {
                // Match the rest using MWPM.
                std::vector<uint8_t> subsyndrome(n_detectors-1, 0);
                for (uint d : detector_vector_register) {
                    if (d == BOUNDARY_INDEX) {
                        continue;
                    }
                    if (!matching.count(d)) {
                        subsyndrome[d] = 1;
                    }
                }
                fp_t base_prob = pow(10, -((fp_t)matching_weight) / ((fp_t)MWPM_INTEGER_SCALE));
                DecoderShotResult hw6res = hw6decoder->decode_error(subsyndrome);
                if (use_mld) {
                    for (uint8_t obs = 0; obs <= 1; obs++) {
                        std::map<uint, uint> obs_matching(matching);
                        auto rep = ((MLDDecoder*)hw6decoder)->get_representative(obs);
                        for (auto pair : rep) {
                            if (obs_matching.count(pair.first) || obs_matching.count(pair.second)) {
                                continue;
                            }
                            obs_matching[pair.first] = pair.second;
                            obs_matching[pair.second] = pair.first;
                        }
                        uint8_t true_obs = hw6decoder->get_correction_from_matching(obs_matching)[0];

                        if (mld_prob_register[true_obs] == 0.0) {
                            mld_bm_register[true_obs] = {obs_matching, 0, detector_vector_register.size()};
                        }
                        fp_t p = base_prob * ((MLDDecoder*)hw6decoder)->get_obs_prob(obs);
                        mld_prob_register[true_obs] += p;
                    }
                } else {
                    auto submatching = hw6res.matching;
                    // Combine hw6 matching with existing matching.
                    for (auto pair : submatching) {
                        if (matching.count(pair.first) || matching.count(pair.second)) {
                            continue;
                        }
                        matching[pair.first] = pair.second;
                        matching[pair.second] = pair.first;
                        auto vdi = graph.get_vertex(pair.first);
                        auto vdj = graph.get_vertex(pair.second);
                        matching_weight += (uint32_t) (path_table[std::make_pair(vdi, vdj)].distance * MWPM_INTEGER_SCALE);
                    }
                    // Update the output register for MWPM.
                    if (matching_weight < best_matching_register.matching_weight) {
                        best_matching_register = {matching, matching_weight, detector_vector_register.size()};
                        // Update stats.
                        cycles_to_converge = cycles_after_last_converge;
                    } 
                }
            }
            latch.proposed_matches.pop_front();
        }
        // Update next latch if we can.
        if (stage < bfu_compute_stages-1) {
            BFUPipelineLatch& next = bfu_pipeline_latches[f][stage+STAGE_OFFSET+1];
            if (latch.proposed_matches.empty()) {
                // Invalidate future stages.
                next.valid = false;
            } else {
                next.stalled = false;
                next.valid = true;
                next.base_entry = latch.base_entry;
                next.proposed_matches = latch.proposed_matches;
            }
        }
    }
}

void
AstreaSimulator::tick_bfu_sort(uint stage) {
    // Get pipeline latch. 
    for (uint f = 0; f < bfu_fetch_width; f++) {
        BFUPipelineLatch& latch = bfu_pipeline_latches[f][stage];
        if (!latch.valid) {
            bfu_pipeline_latches[f][stage+1].valid = false;
            continue;
        }
#ifdef ASTREA_DEBUG2
        if (f == 0) {
            // Check if proposed_matches are sorted.
            std::cout << "\tProposed matches (SORT " << stage << "):";
            for (uint i = 0; i < latch.proposed_matches.size(); i++) {
                auto pair = latch.proposed_matches[i];
                std::cout << " " << pair.first << "( w = " << pair.second << " )";
            }
            std::cout << "\n";
        }
#endif

        bfu_idle = false;
        // Get radix bins from prior round.
        auto proposed_matches = latch.proposed_matches;
        uint8_t nibble_pos = stage;
        std::array<std::vector<uint>, 1<<RADIX_WIDTH> bins;
        for (uint i = 0; i < proposed_matches.size(); i++) {
            uint32_t w = proposed_matches[i].second;
            // Get nibble (4-bits) at relevant position.
            uint8_t nibble = (w >> (nibble_pos << 2)) & 0xf;
            bins[nibble].push_back(i);
        }
        // Merge into one array.
        std::deque<std::pair<uint, uint32_t>> new_matches;
        for (auto bucket : bins) {
            for (auto i : bucket) {
                new_matches.push_back(proposed_matches[i]); 
            }
        }

        BFUPipelineLatch& next = bfu_pipeline_latches[f][stage+1];
        next.proposed_matches = new_matches;
        next.base_entry = latch.base_entry;
        next.stalled = false;
        next.valid = true;
    }
}

void
AstreaSimulator::tick_bfu_fetch() {
    for (uint f = 0; f < bfu_fetch_width; f++) {
        // Determine what reads need to be made to the register file.
        // As we don't ever write to the register file outside of PREFETCH,
        // we can support many parallel reads.

        BFUPipelineLatch& next = bfu_pipeline_latches[f][0];
        auto& pq = hardware_deques[f];
        if (pq.empty()) {
            // If priority queue is empty, invalidate the latch.
            next.valid = false;
            continue;
        }
        bfu_idle = false;
        if (next.stalled) {
            continue;  // We'll just repeat whatever is being done.
        }

        bool invalidate_next = false;

        std::deque<std::pair<uint, uint32_t>> proposed_matches;
        DequeEntry ei = pq.top();
        uint di = detector_vector_register[ei.next_unmatched_index];

        std::vector<std::pair<uint, uint32_t>> match_list;
        uint32_t min_weight = std::numeric_limits<uint32_t>::max();
        for (uint aj = ei.next_unmatched_index+1; 
                aj < detector_vector_register.size(); aj++) 
        {
            uint dj = detector_vector_register[aj];
            if (ei.running_matching.count(dj)) {
                continue;
            }
            addr_t address = std::make_pair(di, dj);
            bool is_hit = access(address, true);
            // If there is no hit on the register file, then we have to go
            // to DRAM. We will stall the pipeline until DRAM gives us
            // the value.
            auto vdi = graph.get_vertex(di);
            auto vdj = graph.get_vertex(dj);
            auto vdi_vdj = std::make_pair(vdi, vdj);
            if (is_hit) {
                // Get weight and add data to proposed_matches.
                fp_t raw_weight = path_table[vdi_vdj].distance;
                // Convert to integer (should be like this in DRAM).
                uint32_t w = (uint32_t)(raw_weight * MWPM_INTEGER_SCALE);
                match_list.push_back(std::make_pair(dj, w));
                if (w < min_weight) {
                    min_weight = w;
                }
            } else {
                invalidate_next = true;
            } 
        }
        
        if (!invalidate_next) {
            for (auto p : match_list) {
                uint dj = p.first;
                uint32_t w = p.second;
                if (w >= best_matching_register.matching_weight) {
                    continue;
                }
                if (detector_vector_register.size() <= HW_CUTOFF
                    || w <= weight_filter_cutoff
                    || w == min_weight)
                {
                    proposed_matches.push_back(p);
                }
            }
        }
        // Update the latch data.
        next.valid = !invalidate_next;
        next.proposed_matches = proposed_matches;
        next.base_entry = ei;
        if (!invalidate_next) {
            pq.pop();
        }
    }
}

bool
AstreaSimulator::access(addr_t address, bool set_evictable_on_hit) {
    // First, check that the data isn't already in the register file.
    // or in the replacement queue.
    bool is_hit = false;
    for (Register& r : register_file) {
        if (r.valid && r.address == address) {
            is_hit = true;
            r.last_use = 0;
            break;
        }
    }
    for (addr_t x : replacement_queue) {
        if (x == address) {
            is_hit = true;
            break;
        }
    }

    if (!is_hit) {
        replacement_queue.push_back(address);
    }

    return is_hit;
}

void
AstreaSimulator::update_state() {
    if (state == State::prefetch) {
        mean_weight_register /= access_counter;
#ifdef ASTREA_DEBUG
        std::cout << "[Mean Weight] " << mean_weight_register << "\n";
        std::cout << "BFU\n";
#endif
        state = State::bfu;
        hardware_deques[0].push((DequeEntry) {std::map<uint, uint>(), 0, 0});
        // Construct greedy initial matching.
        // Note: typically tracked as accesses are performed via DMA,
        // but we'll just do it here because it is easier.
        std::map<uint, uint> init_matching;
        fp_t matching_weight = 0.0;
        // Load the output register with a greedily computed matching
        // found during weight retrieval.
        for (uint i = 0; i < detector_vector_register.size(); i++) {
            uint di = detector_vector_register[i];
            auto vdi = graph.get_vertex(di);
            if (init_matching.count(di)) {
                continue;
            }
            DecodingGraph::Vertex * min_mate;
            fp_t min_weight = std::numeric_limits<fp_t>::max();
            for (uint j = i + 1; j < detector_vector_register.size(); j++) {
                uint dj = detector_vector_register[j];
                if (init_matching.count(dj)) {
                    continue;
                }
                auto vdj = graph.get_vertex(dj);
                fp_t w = path_table[std::make_pair(vdi, vdj)].distance;
                if (w < min_weight) {
                    min_mate = vdj;
                    min_weight = w;
                }
            }
            init_matching[di] = min_mate->detector;
            init_matching[min_mate->detector] = di;
            matching_weight += 
                path_table[std::make_pair(vdi, min_mate)].distance;
        }

        best_matching_register = (DequeEntry) {
            init_matching,
            (uint32_t) (matching_weight * MWPM_INTEGER_SCALE),
            0
        };
    } else if (state == State::bfu) {
#ifdef ASTREA_DEBUG
        std::cout << "IDLE\n";
        std::cout << "\tcycles in prefetch: " << prefetch_cycles << "\n";
        std::cout << "\tcycles in bfu: " << bfu_cycles << "\n";
        std::cout << "\thamming weight: " << detector_vector_register.size() << "\n";
#endif
        state = State::idle;
    }
}

void
AstreaSimulator::clear() {
    detector_vector_register.clear();
    mean_weight_register = 0;
    access_counter = 0;
    
    major_detector_register = 0;
    minor_detector_table.clear();

    for (auto& pq : hardware_deques) {
        pq.clear();
    }

    for (uint f = 0; f < bfu_fetch_width; f++) {
        for (auto& latch : bfu_pipeline_latches[f]) {
            latch.stalled = false;
            latch.valid = false;
        }
    }
    best_matching_register = {
        std::map<uint, uint>(), 
        std::numeric_limits<uint32_t>::max(),
        0 
    };

    mld_bm_register.fill({std::map<uint, uint>(), std::numeric_limits<uint32_t>::max(), 0});
    mld_prob_register.fill(0.0);

    replacement_queue.clear();
    state = State::prefetch;
#ifdef ASTREA_DEBUG
    std::cout << "(clear)PREFETCH\n";
#endif
}

//
// ASTREA SUPPLEMENTAL STRUCTURES
//

AstreaDeque::AstreaDeque(uint max_size)
    :max_size(max_size)
{}

void
AstreaDeque::push(AstreaSimulator::DequeEntry entry) {
    if (backing_array.size() == 0) {
        backing_array.push_back(entry);
    } else {
        backing_array.push_back(entry);
        upheap(backing_array.size()-1);
        if (backing_array.size() > max_size) {
            backing_array.pop_back();
        }
    }
}

AstreaSimulator::DequeEntry
AstreaDeque::top() {
    return backing_array.front();
}

void
AstreaDeque::pop() {
    auto back_entry = backing_array.back();
    backing_array[0] = back_entry;
    backing_array.pop_back();
    downheap(0);
}

uint
AstreaDeque::size() {
    return backing_array.size();
}

bool
AstreaDeque::empty() {
    return size() == 0;
}

void
AstreaDeque::clear() {
    backing_array.clear();
}

void
AstreaDeque::downheap(uint index) {
    if (index >= size()) {
        return;
    }
    auto curr = backing_array[index];
    uint left_index = left(index);
    uint right_index = right(index);
    if (left_index >= size() && right_index >= size()) {
        return;
    }
    if (left_index < size() && right_index < size()) {
        auto left_child = backing_array[left_index];
        auto right_child = backing_array[right_index];

        if (curr < left_child && curr < right_child) {
            return;
        }
        if (left_child < right_child) {
            // Swap parent with left child.
            backing_array[index] = left_child;
            backing_array[left_index] = curr;
            downheap(left_index);
        } else {
            // Swap with right child.
            backing_array[index] = right_child;
            backing_array[right_index] = curr;
            downheap(right_index);
        }
    } else if (left_index < size()) {
        auto left_child = backing_array[left_index];
        if (left_child < curr) {
            // Swap parent with left child.
            backing_array[index] = left_child;
            backing_array[left_index] = curr;
            downheap(left_index);
        } 
    } else {
        auto right_child = backing_array[right_index];
        if (right_child < curr) {
            // Swap with right child.
            backing_array[index] = right_child;
            backing_array[right_index] = curr;
            downheap(right_index);
        }
    }
}

void
AstreaDeque::upheap(uint index) {
    if (index <= 0) {
        return;
    }
    auto curr = backing_array[index];
    uint parent_index = parent(index);
    auto parent = backing_array[parent_index];
    if (curr < parent) {
        // Swap parent and child.
        backing_array[parent_index] = curr;
        backing_array[index] = parent;
        upheap(parent_index);
    }
}

uint
AstreaDeque::parent(uint index) {
    return (index-1) >> 1;
}

uint
AstreaDeque::left(uint index) {
    return ((index+1) << 1) + 1;
}

uint
AstreaDeque::right(uint index) {
    return (index+1) << 1;
}

uint
bound_detector(uint d, uint n_detectors) {
    return d == BOUNDARY_INDEX ? n_detectors-1 : d;
}

uint
unbound_detector(uint d, uint n_detectors) {
    return d == n_detectors-1 ? BOUNDARY_INDEX : d;
}

}  // astrea
}  // qrc
