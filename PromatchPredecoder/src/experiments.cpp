/*
 *  author: Narges Alavisamani
 *  date:   6 Feburary 2023
 * */

#include "experiments.h"

void error_chain_distr(std::ofstream& out, 
    uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er){
        
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }

    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );

    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);
    auto path_table = qrc::compute_path_table(decoding_graph);

    std::array<uint64_t, arr_N> ec_all_hw;
    std::array<uint64_t, arr_N> ec_only_hhw;

    std::array<uint64_t, arr_N> ec_all_hw_acc;
    std::array<uint64_t, arr_N> ec_only_hhw_acc;

    std::array<uint64_t, arr_N> all_hhw;
    std::array<uint64_t, arr_N> all_hhw_acc;

    std::array<uint64_t, arr_N> hw_for_ec_one;
    std::array<uint64_t, arr_N> hw_for_ec_one_acc;

    ec_all_hw.fill(0);
    ec_only_hhw.fill(0);

    ec_all_hw_acc.fill(0);
    ec_only_hhw_acc.fill(0);

    all_hhw.fill(0);
    all_hhw_acc.fill(0);

    hw_for_ec_one.fill(0);
    hw_for_ec_one_acc.fill(0);


    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            auto res = decoder.decode_error(syndrome);
            for (auto pair : res.matching) {
                auto v = decoding_graph.get_vertex(pair.first);
                auto w = decoding_graph.get_vertex(pair.second);
                auto v_w = std::make_pair(v, w);
                // Finding error chains from the path table of decoding graph 
                uint ec = path_table[v_w].path.size() - 1;
                // count the number of error chains. Each cell of array i in ec_all_hw
                // contains the number of error chains with length i. 
                ec_all_hw[ec]++;
                if (hw > 10) {
                    ec_only_hhw[ec]++;
                    all_hhw[hw] ++;
                    if (ec == 1){
                        hw_for_ec_one[hw] ++;
                    }
                }

            }
        }
        shots -= shots_this_round;
    }
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&ec_all_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        ec_all_hw_acc[i] = sum;

        MPI_Reduce(&ec_only_hhw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        ec_only_hhw_acc[i] = sum;

        MPI_Reduce(&all_hhw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        all_hhw_acc[i] = sum;

        MPI_Reduce(&hw_for_ec_one[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        hw_for_ec_one_acc[i] = sum;
    }

    if (world_rank == 0) {
        out << "Error chain length, number of ec with this length.\n";
        for (uint i = 0; i < arr_N; i++) {
            if (ec_all_hw_acc[i]) {
                out << i << "," << ec_all_hw_acc[i] << "\n";
            }
        }
        out << "Error chain length, number of ec with this length (high hamming weights).\n";
        for(uint i = 0; i < arr_N; i++){
            if (ec_only_hhw_acc[i]){
                out << i << "," << ec_only_hhw_acc[i] << "\n";
            }
        }

        out << "High hw value, number.\n";
        for(uint i = 0; i < arr_N; i++){
            if (all_hhw_acc[i]){
                out << i << "," << all_hhw_acc[i] << "\n";
            }
        }

        out << "(For error chain length 1)High hw value, number.\n";
        for(uint i = 0; i < arr_N; i++){
            if (hw_for_ec_one[i]){
                out << i << "," << hw_for_ec_one_acc[i] << "\n";
            }
        }
    }
}

void error_chain_distr_v2(std::ofstream& out, 
    uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er, FILE *file){
        
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }
    const stim::Circuit circ = (file == nullptr) ? qrc::build_circuit(
            distance,
            0,
            0,
            true,
            true,
            false,
            0,
            distance,
            physcial_error,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            meas_er //measurement error rates
        ):stim::Circuit::from_file(file);

    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);
    auto path_table = qrc::compute_path_table(decoding_graph);

    std::array<std::array<uint64_t, arr_N>,arr_N> ec_hw;
    std::array<std::array<uint64_t, arr_N>,arr_N> ec_hw_acc;

    ec_hw.fill(std::array<uint64_t, arr_N>{0});
    ec_hw_acc.fill(std::array<uint64_t, arr_N>{0});

    std::array<uint64_t, arr_N> ec_length;
    std::array<uint64_t, arr_N> ec_length_acc;

    ec_length.fill(0);
    ec_length_acc.fill(0);

    std::array<uint64_t, arr_N> hw_;
    std::array<uint64_t, arr_N> hw_acc;

    hw_.fill(0);
    hw_acc.fill(0);

    uint64_t ff = 0;
    uint64_t ff_c = 0;
    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            auto res = decoder.decode_error(syndrome);
            for (auto pair : res.matching) {
                auto v = decoding_graph.get_vertex(pair.first);
                auto w = decoding_graph.get_vertex(pair.second);
                auto v_w = std::make_pair(v, w);
                // Finding error chains from the path table of decoding graph 
                uint ec = path_table[v_w].path.size() - 1;
                // count the number of error chains. Each cell of array i in ec_all_hw
                // contains the number of error chains with length i. 

                ec_hw[hw][ec] ++;
                ec_length[ec]++;
                ff ++;
                hw_[hw]++;

            }
        }
        shots -= shots_this_round;
    }
    for (uint j = 0; j < arr_N; j++){
        for ( uint i = 0; i< arr_N;  i++){
            uint64_t sum = 0;
            MPI_Reduce(&ec_hw[j][i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
            ec_hw_acc[j][i] = sum;
        }
    }
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&ec_length[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        ec_length_acc[i] = sum;
    }
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&hw_[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        hw_acc[i] = sum;
    }
    uint64_t sum = 0;
    MPI_Reduce(&ff, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        bool first_line_printed = false;
        if(file != NULL)
            out << "Sycamore Simulation\n";
        else
            out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        for (uint j = 0; j<arr_N; j ++){
            first_line_printed = false;
            if(hw_acc[j]){
                for (uint i = 0; i < arr_N; i++) {
                    if (ec_hw_acc[j][i]) {
                        if(!first_line_printed){
                            first_line_printed = true;
                            out << "HW = " << j <<  "| Error chain length | # of ec with this length" << std::endl;
                            out << "Freq = "<< ((double)hw_acc[j]/sum) << std::endl;
                        }
                        out << i << " | " << ec_hw_acc[j][i]  <<std::endl;;
                    }
                }
            }
        }

        out << "**********_______________________________**********" << std::endl;
        out << "EC| | # of ec with this length | freq " << std::endl;
        for (uint j = 0; j<arr_N; j ++){
            if(ec_length_acc[j])
                out << j << " | " << ec_length_acc[j] << " | "<< ((double)ec_length_acc[j]/sum) << std::endl;
        }

        out << sum;
        
        
        
    }
}

void degree_one_distr(std::ofstream& out,
uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er){


       
    const uint64_t SHOTS_PER_BATCH = 100'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }

    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );

    qpd::Predecoder predecoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    std::array<uint64_t, arr_N> d1_hw;
    std::array<uint64_t, arr_N> d1_hw_acc;

    d1_hw.fill(0);
    d1_hw_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hw;
    number_of_hw.fill(0);

    std::array<uint64_t, arr_N> number_of_hw_acc;
    number_of_hw_acc.fill(0);

    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            predecoder.update_adjacency_matrix_and_detector_list(syndrome, 1, true);
            uint n_1s = predecoder.count_degree_one();
                // count the number of error chains. Each cell of array i in ec_all_hw
            // contains the number of error chains with length i. 
            d1_hw[hw] += n_1s;
            number_of_hw[hw] ++;

            

        }
        shots -= shots_this_round;
    }
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&d1_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        d1_hw_acc[i] = sum;

        MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hw_acc[i] = sum;
    }

    if (world_rank == 0) {
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "HW | #degree One | #syndromes with that HW | Average #degree 1\n";
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hw_acc[i]) {
                out << i << " | " << (d1_hw_acc[i])<< " | " << number_of_hw_acc[i] << " | "<<
                ((double)d1_hw_acc[i])/number_of_hw_acc[i] << "\n";
            }
        }
    }
}

std::string double_to_scientific_notation_string(double number){
    double exponent = 0;
    if (number != 0.0) {
        exponent = floor(std::log10(number));
    }
    double mantissa = number / std::pow(10, exponent);

    std::string exponent_string  = std::to_string(exponent);
    std::string mantissa_string = std::to_string(mantissa);

    size_t dot_pos = mantissa_string.find(".");
    if (dot_pos != std::string::npos)
        mantissa_string.erase(dot_pos+2);

    dot_pos = exponent_string.find(".");
    if (dot_pos != std::string::npos)
        exponent_string.erase(dot_pos);
    
    std::string result = mantissa_string + "e" + exponent_string;

    return result;
}

void HW_remained_distr(std::ofstream& out,uint64_t total_shots, 
                uint distance, fp_t physcial_error, fp_t meas_er, uint ec_length){
       
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }

    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );

    qpd::Predecoder predecoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    std::array<uint64_t, arr_N> acc_hw_after_prematch;
    std::array<uint64_t, arr_N> acc_hw_after_prematch_acc;

    acc_hw_after_prematch.fill(0);
    acc_hw_after_prematch_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hw;
    number_of_hw.fill(0);

    std::array<uint64_t, arr_N> number_of_hw_acc;
    number_of_hw_acc.fill(0);

    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            qpd::PredecoderShotResult prematching_results, ec_1_prematching_results;
            if(1 < ec_length ){
                ec_1_prematching_results = predecoder.matching_ec_length_one(syndrome);
                prematching_results = predecoder.matching_ec(ec_length, ec_1_prematching_results.post_syndrome);
            }
            else
                prematching_results = predecoder.matching_ec_length_one(syndrome);
                // count the number of error chains. Each cell of array i in ec_all_hw
            // contains the number of error chains with length i. 
            uint new_hw = std::accumulate(prematching_results.post_syndrome.begin(), prematching_results.post_syndrome.begin() + circ.count_detectors(), 0);
            acc_hw_after_prematch[hw] += new_hw;
            number_of_hw[hw] ++;

            

        }
        shots -= shots_this_round;
    }
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&acc_hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        acc_hw_after_prematch_acc[i] = sum;

        MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hw_acc[i] = sum;
    }

    if (world_rank == 0) {
        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hw_acc[i];
        std::cout << "t: " << t << std::endl;
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "HW | Accumulated HW After Prematching | #syndromes with that original HW | HW Frequency | Average HW after Prematching\n";
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hw_acc[i]) {
                out << i << " | " << (acc_hw_after_prematch_acc[i])<< " | " << number_of_hw_acc[i] << " | "<<
                ((double)number_of_hw_acc[i]/t) << " | "<< ((double)acc_hw_after_prematch_acc[i])/number_of_hw_acc[i] << "\n";
            }
        }
    }


}

void HW_distr(std::ofstream& out,uint64_t total_shots, 
                uint distance, fp_t physcial_error, fp_t meas_er){
       
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }

    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    //qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    std::array<uint64_t, arr_N> hw_n;
    hw_n.fill(0);

    std::array<uint64_t, arr_N> hw_n_acc;
    hw_n_acc.fill(0);

    // auto decoding_graph = qrc::to_decoding_graph(circ);
    // auto path_table = qrc::compute_path_table(decoding_graph);
    // bool found = false;

    // std::cout << "__________print_graph___________" << std::endl;
    // std::vector<qrc::DecodingGraph::Vertex*> _vertices = decoding_graph.vertices();
    // for(qrc::DecodingGraph::Vertex* v : _vertices){
    //     for(qrc::DecodingGraph::Vertex* w : _vertices){
    //         auto v_w = std::make_pair(v, w);
    //         // Finding error chains from the path table of decoding graph 
    //         uint ec = path_table[v_w].path.size() - 1;
    //         // std::vector<qrc::DecodingGraph::Vertex*> ad_list = decoding_graph.adjacency_list(v);
    //         if(!found && ec == 1){
    //             std::cout << "\n___Vertex__ " << v->detector << std::endl;
    //             std::cout << w->detector << ", ";
    //             found = true;
    //         }
    //         if (ec == 1){  
    //             std::cout << w->detector << ", ";
    //         }
    //     }
    //     found = false;
    // }

    // return;


    while (shots > 0) {
        if(shots % 2'000'000 == 0 && world_rank == 0)
            std::cout <<  world_rank << ": "<< shots << std::endl;
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            hw_n[hw] ++;

        }
        shots -= shots_this_round;
    }
    
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&hw_n[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        hw_n_acc[i] = sum;
    }

    if (world_rank == 0) {
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "HW | #Of HW" << std::endl;
        for (uint i = 0; i < arr_N; i++) {
            if (hw_n_acc[i]) {
                out << i << " | " << hw_n_acc[i] << std::endl;
            }
        }
    }
}

// This function calculates the frequnecy of HWs after prematching regardless of their original HW
void HW_remained_distr_v2(std::ofstream& out,uint64_t total_shots, 
                uint distance, fp_t physcial_error, fp_t meas_er, uint ec_length){
       
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }

    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );

    qpd::Predecoder predecoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    std::array<uint64_t, arr_N> number_hw_after_prematch;
    std::array<uint64_t, arr_N> number_hw_after_prematch_acc;

    number_hw_after_prematch.fill(0);
    number_hw_after_prematch_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hw;
    number_of_hw.fill(0);

    std::array<uint64_t, arr_N> number_of_hw_acc;
    number_of_hw_acc.fill(0);

    while (shots > 0) {
        if(shots % 1000 == 0)
            std::cout << shots << std::endl;
        
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            
            qpd::PredecoderShotResult prematching_results, ec_1_prematching_results;
            if(1 < ec_length ){
                ec_1_prematching_results = predecoder.matching_ec_length_one(syndrome);

                prematching_results = predecoder.matching_ec(ec_length, ec_1_prematching_results.post_syndrome);
            }
            else
                prematching_results = predecoder.matching_ec_length_one(syndrome);
                // count the number of error chains. Each cell of array i in ec_all_hw
            // contains the number of error chains with length i. 
            uint new_hw = std::accumulate(prematching_results.post_syndrome.begin(), prematching_results.post_syndrome.begin() + circ.count_detectors(), 0);
            number_hw_after_prematch[new_hw]++ ;
            number_of_hw[hw] ++;

            

        }
        shots -= shots_this_round;
    }
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&number_hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_hw_after_prematch_acc[i] = sum;

        MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hw_acc[i] = sum;
    }

    if (world_rank == 0) {
        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hw_acc[i];
        std::cout << "t: " << t << std::endl;
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "HW | #syndromes with that original HW | HW Frequency | #syndromes with that HW after prematching | After prematching HW Frequency \n";
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hw_acc[i]) {
                out << i << " | " << number_of_hw_acc[i] << " | "<<
                ((double)number_of_hw_acc[i]/t) << " | "<< number_hw_after_prematch_acc[i] << " | "<<
                ((double)number_hw_after_prematch_acc[i]/t) << std::endl;
            }
        }
    }


}

std::string number_of_shots_to_string(uint64_t shots_number){
    std::string result;
    std::string postfix;
    uint32_t  divisor = shots_number/ 1'000;
    postfix = "K";

    if(divisor>999){
        divisor = shots_number/ 1'000'000;
        postfix = "M";
    }

    if(divisor>999){
        divisor = shots_number/ 1'000'000'000;
        postfix = "B";
    }
    if(divisor>999){
        divisor = shots_number/ 1'000'000'000'000;
        postfix = "T";
    }
    result = std::to_string(divisor) + postfix;


    return result;
}

std::string generate_experiment_name(std::string exp_type, uint distance, 
fp_t physical_er, fp_t m_error, uint64_t total_shots, bool mention_shots){
    
    std::string result = exp_type + 
    "_p=" + double_to_scientific_notation_string(physical_er)+
    "_m=" + double_to_scientific_notation_string(m_error)+
    "_d=" + std::to_string(distance);
    if(mention_shots){
        result += "_s=" + number_of_shots_to_string(total_shots);
    }
    result +=".txt";

    return result;
}

std::string generate_matching_file_name(std::string exp_type, uint distance, 
fp_t physical_er, fp_t m_error, uint64_t total_shots, bool mention_shots){
    
    std::string result = "matchings_"+ exp_type + 
    "_p=" + double_to_scientific_notation_string(physical_er)+
    "_m=" + double_to_scientific_notation_string(m_error)+
    "_d=" + std::to_string(distance);
    
    if(exp_type == "exp4" or exp_type == "exp5" or mention_shots){
        result += "_s=" + number_of_shots_to_string(total_shots);
    }
    result +=".txt";

    return result;
}

void print_info(const std::vector<uint8_t>& syndrome,
                qpd::PredecoderShotResult predecoder_result, uint n_detectors
                , bool print_syndrome, bool check_syndrome_correctness){
    
    if (print_syndrome){
        std::cout << "Syndrome:" << std::endl;
        for (int id =0; id<syndrome.size();id++) {
            std::cout<< id<<":"<< ((int)syndrome[id] )<< ", ";
        }
        std::cout << std::endl;
        std::cout << "Predecoder Syndrome:" << std::endl;
        for (int id =0; id<predecoder_result.prematch_syndrome.size();id++){
            std::cout<< id<<":" << ((int)predecoder_result.prematch_syndrome[id] ) << ", ";
        }
        std::cout << std::endl;
        std::cout << "Post Syndrome:" << std::endl;
        for (int id =0; id<predecoder_result.post_syndrome.size();id++){
            std::cout<< id<<":" << ((int)predecoder_result.post_syndrome[id] ) << ", ";
        }
    }
    std::cout << std::endl;

    std::cout << "_____________________Prematching MAP_____________________" << std::endl;
    for (const auto& kv : predecoder_result.pre_matches) {
        std::cout << "V1: " << kv.first << ", V2: " << kv.second << std::endl;
    }
    std::cout << "_____________________End of MAP_____________________" << std::endl;

    if(check_syndrome_correctness){
        for (uint i = 0; i < n_detectors; i++){
            if( predecoder_result.post_syndrome[i] != syndrome[i]){
                std::cout << "Diff between current and post syndrome. Detector id: " << i  << " prematch bit: " <<
                ((int)predecoder_result.prematch_syndrome[i]) << std::endl;
            }
            if(predecoder_result.prematch_syndrome[i]!=((uint)0))
                std::cout<< "None zero bit found at "<< i << endl;

        }
    }

}



void check_subgraphs(uint distance, 
fp_t physcial_error, fp_t meas_er, uint ec_length, bool print_l1_prematches){


    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );

    qrc::Decoder* decoder  = new qrc::MWPMDecoder(circ);

    
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);
    predecoder->set_decoder(decoder);
    
    std::cout << "Enter the Syndrome: " << std::endl;
    std::string numbers;
    std::getline(std::cin, numbers);

    std::vector<uint16_t> result;
    std::istringstream iss(numbers);
    uint16_t num;

    while (iss >> num) {
        result.push_back(num);
    }

    const uint n_detectors = circ.count_detectors();
    std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(result, distance);
    std::cout << std::endl << "HW = " <<std::accumulate(syndrome.begin(), syndrome.begin()+predecoder->n_detectors,0) << std::endl;

    predecoder->print_subgraph(syndrome);

}

void predecoder_plus_MWPM_decoder(std::ofstream& out,uint64_t total_shots, 
uint distance, fp_t physcial_error, fp_t meas_er){
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }
    const stim::Circuit circ = qrc::build_circuit(
            distance,
            0,
            0,
            true,
            true,
            false,
            0,
            distance,
            physcial_error,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            meas_er //measurement error rates
        );

    std::array<uint64_t, arr_N> acc_hw_after_prematch;
    std::array<uint64_t, arr_N> acc_hw_after_prematch_acc;

    acc_hw_after_prematch.fill(0);
    acc_hw_after_prematch_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hw;
    number_of_hw.fill(0);

    std::array<uint64_t, arr_N> number_of_hw_acc;
    number_of_hw_acc.fill(0);

    qpd::Predecoder predecoder(circ);
    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    const uint n_observables = circ.count_observables();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    uint64_t n_predecoder_decoder_logical_error = 0;
    uint64_t n_decoder_logcial_error = 0;
    uint64_t n_better_with_predec = 0;
    uint64_t n_better_with_predec_acc = 0;
    uint64_t n_predecoder_decoder_logical_error_acc = 0;
    uint64_t n_decoder_logcial_error_acc = 0;
    bool predecoder_decoder_error = false;
    bool decoder_error = false;

    //std :: cout << "HERE"<< std::endl;
    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
            //predecoder.print_subgraph(syndrome);
            
            //Predecoder + MWPM decoder
            auto prematch_results = predecoder.matching_ec_length_one(syndrome);
            //std::cout << "1, "  << prematch_results.pre_matches.size() << std::endl;
            auto predecoder_results = predecoder.predecode_error(prematch_results);
            //std::cout << "2" << std::endl;
            auto predecoder_decoder_results = decoder.decode_error(predecoder_results.post_syndrome);
            uint new_hw = std::accumulate(predecoder_results.post_syndrome.begin(), predecoder_results.post_syndrome.begin() + circ.count_detectors(), 0);
            acc_hw_after_prematch[hw] += new_hw;
            number_of_hw[hw] ++;
            //std::cout << "3" << std::endl;
            std::vector<uint8_t> final_correction = predecoder_decoder_results.correction;
            // std::cout << "4" << " n_observables" << n_observables << " " << circ.count_observables() << " n_detectors" << n_detectors<< 
            // "  "  << circ.count_detectors() << "-- "<< syndrome.size()<< " CSIZE"<<final_correction.size() << std::endl;
            for (uint i = 0; i < n_observables; i++) {
                final_correction[i] ^= predecoder_results.correction[i];
            }
            
            //std::cout << "5" << std::endl;
            predecoder_decoder_error = qrc::is_logical_error(final_correction,syndrome,n_detectors, n_observables);
            if(predecoder_decoder_error){
                n_predecoder_decoder_logical_error++;
                //std::cout << "ERROR"<< std::endl;
                // qpd::print_syndrome(syndrome);
                //std::cout << "size of all "<< final_correction.size()<< ", xx: "<< ((int) final_correction[0])<<std::endl;
                // std::cout << "Size of pred: "<< ", CVVVV: "<< predecoder_results.correction.size()<< ", xx: "<< ((int) predecoder_results.correction[0])<<std::endl;
            }
            else{
                // std::cout << " -- SUCCESS"<< std::endl;
                //qpd::print_syndrome(syndrome);
                // std::cout << "size of all "<< final_correction.size()<< ", xx: "<< ((int) final_correction[0])<<std::endl;
                // std::cout << "Size of pred: "<< predecoder_results.correction.size()<< ", xx: "<< ((int) predecoder_results.correction[0])<<std::endl;
            

            }
            
            
            //MWPM Decoder
            //std::cout << "6" << std::endl;
            auto decoder_results = decoder.decode_error(syndrome);
            //std::cout << "7" << std::endl;
            decoder_error = qrc::is_logical_error(decoder_results.correction,syndrome,n_detectors, n_observables);
            //std::cout << "8" << std::endl;
            if(decoder_error)
                n_decoder_logcial_error++;
            //std::cout << "Just before condition " << std::endl;
            if(decoder_error==false && predecoder_decoder_error==true){
                uint number_of_degree_1 = predecoder.count_degree_one();
                // predecoder.print_shot_info(syndrome, prematch_results, 
                //                 predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
                //std::cout << "end of the if " << std::endl;
            }
            //std::cout << "Before starting next round " << std::endl;
            if(decoder_error==true && predecoder_decoder_error==false){
                n_better_with_predec++;
                //std::cout << "**********BETTER MATCHING WITH PREDECODER!*********"<<std::endl;
                uint number_of_degree_1 = predecoder.count_degree_one();
                // predecoder.print_shot_info(syndrome, prematch_results, 
                //                predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
                //std::cout << "***************************************************"<<std::endl;

            }

            
        }
        shots -= shots_this_round;
    }

    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&acc_hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        acc_hw_after_prematch_acc[i] = sum;

        MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hw_acc[i] = sum;
    }

    uint64_t sum = 0;
    MPI_Reduce(&n_predecoder_decoder_logical_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_predecoder_decoder_logical_error_acc = sum;

    MPI_Reduce(&n_decoder_logcial_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_decoder_logcial_error_acc = sum;

    MPI_Reduce(&n_better_with_predec, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_better_with_predec_acc = sum;
    //std::cout << "Before gathter" << std:: endl;
    // MPI_Gather(output_str.c_str(), output_str.size() + 1, MPI_CHAR, 
    //        output_strings.data(), output_str.size() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    //std::cout << "after" << std:: endl;

    if (world_rank == 0) {
    
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        
        out << "Predecoder + MWPM Decoder # logical error: " << n_predecoder_decoder_logical_error_acc << std::endl;
        out << "MWPM Decoder # logical error: " << n_decoder_logcial_error_acc << std::endl;
        out << "# Times MWPM fails but Predecoder + MWPM deos not: " << n_better_with_predec_acc << std::endl;

        out << "***********************************************************" <<std::endl;

        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hw_acc[i];
        out << "HW | Accumulated HW After Prematching | #syndromes with that original HW | HW Frequency | Average HW after Prematching\n";
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hw_acc[i]) {
                out << i << " | " << (acc_hw_after_prematch_acc[i])<< " | " << number_of_hw_acc[i] << " | "<<
                ((double)number_of_hw_acc[i]/t) << " | "<< ((double)acc_hw_after_prematch_acc[i])/number_of_hw_acc[i] << "\n";
            }
        }
        
        // std::ostringstream final_output;
        // for (const auto& output : output_strings) {
        //     final_output << output;
        // }
        // out << final_output.str();

    }
}

void decoding_cost_experiments(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er){
        
       
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }

    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );

    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    std::array<uint64_t, arr_N> acc_time_hw;
    std::array<uint64_t, arr_N> acc_time_hw_acc;

    acc_time_hw.fill(0);
    acc_time_hw_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hw;
    number_of_hw.fill(0);

    std::array<uint64_t, arr_N> number_of_hw_acc;
    number_of_hw_acc.fill(0);

    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            auto res = decoder.decode_error(syndrome);

            acc_time_hw[hw] += res.execution_time;
            number_of_hw[hw] ++;

        }
        shots -= shots_this_round;
    }
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&acc_time_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        acc_time_hw_acc[i] = sum;

        MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hw_acc[i] = sum;
    }

    if (world_rank == 0) {
        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hw_acc[i];
        std::cout << "t: " << t << std::endl;
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "HW | Mean Decoding Execution Time | # Times for the HW" << std::endl;
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hw_acc[i]) {
                out << i << " | "<< ((double)acc_time_hw_acc[i])/number_of_hw_acc[i] << " | " << 
                number_of_hw_acc[i] << std::endl;
            }
        }
    }

}

// void measurements_stacked(uint distance, fp_t physcial_error, fp_t meas_er, uint ec_length = 1, bool print_l1_prematches = false){
//     const uint64_t SHOTS_PER_BATCH = 1'000'000;

//     int world_size, world_rank;
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

//     uint64_t shots = total_shots / world_size;
//     if (world_rank == world_size - 1) {
//         shots += total_shots % world_size;
//     }

//     const stim::Circuit circ = qrc::build_circuit(
//         distance,
//         0,
//         0,
//         true,
//         true,
//         false,
//         0,
//         distance,
//         physcial_error,
//         0,
//         0,
//         0,
//         -1,
//         -1,
//         -1,
//         meas_er //measurement error rates
//     );

//     qpd::Predecoder predecoder(circ);
//     std::mt19937_64 rng(world_rank);

//     const uint n_detectors = circ.count_detectors();
//     qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

//     std::array<uint64_t, arr_N> number_hw_after_prematch;
//     std::array<uint64_t, arr_N> number_hw_after_prematch_acc;

//     number_hw_after_prematch.fill(0);
//     number_hw_after_prematch_acc.fill(0);

//     std::array<uint64_t, arr_N> number_of_hw;
//     number_of_hw.fill(0);

//     std::array<uint64_t, arr_N> number_of_hw_acc;
//     number_of_hw_acc.fill(0);

//     while (shots > 0) {
//         uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
//         // Make syndromes
//         auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
//         sample_buffer = sample_buffer.transposed();

//         for (uint64_t i = 0; i < shots_this_round; i++) {
//             auto syndrome = qrc::_to_vector(sample_buffer[i], 
//                                     circ.count_detectors(),
//                                     circ.count_observables());
//             // Calculate the hamming weights
//             uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            
//             // qpd::PredecoderShotResult prematching_results, ec_1_prematching_results;
//             // if(1 < ec_length ){
//             //     ec_1_prematching_results = predecoder.matching_ec_length_one(syndrome);

//             //     prematching_results = predecoder.matching_ec(ec_length, ec_1_prematching_results.post_syndrome);
//             // }
//             // else
//             //     prematching_results = predecoder.matching_ec_length_one(syndrome);
//             //     // count the number of error chains. Each cell of array i in ec_all_hw
//             // // contains the number of error chains with length i. 
//             // uint new_hw = std::accumulate(prematching_results.post_syndrome.begin(), prematching_results.post_syndrome.begin() + circ.count_detectors(), 0);
//             // number_hw_after_prematch[new_hw]++ ;
//             // number_of_hw[hw] ++;

//         }
//         shots -= shots_this_round;
//     }
//     for ( uint i = 0; i< arr_N;  i++){
//         uint64_t sum = 0;
//         MPI_Reduce(&number_hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//         number_hw_after_prematch_acc[i] = sum;

//         MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//         number_of_hw_acc[i] = sum;
//     }

//     if (world_rank == 0) {
//         uint64_t t = 0;
//         for (uint i=0; i<arr_N; i++)
//             t+=number_of_hw_acc[i];
//         std::cout << "t: " << t << std::endl;
//         out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
//         out << "HW | #syndromes with that original HW | HW Frequency | #syndromes with that HW after prematching | After prematching HW Frequency \n";
//         for (uint i = 0; i < arr_N; i++) {
//             if (number_of_hw_acc[i]) {
//                 out << i << " | " << number_of_hw_acc[i] << " | "<<
//                 ((double)number_of_hw_acc[i]/t) << " | "<< number_hw_after_prematch_acc[i] << " | "<<
//                 ((double)number_hw_after_prematch_acc[i]/t) << std::endl;
//             }
//         }
//     }
// }

void print_decoding_results(const std::vector<uint8_t>& syndrome, 
                            qrc::DecoderShotResult results, uint n_detectors){
    
    std::map<uint8_t,bool> flipped_bits;

    std::cout << std::endl << "_";
    for (int id =0; id< n_detectors; id++) {
        if((int)syndrome[id] )
            std::cout<<id<< ", ";
            flipped_bits[id] = false;
    }
    std::cout << std::endl << results.is_logical_error << std::endl;
    for (int id =0; id<results.correction.size();id++) {
            std::cout<< ((int)results.correction[id] )<< ",";
    }

    std::cout << std::endl;
    for(auto matches : results.matching){
        if(!flipped_bits[matches.first])
            std::cout<< matches.first << "," << matches.second << "|";
            flipped_bits[matches.first] = true;
            flipped_bits[matches.second] = true;

    }
}

void high_distance_experiments(std::ofstream& out,uint64_t total_shots, uint distance, fp_t physcial_error, fp_t meas_er){
        
       
    const uint64_t SHOTS_PER_BATCH = 10;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }

    const stim::Circuit circ = qrc::build_circuit(
        distance,
        0,
        0,
        true,
        true,
        false,
        0,
        distance,
        physcial_error,
        0,
        0,
        0,
        -1,
        -1,
        -1,
        meas_er //measurement error rates
    );

    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    std::array<uint64_t, arr_N> acc_time_hw;
    std::array<uint64_t, arr_N> acc_time_hw_acc;

    acc_time_hw.fill(0);
    acc_time_hw_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hw;
    number_of_hw.fill(0);

    std::array<uint64_t, arr_N> number_of_hw_acc;
    number_of_hw_acc.fill(0);

    printf("Hello world , rank %d out of %d processors: %d\n"
           , world_rank, world_size, shots);

    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());

            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            auto res = decoder.decode_error(syndrome);

            acc_time_hw[hw] += res.execution_time;
            number_of_hw[hw] ++;
            if (world_rank == 0) {
                std::cout << world_rank << std::endl;
                print_decoding_results(syndrome, res, n_detectors);
            }
        }
        shots -= shots_this_round;
    }
    
    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&acc_time_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        acc_time_hw_acc[i] = sum;

        MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hw_acc[i] = sum;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hw_acc[i];
        std::cout << "t: " << t << std::endl;
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "HW | Mean Decoding Execution Time | # Times for the HW" << std::endl;
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hw_acc[i]) {
                out << i << " | "<< ((double)acc_time_hw_acc[i])/number_of_hw_acc[i] << " | " << 
                number_of_hw_acc[i] << std::endl;
            }
        }
    }

}

void print_decoding_results_v2(const std::vector<uint8_t>& syndrome, qrc::DecoderShotResult results,
 uint n_detectors, std::map<uint, uint> pre_matches, qrc::DecoderShotResult ppost_results, std::ofstream& out){
    std::map<uint8_t,bool> flipped_bits;
    std::map<uint8_t,bool> flipped_bits2;

    out << std::endl << "_";
    for (int id =0; id< n_detectors; id++) {
        if((int)syndrome[id] )
            out<<id<< ", ";
            flipped_bits[id] = false;
            flipped_bits2[id] = false;
    }
    out << std::endl << results.is_logical_error << std::endl;
    for (int id =0; id<results.correction.size();id++) {
            out<< ((int)results.correction[id] )<< ",";
    }

    out << std::endl;
    for(auto matches : results.matching){
        if(!flipped_bits[matches.first])
            out<< matches.first << "," << matches.second << "|";
            flipped_bits[matches.first] = true;
            flipped_bits[matches.second] = true;

    }
    out << std::endl << "*" <<ppost_results.is_logical_error << std::endl;
    for(auto matches : pre_matches){
        if(!flipped_bits2[matches.first])
            out<< matches.first << "," << matches.second << ")";
            flipped_bits2[matches.first] = true;
            flipped_bits2[matches.second] = true;
    }
    for(auto matches : ppost_results.matching){
        if(!flipped_bits2[matches.first])
            out<< matches.first << "," << matches.second << "|";
            flipped_bits2[matches.first] = true;
            flipped_bits2[matches.second] = true;
    }
 }


void high_distance_predecoder_plus_MWPM_decoder(std::ofstream& out,uint64_t total_shots, 
uint distance, fp_t physcial_error, fp_t meas_er, std::ofstream& out1){
    const uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }
    const stim::Circuit circ = qrc::build_circuit(
            distance,
            0,
            0,
            true,
            true,
            false,
            0,
            distance,
            physcial_error,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            meas_er //measurement error rates
        );

    std::array<uint64_t, arr_N> acc_hw_after_prematch;
    std::array<uint64_t, arr_N> acc_hw_after_prematch_acc;

    acc_hw_after_prematch.fill(0);
    acc_hw_after_prematch_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hw;
    number_of_hw.fill(0);

    std::array<uint64_t, arr_N> number_of_hw_acc;
    number_of_hw_acc.fill(0);

    qpd::Predecoder predecoder(circ);
    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    const uint n_observables = circ.count_observables();
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    uint64_t n_predecoder_decoder_logical_error = 0;
    uint64_t n_decoder_logcial_error = 0;
    uint64_t n_better_with_predec = 0;
    uint64_t n_better_with_predec_acc = 0;
    uint64_t n_predecoder_decoder_logical_error_acc = 0;
    uint64_t n_decoder_logcial_error_acc = 0;
    bool predecoder_decoder_error = false;
    bool decoder_error = false;

    //std :: cout << "HERE"<< std::endl;
    while (shots > 0) {
        std::cout << "S: " << shots << std::endl;
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            if (i % 50'000 == 0)
                std::cout << i << std::endl;
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
            //predecoder.print_subgraph(syndrome);
            number_of_hw[hw] ++ ;
            //Predecoder + MWPM decoder (matching_ec_length_one finds the matching
            // and send the results to predecode to find the corrections)
            auto prematch_results = predecoder.matching_ec_length_one(syndrome);
            //std::cout << "1, "  << prematch_results.pre_matches.size() << std::endl;
            auto predecoder_results = predecoder.predecode_error(prematch_results);
            //std::cout << "2" << std::endl;
            auto predecoder_decoder_results = decoder.decode_error(predecoder_results.post_syndrome);
            uint new_hw = std::accumulate(predecoder_results.post_syndrome.begin(), predecoder_results.post_syndrome.begin() + circ.count_detectors(), 0);
            acc_hw_after_prematch[new_hw] ++;

            //std::cout << "3" << std::endl;
            std::vector<uint8_t> final_correction = predecoder_decoder_results.correction;
            // std::cout << "4" << " n_observables" << n_observables << " " << circ.count_observables() << " n_detectors" << n_detectors<< 
            // "  "  << circ.count_detectors() << "-- "<< syndrome.size()<< " CSIZE"<<final_correction.size() << std::endl;
            for (uint i = 0; i < n_observables; i++) {
                final_correction[i] ^= predecoder_results.correction[i];
            }
            //std::cout << "5" << std::endl;
            predecoder_decoder_error = qrc::is_logical_error(final_correction,syndrome,n_detectors, n_observables);
            if(predecoder_decoder_error)
                n_predecoder_decoder_logical_error++;
            
            //MWPM Decoder
            //std::cout << "6" << std::endl;
            auto decoder_results = decoder.decode_error(syndrome);
            print_decoding_results_v2(syndrome, decoder_results,n_detectors, 
            prematch_results.pre_matches, predecoder_decoder_results,out1);

            //std::cout << "7" << std::endl;
            //decoder_error = qrc::is_logical_error(decoder_results.correction,syndrome,n_detectors, n_observables);
            //std::cout << "8" << std::endl;
            if(decoder_results.is_logical_error )
                n_decoder_logcial_error++;            
        }
        shots -= shots_this_round;
    }

    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&acc_hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        acc_hw_after_prematch_acc[i] = sum;

        MPI_Reduce(&number_of_hw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hw_acc[i] = sum;
    }

    uint64_t sum = 0;
    MPI_Reduce(&n_predecoder_decoder_logical_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_predecoder_decoder_logical_error_acc = sum;

    MPI_Reduce(&n_decoder_logcial_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_decoder_logcial_error_acc = sum;
    //std::cout << "Before gathter" << std:: endl;
    // MPI_Gather(output_str.c_str(), output_str.size() + 1, MPI_CHAR, 
    //        output_strings.data(), output_str.size() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    //std::cout << "after" << std:: endl;

    if (world_rank == 0) {
    
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        
        out << "Predecoder + MWPM Decoder # logical error: " << n_predecoder_decoder_logical_error_acc << std::endl;
        out << "MWPM Decoder # logical error: " << n_decoder_logcial_error_acc << std::endl;
        out << "# Times MWPM fails but Predecoder + MWPM deos not: " << n_better_with_predec_acc << std::endl;

        out << "***********************************************************" <<std::endl;

        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hw_acc[i];
        out << "HW | #syndromes with that original HW | HW Frequency | #Post Pre | Post HW Frequency\n";
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hw_acc[i]) {
                out << i << " | " << number_of_hw_acc[i] << " | " << ((double)number_of_hw_acc[i]/t) << " | "<<
                (acc_hw_after_prematch_acc[i])<< " | " <<
                ((double)acc_hw_after_prematch_acc[i])/t << "\n";
            }
        }
        
        // std::ostringstream final_output;
        // for (const auto& output : output_strings) {
        //     final_output << output;
        // }
        // out << final_output.str();

    }
}

std::string syndrome_file_name(std::string experiment_name, uint distance, fp_t physical_er, 
                        fp_t m_error, uint64_t total_shots, uint8_t core_number, uint64_t round, bool partially){
        

    std::string result = experiment_name + 
    "_p=" + double_to_scientific_notation_string(physical_er)+
    "_m=" + double_to_scientific_notation_string(m_error)+
    "_d=" + std::to_string(distance)+
    "_s=" + number_of_shots_to_string(total_shots);
    
    if(partially)
        return result;
    result = result +
    "_c=" + std::to_string(core_number)+
    "_r=" + std::to_string(round)+
    ".txt";

    return result;  

}

void saving_hhw_syndromes(std::ofstream& out, uint64_t total_shots, 
uint distance, fp_t physcial_error, fp_t meas_er, std::string syndrome_folder_name, bool print_logs){
    const uint64_t SHOTS_PER_BATCH = 1'000'000;
    uint64_t hhw_n = 0;
    uint64_t hhw_n_acc = 0;

    uint64_t n_saved = 0;
    uint64_t n_saved_acc = 0;

    uint64_t n_ler = 0;
    uint64_t n_ler_acc = 0;
    

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); 

    if (world_rank == 0){
        if (opendir(syndrome_folder_name.c_str()) == NULL) {
            if (mkdir(syndrome_folder_name.c_str(), 0777) == 0) {
                std::cout << "Directory created successfully\n";
            } else {
                std::cout << "Failed to create directory\n";
            }
        } 
        else {
            std::cout << "Directory already exists\n";
        } 

    // To check if we are reading correct files
    }

    uint64_t shots = total_shots / world_size;
    if (world_rank == world_size - 1) {
        shots += total_shots % world_size;
    }
    uint64_t round = 1;
    const stim::Circuit circ = qrc::build_circuit(
            distance,
            0,
            0,
            true,
            true,
            false,
            0,
            distance,
            physcial_error,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            meas_er //measurement error rates
        );

    std::mt19937_64 rng(world_rank);
    const uint n_detectors = circ.count_detectors();
    std::vector<vector<vector<uint16_t>>> saved_hhw_syndromes;
    qrc::MWPMDecoder mwpm_decoder(circ);


    //std::vector<uint16_t>> round;
    
    saved_hhw_syndromes.resize(world_size);
    //round.resize();
     std::vector<uint16_t> flipped_bits;
    //std :: cout << "HERE"<< std::endl;
    while (shots > 0) {
        if(shots % 2'000'000 == 0 && world_rank == 0 && print_logs) //2500
            std::cout << shots << std::endl;
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            // Calculate the hamming weights
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
            if( hw < 11)
                continue;
            // print_graph(graph);
            hhw_n ++;
            if(mwpm_decoder.decode_error(syndrome).is_logical_error){
                n_ler ++;
            }
            flipped_bits = qpd::syndrome_compressed(syndrome);
            saved_hhw_syndromes[world_rank].push_back(flipped_bits);
            

        }
        shots -= shots_this_round;

        // if(10'000'000 < saved_hhw_syndromes[world_rank].size()){
        //     std::string file_name = directory_addr + syndrome_file_name("sgen",distance,physcial_error,
        //     meas_er, total_shots, world_rank, round);
        //     qpd::write_vector_of_vectors(saved_hhw_syndromes[world_rank], file_name);
        //     n_saved += (uint64_t)saved_hhw_syndromes[world_rank].size();
        //     saved_hhw_syndromes[world_rank].clear();
        //     round++;
        // }
    }


    if(int(saved_hhw_syndromes[world_rank].size()) != 0){
        std::string file_name = syndrome_folder_name + syndrome_file_name("sgen",distance,physcial_error,
        meas_er, total_shots, world_rank, round);
        qpd::write_vector_of_vectors(saved_hhw_syndromes[world_rank], file_name);
        n_saved += (uint64_t)saved_hhw_syndromes[world_rank].size();
        saved_hhw_syndromes[world_rank].clear();
        round++;
    }

    uint64_t sum = 0;
    MPI_Reduce(&hhw_n, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    hhw_n_acc = sum;
    sum = 0;
    MPI_Reduce(&n_saved, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_saved_acc = sum;
    sum = 0;
    MPI_Reduce(&n_ler, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_ler_acc = sum;

    if (world_rank == 0) {
        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << " n_detectors:" << n_detectors << std::endl;
        out << "# Of High HW: " << hhw_n_acc << ", " <<n_saved_acc << std::endl;
        out << "Astrea LER: " << ((double)hhw_n_acc/ total_shots) << std::endl << "MWPM LER: " << ((double)n_ler_acc/total_shots);
    }


    //std::cout << "Before gathter" << std:: endl;
    // MPI_Gather(output_str.c_str(), output_str.size() + 1, MPI_CHAR, 
    //        output_strings.data(), output_str.size() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    //std::cout << "after" << std:: endl;

}

void predecoding_hhw_syndromes(std::ofstream& out, std::string syndromes_directory, uint64_t total_shots, uint distance, 
fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry){

    //data_folder_name += ("syndromes/" + number_of_shots_to_string(total_shots));
    
    // for (const auto & entry : std::filesystem::directory_iterator(data_folder_name))
    //     synd_files.push_back(entry.path());

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    bool abort_f = false;     
    std::vector<std::vector<uint16_t>> buffer;

    if(world_rank == 0){
        std::cout << "Running predecoding_hhw_syndromes: Reading from " << syndromes_directory<<std::endl;
        // while ((dp = readdir(dirp)) != NULL) {
        //     if (dp->d_type == DT_REG) { // If the entry is a regular file
        //         std::string file_name = dp->d_name;
        //         std::string file_path = data_folder_name + "/" + file_name;
        //         synd_files.push_back(file_path);
        //     }
        // }
        // closedir(dirp);


        DIR* dirp = opendir(syndromes_directory.c_str());
        if (dirp == nullptr) {
            std::cout << "Error opening directory" << std::endl;
            abort_f = true;
        }
        std::cout << "These are two of the files in the directory " << syndromes_directory << ": ";
        int count = 0;
        struct dirent* dp;
        while ((dp = readdir(dirp)) != nullptr) {
            if (count == 2) { // Found two files, break out of the loop
                break;
            }

            if (dp->d_type == DT_REG && std::strstr(dp->d_name, ".txt") != nullptr) { // Check if the file is a regular file with a .txt extension
                std::cout << dp->d_name << std::endl;
                count++;
            }
        }
        closedir(dirp);

        std::cout <<"Are they match the info for this function: total_shots: " << total_shots << ", distance: " <<
             distance << ", physcial_error: " << physcial_error << ", meas_er:" <<  meas_er <<"?[y/n]";
        std::string answer;
        std::getline(std::cin, answer);

        for (int i = 1; i < world_size; i++) {
            MPI_Send(answer.c_str(), answer.size() + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        } 

        if (answer[0] == 'N' || answer[0] == 'n') {
            return;
        } 
    } else {
        // Receive the answer from process 0
        char answer[2];
        MPI_Recv(answer, 2, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // If answer is 'N' or 'n', terminate program
        if (answer[0] == 'N' || answer[0] == 'n') {
            return;
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if(abort_f)
        return;
    syndromes_directory += "aggregate.txt";
    buffer = qpd::read_vector_of_vectors(syndromes_directory);
  
    // Broadcast the size of the buffer to all processes
    int buffer_size = buffer.size();
    MPI_Bcast(&buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize the buffer on all processes to match the size on process 0
    buffer.resize(buffer_size);

    // Broadcast the buffer data to all processes
    for (int i = 0; i < buffer_size; i++) {
        int vector_size = buffer[i].size();
        MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        buffer[i].resize(vector_size);
        MPI_Bcast(buffer[i].data(), vector_size, MPI_UINT16_T, 0, MPI_COMM_WORLD);
    }

    uint64_t start, end, syndromes_per_proc, num_syndromes;
    num_syndromes = ((uint64_t)buffer_size);
    syndromes_per_proc = num_syndromes/world_size;
    start = world_rank * syndromes_per_proc; 
    end = (world_rank == world_size - 1) ? num_syndromes : (world_rank + 1) * syndromes_per_proc;
    if(world_rank == 0 && print_logs){
        std::cout << "Total Number of HHWs: " << num_syndromes << std::endl;
    }
    const stim::Circuit circ = qrc::build_circuit(
            distance,
            0,
            0,
            true,
            true,
            false,
            0,
            distance,
            physcial_error,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            meas_er //measurement error rates
        );

    std::array<uint64_t, arr_N> hw_after_prematch;
    std::array<uint64_t, arr_N> hw_after_prematch_acc;

    hw_after_prematch.fill(0);
    hw_after_prematch_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hhw;
    number_of_hhw.fill(0);

    std::array<uint64_t, arr_N> number_of_hhw_acc;
    number_of_hhw_acc.fill(0);

    qpd::Predecoder predecoder(circ);
    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    const uint n_observables = circ.count_observables();

    uint64_t n_predecoder_decoder_logical_error = 0;
    uint64_t n_decoder_logcial_error = 0;
    uint64_t n_better_with_predec = 0;
    uint64_t n_better_with_predec_acc = 0;
    uint64_t n_predecoder_decoder_logical_error_acc = 0;
    uint64_t n_decoder_logcial_error_acc = 0;
    uint64_t n_remained_hhw = 0;
    uint64_t n_remained_hhw_acc = 0;
    bool predecoder_decoder_error = false;
    bool decoder_error = false;

    //std :: cout << "HERE"<< std::endl;
    for (int i = start; i < end; i++){
        //std::cout << "6" << std::endl;
        if(i%2'000'000 == 0 && print_logs && world_rank == 0){
            std::cout << world_rank << ": " << i << " out of " << end-start<< std::endl;
        }
        //std::cout << "Going to decompress" << std::endl;
        std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(buffer[i], distance);
        // std::cout << "next syndrome?" << std::endl;
        //qpd::print_syndrome(syndrome);
        // Calculate the hamming weights
        //std::cout << "5" << std::endl;
        uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
        // std::cout << "11" << std::endl;
        //predecoder.print_subgraph(syndrome);
        number_of_hhw[hw]++;
        // std::cout << "12" << std::endl;
        //Predecoder + MWPM decoder
        auto prematch_results = predecoder.matching_ec_length_one(syndrome, adding_boundry);
        // std::cout << "1, "  << prematch_results.pre_matches.size() << std::endl;
        auto predecoder_results = predecoder.predecode_error(prematch_results);
        // std::cout << "3" << std::endl;
        auto predecoder_decoder_results = decoder.decode_error(predecoder_results.post_syndrome);
        hw = std::accumulate(predecoder_results.post_syndrome.begin(), predecoder_results.post_syndrome.begin() + n_detectors, 0);
        hw_after_prematch[hw] ++;
        if( 10 < hw)
            n_remained_hhw++;

        std::vector<uint8_t> final_correction = predecoder_decoder_results.correction;
        // std::cout << "4" << " n_observables" << n_observables << " " << circ.count_observables() << " n_detectors" << n_detectors<< 
        // "  "  << circ.count_detectors() << "-- "<< syndrome.size()<< " CSIZE"<<final_correction.size() << std::endl;
        for (uint i = 0; i < n_observables; i++) {
            final_correction[i] ^= predecoder_results.correction[i];
        }
        //std::cout << "5" << std::endl;
        predecoder_decoder_error = qrc::is_logical_error(final_correction,syndrome,n_detectors, n_observables);
        if(predecoder_decoder_error)
            n_predecoder_decoder_logical_error++;
        
        // MWPM Decoder
        //std::cout << "6" << std::endl;
        // auto decoder_results = decoder.decode_error(syndrome);
        // //std::cout << "7" << std::endl;
        // decoder_error = qrc::is_logical_error(decoder_results.correction,syndrome,n_detectors, n_observables);
        // //std::cout << "8" << std::endl;
        // if(decoder_error)
        //     n_decoder_logcial_error++;
        //std::cout << "Just before condition " << std::endl;
        // if(decoder_error==false && predecoder_decoder_error==true){
        //     uint number_of_degree_1 = predecoder.count_degree_one();
        //     predecoder.print_shot_info(syndrome, prematch_results, 
        //                     predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
        //     //std::cout << "end of the if " << std::endl;
        // }
        //std::cout << "Before starting next round " << std::endl;
        if(decoder_error==true && predecoder_decoder_error==false){
            n_better_with_predec++;
            // std::cout << "**********BETTER MATCHING WITH PREDECODER!*********"<<std::endl;
            // uint number_of_degree_1 = predecoder.count_degree_one();
            // predecoder.print_shot_info(syndrome, prematch_results, 
            //                predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
            // std::cout << "***************************************************"<<std::endl;

        }

            
    
    }

    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        hw_after_prematch_acc[i] = sum;

        sum = 0;
        MPI_Reduce(&number_of_hhw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hhw_acc[i] = sum;
    }

    uint64_t sum = 0;
    MPI_Reduce(&n_predecoder_decoder_logical_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_predecoder_decoder_logical_error_acc = sum;
    sum = 0;
    MPI_Reduce(&n_remained_hhw, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_remained_hhw_acc = sum;

    // MPI_Reduce(&n_decoder_logcial_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    // n_decoder_logcial_error_acc = sum;
    sum = 0;
    MPI_Reduce(&n_better_with_predec, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_better_with_predec_acc = sum;
    //std::cout << "Before gathter" << std:: endl;
    // MPI_Gather(output_str.c_str(), output_str.size() + 1, MPI_CHAR, 
    //        output_strings.data(), output_str.size() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    //std::cout << "after" << std:: endl;

    if (world_rank == 0) {
        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hhw_acc[i];

        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "Astrea LER: " <<((double)t/total_shots) <<" - Number of HHW Syndromes: " << t << std::endl;
        out << "Predecoder + MWPM Decoder LER: " << ((double)n_predecoder_decoder_logical_error_acc/total_shots) << std::endl;
        out << "Predecoder + Astrea LER: " << ((double)n_remained_hhw_acc/total_shots) << " - Number of Remained HHW Syndromes: " << n_remained_hhw_acc << std::endl;
        //out << "MWPM Decoder # logical error: " << n_decoder_logcial_error_acc << std::endl;
        out << "# Times MWPM fails but Predecoder + MWPM deos not: " << n_better_with_predec_acc << std::endl;
        out << "***********************************************************" <<std::endl;
        out << "HW | #syndromes with that original HW | HW Frequency | #syndromes with that HW after prematching | After prematching HW Frequency \n";
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hhw_acc[i] || hw_after_prematch_acc[i]) {
                out << i << " | " << number_of_hhw_acc[i] << " | "<<
                ((double)number_of_hhw_acc[i]/t) << " | "<< hw_after_prematch_acc[i] << " | "<<
                ((double)hw_after_prematch_acc[i]/t) << std::endl;
            }
        }

    }

}

void predecoding_hhw_syndromes_astreaG(std::ofstream& out, std::string syndromes_directory, uint64_t total_shots, uint distance, 
        fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry){

    //data_folder_name += ("syndromes/" + number_of_shots_to_string(total_shots));
    
    // for (const auto & entry : std::filesystem::directory_iterator(data_folder_name))
    //     synd_files.push_back(entry.path());

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    bool abort_f = false;   
    std::vector<std::vector<uint16_t>> buffer;

    // This part is for checking if we are running correct experiments
    // based on the info of the main and based on the directory of the data.
    if(world_rank == 0){
        
        std::cout << "Running predecoding_hhw_syndromes_astreaG: Reading from " << syndromes_directory<<std::endl;
        // while ((dp = readdir(dirp)) != NULL) {
        //     if (dp->d_type == DT_REG) { // If the entry is a regular file
        //         std::string file_name = dp->d_name;
        //         std::string file_path = data_folder_name + "/" + file_name;
        //         synd_files.push_back(file_path);
        //     }
        // }
        // closedir(dirp);


        DIR* dirp = opendir(syndromes_directory.c_str());
        if (dirp == nullptr) {
            std::cout << "Error opening directory" << std::endl;
            abort_f = true;
        }
        std::cout << "These are two of the files in the directory " << syndromes_directory << ": ";
        int count = 0;
        struct dirent* dp;
        while ((dp = readdir(dirp)) != nullptr) {
            if (count == 2) { // Found two files, break out of the loop
                break;
            }

            if (dp->d_type == DT_REG && std::strstr(dp->d_name, ".txt") != nullptr) { // Check if the file is a regular file with a .txt extension
                std::cout << dp->d_name << std::endl;
                count++;
            }
        }
        closedir(dirp);

        std::cout <<"Are they match the info for this function: total_shots: " << total_shots << ", distance: " <<
             distance << ", physcial_error: " << physcial_error << ", meas_er:" <<  meas_er <<"?[y/n]";
        std::string answer;
        std::getline(std::cin, answer);

        for (int i = 1; i < world_size; i++) {
            MPI_Send(answer.c_str(), answer.size() + 1, MPI_CHAR, i, 0, MPI_COMM_WORLD);
        } 

        if (answer[0] == 'N' || answer[0] == 'n') {
            return;
        } 
    } else {
        // Receive the answer from process 0
        char answer[2];
        MPI_Recv(answer, 2, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // If answer is 'N' or 'n', terminate program
        if (answer[0] == 'N' || answer[0] == 'n') {
            return;
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    if(abort_f)
        return;
    syndromes_directory += "/aggregate.txt";
    buffer = qpd::read_vector_of_vectors(syndromes_directory);
  
    // Broadcast the size of the buffer to all processes
    int buffer_size = buffer.size();
    MPI_Bcast(&buffer_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Resize the buffer on all processes to match the size on process 0
    buffer.resize(buffer_size);

    // Broadcast the buffer data to all processes
    for (int i = 0; i < buffer_size; i++) {
        int vector_size = buffer[i].size();
        MPI_Bcast(&vector_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        buffer[i].resize(vector_size);
        MPI_Bcast(buffer[i].data(), vector_size, MPI_UINT16_T, 0, MPI_COMM_WORLD);
    }

    uint64_t start, end, syndromes_per_proc, num_syndromes;
    num_syndromes = ((uint64_t)buffer_size);
    syndromes_per_proc = num_syndromes/world_size;
    start = world_rank * syndromes_per_proc; 
    end = (world_rank == world_size - 1) ? num_syndromes : (world_rank + 1) * syndromes_per_proc;
    if(world_rank == 0 && print_logs){
        std::cout << "Total Number of HHWs: " << num_syndromes << std::endl;
    }
    const stim::Circuit circ = qrc::build_circuit(
            distance,
            0,
            0,
            true,
            true,
            false,
            0,
            distance,
            physcial_error,
            0,
            0,
            0,
            -1,
            -1,
            -1,
            meas_er //measurement error rates
        );

    std::array<uint64_t, arr_N> hw_after_prematch;
    std::array<uint64_t, arr_N> hw_after_prematch_acc;

    hw_after_prematch.fill(0);
    hw_after_prematch_acc.fill(0);

    std::array<uint64_t, arr_N> number_of_hhw;
    number_of_hhw.fill(0);

    std::array<uint64_t, arr_N> number_of_hhw_acc;
    number_of_hhw_acc.fill(0);
    qrc::AstreaParams astreaG_param= {};
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -1*log10(0.03*pow((physcial_error/0.0054), (((double)(distance+1)/2)))); // -1*log10(0.01*1e-9);
    qpd::Predecoder predecoder(circ);
    qrc::Astrea decoder(circ,
        n_detectors_per_round,
        weight_filter_cutoff,
        astreaG_param);

    std::mt19937_64 rng(world_rank);

    const uint n_detectors = circ.count_detectors();
    const uint n_observables = circ.count_observables();

    uint64_t n_predecoder_decoder_logical_error = 0;
    uint64_t n_decoder_logcial_error = 0;
    uint64_t n_better_with_predec = 0;
    uint64_t n_better_with_predec_acc = 0;
    uint64_t n_predecoder_decoder_logical_error_acc = 0;
    uint64_t n_decoder_logcial_error_acc = 0;
    uint64_t n_remained_hhw = 0;
    uint64_t n_remained_hhw_acc = 0;
    bool predecoder_decoder_error = false;
    bool decoder_error = false;

    //std :: cout << "HERE"<< std::endl;
    for (int i = start; i < end; i++){
        //std::cout << "6" << std::endl;
        if(i%5000 == 0 && print_logs && world_rank == 0){
            std::cout << world_rank << ": " << i << " out of " << end-start<< std::endl;
        }
        //std::cout << "Going to decompress" << std::endl;
        std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(buffer[i], distance);
        // std::cout << "next syndrome?" << std::endl;
        //qpd::print_syndrome(syndrome);
        // Calculate the hamming weights
        //std::cout << "5" << std::endl;
        uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
        // std::cout << "11" << std::endl;
        //predecoder.print_subgraph(syndrome);
        number_of_hhw[hw]++;
        // std::cout << "12" << std::endl;
        //Predecoder + MWPM decoder
        auto prematch_results = predecoder.matching_ec_length_one(syndrome, adding_boundry);
        // std::cout << "1, "  << prematch_results.pre_matches.size() << std::endl;
        auto predecoder_results = predecoder.predecode_error(prematch_results);
        // std::cout << "3" << std::endl;
        auto predecoder_decoder_results = decoder.decode_error(predecoder_results.post_syndrome);
        hw = std::accumulate(predecoder_results.post_syndrome.begin(), predecoder_results.post_syndrome.begin() + n_detectors, 0);
        hw_after_prematch[hw] ++;
        if( 10 < hw)
            n_remained_hhw++;

        std::vector<uint8_t> final_correction = predecoder_decoder_results.correction;
        // std::cout << "4" << " n_observables" << n_observables << " " << circ.count_observables() << " n_detectors" << n_detectors<< 
        // "  "  << circ.count_detectors() << "-- "<< syndrome.size()<< " CSIZE"<<final_correction.size() << std::endl;
        for (uint i = 0; i < n_observables; i++) {
            final_correction[i] ^= predecoder_results.correction[i];
        }
        //std::cout << "5" << std::endl;
        predecoder_decoder_error = qrc::is_logical_error(final_correction,syndrome,n_detectors, n_observables);
        if(predecoder_decoder_error){
            n_predecoder_decoder_logical_error++;
            std::cout << "Could not decode hw of " << hw << std::endl;
        }
        // MWPM Decoder
        // std::cout << "6" << std::endl;
        auto decoder_results = decoder.decode_error(syndrome);
        //std::cout << "7" << std::endl;
        // decoder_error = qrc::is_logical_error(decoder_results.correction,syndrome,n_detectors, n_observables);
        //std::cout << "8" << std::endl;
        if(decoder_results.is_logical_error)
            n_decoder_logcial_error++;
        if(predecoder_decoder_error){
            uint ex_hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
            std::cout << "Could not decode hw of " <<ex_hw<<","<< hw << std::endl;
            uint number_of_degree_1 = predecoder.count_degree_one();
            predecoder.print_shot_info(syndrome, prematch_results, 
                            predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
            
        }
        // std::cout << "Just before condition " << std::endl;
        if(decoder_error==false && predecoder_decoder_error==true){
            std::cout << "((()))"<< std::endl;
            uint number_of_degree_1 = predecoder.count_degree_one();
            predecoder.print_shot_info(syndrome, prematch_results, 
                            predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
            //std::cout << "end of the if " << std::endl;
        }
        // // std::cout << "Before starting next round " << std::endl;
        // if(decoder_error==true && predecoder_decoder_error==false){
        //     n_better_with_predec++;
        //     // std::cout << "**********BETTER MATCHING WITH PREDECODER!*********"<<std::endl;
        //     // uint number_of_degree_1 = predecoder.count_degree_one();
        //     // predecoder.print_shot_info(syndrome, prematch_results, 
        //     //                predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
        //     // std::cout << "***************************************************"<<std::endl;

        // }

            
    
    }

    for ( uint i = 0; i< arr_N;  i++){
        uint64_t sum = 0;
        MPI_Reduce(&hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        hw_after_prematch_acc[i] = sum;

        sum = 0;
        MPI_Reduce(&number_of_hhw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        number_of_hhw_acc[i] = sum;
    }

    uint64_t sum = 0;
    MPI_Reduce(&n_predecoder_decoder_logical_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_predecoder_decoder_logical_error_acc = sum;
    sum = 0;
    MPI_Reduce(&n_remained_hhw, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_remained_hhw_acc = sum;
    sum = 0;
    MPI_Reduce(&n_decoder_logcial_error, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_decoder_logcial_error_acc = sum;
    sum = 0;
    MPI_Reduce(&n_better_with_predec, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    n_better_with_predec_acc = sum;
    //std::cout << "Before gathter" << std:: endl;
    // MPI_Gather(output_str.c_str(), output_str.size() + 1, MPI_CHAR, 
    //        output_strings.data(), output_str.size() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    //std::cout << "after" << std:: endl;

    if (world_rank == 0) {
        uint64_t t = 0;
        for (uint i=0; i<arr_N; i++)
            t+=number_of_hhw_acc[i];

        out << "Distance: "<< distance << " Total Shots: " << total_shots << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
        out << "Astrea LER: " <<((double)t/total_shots) <<" - Number of HHW Syndromes: " << t << std::endl;
        out << "Predecoder + Astrea LER: " << ((double)n_remained_hhw_acc/total_shots) << " - Number of Remained HHW Syndromes: " << n_remained_hhw_acc << std::endl;
        out << "Astrea-G Decoder LER: " <<  ((double)n_decoder_logcial_error_acc/total_shots) << std::endl;
        out << "Predecoder + Astrea-G Decoder LER: " << ((double)n_predecoder_decoder_logical_error_acc/total_shots) << std::endl;        //out << "# Times MWPM fails but Predecoder + MWPM deos not: " << n_better_with_predec_acc << std::endl;
        out << "***********************************************************" <<std::endl;
        out << "HW | #syndromes with that original HW | HW Frequency | #syndromes with that HW after prematching | After prematching HW Frequency \n";
        for (uint i = 0; i < arr_N; i++) {
            if (number_of_hhw_acc[i] || hw_after_prematch_acc[i]) {
                out << i << " | " << number_of_hhw_acc[i] << " | "<<
                ((double)number_of_hhw_acc[i]/t) << " | "<< hw_after_prematch_acc[i] << " | "<<
                ((double)hw_after_prematch_acc[i]/t) << std::endl;
            }
        }

    }

}

uint64_t parseMagnitude(const std::string& input){
    std::unordered_map<char, uint64_t> magnitudeMap = {
        {'K', 1000},
        {'M', 1000000},
        {'B', 1000000000},
    };
    uint64_t value = std::stoull(input);
    char magnitudeChar = input.back();

    if (magnitudeMap.find(magnitudeChar) != magnitudeMap.end()) {
        uint64_t magnitude = magnitudeMap[magnitudeChar];
        value *= magnitude;
    }

    return value;
}



