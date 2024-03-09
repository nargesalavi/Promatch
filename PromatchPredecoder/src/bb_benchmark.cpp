#include "bb_benchmark.h"

fp_t poisson_probability(fp_t lambda, int k){
    fp_t p = std::exp(-lambda);
    for (uint i = 1; i <= k; ++i) {
        p *= lambda / i;
    }
    return p;
}

fp_t binomial_probability(uint k, uint n, fp_t p) {
    fp_t q = 1.0 - p;
    fp_t prob = std::exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)) * std::pow(p, k) * std::pow(q, n-k);
    return prob;
}



void test_sim(std::ofstream& out, std::string syndromes_directory, uint64_t total_shots, uint distance, 
fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry){
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
    uint k_max = 50;
    bbsim::BucketBasedSim simulator(circ,k_max);
    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, 500'000'000);
    std::cout << "expected_ler = " << expected_ler << " physcial_error = " << physcial_error << std::endl;
    for(uint i = 0 ; i< simulator.shots_per_k.size(); i++){
        std::cout << i << " - " << simulator.prob_k[i] << " - "<< 
        simulator.shots_per_k[i] << std::endl;
    }
    // simulator.print_decoding_graph();
    // simulator.print_sorted_event_vector();
    // simulator.print_vertices();
    simulator.print_path_list();
    // simulator.print_path_matrix();
    simulator.print_weight_matrix();
    // std::array<std::array<fp_t,721>, 721> m;
    // m.fill(std::array<fp_t,721>{0});
    // auto path_table = qrc::compute_path_table(simulator.decoding_graph);
    // uint i,j;
    // for (auto v :simulator.decoding_graph.vertices()){
    //     for(auto w: simulator.decoding_graph.vertices()){
    //         auto v_w = std::make_pair(v,w);
    //         auto p = path_table[v_w].distance;
    //         i = v->detector;
    //         j=w->detector;
    //         if(v->detector == BOUNDARY_INDEX){
    //             i = 720;
    //         }
    //         if(w->detector == BOUNDARY_INDEX){
    //             j = 720;
    //         }
    //         // std::cout << i << " "<<j << std::endl;
    //         fp_t x = pow(10,-1*p);
    //         m[i][j] = x;;;
    //     }
    // }


    // for(uint i = 0; i<721; i++){
    //     for(uint j = 0; j< 721; j++){
    //         std::cout<< m[i][j] << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "Distance = " << distance << std::endl;
    // std::cout << "Size of edges = " << simulator.list_of_error_probabs.size() << std::endl;
    // for(uint i = 0; i<k_max+1; i++){
    //     std::cout << "K = " << i << " :  p = " << simulator.prob_k[i] << std::endl;
    // }
    
}

void bb_ler_calculation_and_stats(std::ofstream& out, std::string syndromes_directory, uint64_t shots_per_k, uint distance, 
fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry){
    uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    
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
    const uint n_buckets = 20;
    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    bbsim::BucketBasedSim simulator(circ,n_buckets);
    auto path_table = qrc::compute_path_table(simulator.decoding_graph);

    qpd::Predecoder predecoder(circ);
    qrc::MWPMDecoder decoder(circ);
    std::mt19937_64 rng(world_rank);

    
    bool predecoder_decoder_error = false;
    std::array<uint64_t, n_buckets+1> n_predecoder_decoder_logical_error;
    std::array<uint64_t, n_buckets+1> n_decoder_logcial_error;
    n_predecoder_decoder_logical_error.fill(0);
    n_decoder_logcial_error.fill(0);
    fp_t mwpm_ler = 0;
    fp_t p_mwpm_ler = 0;
    fp_t astrea_ler = 0;
    fp_t p_astrea_ler = 0;
    uint n_observables = circ.count_observables();
    uint n_detectors = circ.count_detectors();

    std::array<uint64_t, n_buckets+1> n_predecoder_decoder_logical_error_acc;
    std::array<uint64_t, n_buckets+1> n_decoder_logcial_error_acc;
    n_predecoder_decoder_logical_error_acc.fill(0);
    n_decoder_logcial_error_acc.fill(0);

    std::array<uint64_t, n_buckets+1> n_predecoder_astrea_logical_error;
    std::array<uint64_t, n_buckets+1> n_predecoder_astrea_logical_error_acc;
    n_predecoder_astrea_logical_error.fill(0);
    n_predecoder_astrea_logical_error_acc.fill(0);
    ///#######
    std::array<uint64_t, n_buckets+1> hhw_er;
    std::array<uint64_t, n_buckets+1> hhw_er_acc;
    hhw_er.fill(0);
    hhw_er_acc.fill(0);

    std::array<uint64_t, n_buckets+1> lhw_er;
    std::array<uint64_t, n_buckets+1> lhw_er_acc;
    lhw_er.fill(0);
    lhw_er_acc.fill(0);

    std::array<uint64_t, n_buckets+1> hhw_er0;
    std::array<uint64_t, n_buckets+1> hhw_er0_acc;
    hhw_er0.fill(0);
    hhw_er0_acc.fill(0);

    std::array<uint64_t, n_buckets+1> lhw_er0;
    std::array<uint64_t, n_buckets+1> lhw_er0_acc;
    lhw_er0.fill(0);
    lhw_er0_acc.fill(0);
    ///#######

    std::array<uint64_t, n_buckets+1> n_astrea_logical_error;
    std::array<uint64_t, n_buckets+1> n_astrea_logical_error_acc;
    n_astrea_logical_error.fill(0);
    n_astrea_logical_error_acc.fill(0);

    std::array<uint64_t, n_buckets+1> t;
    std::array<uint64_t, n_buckets+1> t_acc;
    t.fill(0);
    t_acc.fill(0);

    std::array<uint64_t, n_buckets+1> t_hhw;
    std::array<uint64_t, n_buckets+1> t_hhw_acc;
    t_hhw.fill(0);
    t_hhw_acc.fill(0);

    std::array<uint64_t, n_buckets+1> n_hhw;
    std::array<uint64_t, n_buckets+1> n_hhw_acc;
    n_hhw.fill(0);
    n_hhw_acc.fill(0);

    std::vector<uint8_t> syndrome;
    uint max_effective_k = 0;

    std::array<uint64_t, n_buckets*2> hw_errors;
    std::array<uint64_t, n_buckets*2> hw_errors_acc;
    hw_errors.fill(0);
    hw_errors_acc.fill(0);

    std::array<std::array<uint64_t, n_buckets*2>,n_buckets*2> hw_errors_pred;
    std::array<std::array<uint64_t, n_buckets*2>,n_buckets*2>  hw_errors_pred_acc;
    hw_errors_pred.fill(std::array<uint64_t, n_buckets*2>{0});
    hw_errors_pred_acc.fill(std::array<uint64_t, n_buckets*2>{0});

    //std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1> ec_all_hw;
    std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1>  ec_only_hhw;
    //ec_all_hw.fill(std::array<uint64_t, n_buckets*2>{0});
    ec_only_hhw.fill(std::array<uint64_t, n_buckets*2>{0});

    //std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1> ec_all_hw_acc;
    std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1>  ec_only_hhw_acc;
    //ec_all_hw_acc.fill(std::array<uint64_t, n_buckets*2>{0});
    ec_only_hhw_acc.fill(std::array<uint64_t, n_buckets*2>{0});

    // std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1> hhws;
    // std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1> new_hws;
    // hhws.fill(std::array<uint64_t, n_buckets*2>{0});
    // new_hws.fill(std::array<uint64_t, n_buckets*2>{0});

    // std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1> hhws_acc;
    // std::array<std::array<uint64_t, n_buckets*2>,n_buckets+1> new_hws_acc;
    // hhws_acc.fill(std::array<uint64_t, n_buckets*2>{0});
    // new_hws_acc.fill(std::array<uint64_t, n_buckets*2>{0});

    bool predecoded = false;
    if(world_rank == 0){
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Shots per bucket = " << shots_per_k << std::endl;
    }
    
    for(uint k = floor(distance/2); k < n_buckets+1; k++){
        if(simulator.prob_k[k] < expected_ler){
            max_effective_k = k - 1;
            break;
        }
        uint64_t shots = shots_per_k / world_size;
        if (world_rank == world_size - 1) {
            shots += shots_per_k % world_size;
        }
        if(world_rank == 0 ){
                std::cout << "--k: " << k <<std::endl;//qpd::print_syndrome(syndrome);
        }
         
        while(shots > 0){
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            std::vector<std::vector<uint8_t>> buffer = simulator.create_syndromes( k, shots_this_round, rng, print_logs);


            if(world_rank == 0 ){
                std::cout << "k: " << k<<", "<<buffer.size()<< " " <<simulator.prob_k[k]<<std::endl;//qpd::print_syndrome(syndrome);
            }
            for(auto syndrome : buffer){
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
            predecoded = false;
            if(10 < hw){
                n_hhw[k]++;
                //hhws[k][hw]++;
                predecoded = true;
                auto prematch_results = predecoder.matching_ec_length_one(syndrome, adding_boundry);
                
                auto predecoder_results = predecoder.predecode_error(prematch_results);
                
                auto predecoder_decoder_results = decoder.decode_error(predecoder_results.post_syndrome);
                uint new_hw = std::accumulate(predecoder_results.post_syndrome.begin(), predecoder_results.post_syndrome.begin() + n_detectors, 0);
                
                std::vector<uint8_t> final_correction = predecoder_decoder_results.correction;
                //new_hws[k][new_hw]++;
                for (uint i = 0; i < n_observables; i++) {
                    final_correction[i] ^= predecoder_results.correction[i];
                }
                
                predecoder_decoder_error = qrc::is_logical_error(final_correction,syndrome,n_detectors, n_observables);
                if(predecoder_decoder_error){
                    n_predecoder_decoder_logical_error[k]++;
                    hw_errors_pred[hw][new_hw]++;
                }
                if( 10 < new_hw){
                    n_predecoder_astrea_logical_error[k]++;
                    hhw_er[k]++;
                }
                else{
                    if(predecoder_decoder_error){
                        n_predecoder_astrea_logical_error[k]++;
                        lhw_er[k]++;
                    }
                }
            }

            auto decoder_results = decoder.decode_error(syndrome);
            
            if(decoder_results.is_logical_error){
                n_decoder_logcial_error[k]++;
                hw_errors[hw]++;
                if(!predecoded){
                    n_predecoder_decoder_logical_error[k]++;
                }
            }
            if( 10 < hw){
                n_astrea_logical_error[k]++;
                hhw_er0[k]++;
            }
            else{
                if(decoder_results.is_logical_error){
                    n_astrea_logical_error[k]++;
                    lhw_er0[k]++;
                }
            }

            for (auto pair : decoder_results.matching) {
                auto v = simulator.decoding_graph.get_vertex(pair.first);
                auto w = simulator.decoding_graph.get_vertex(pair.second);
                auto v_w = std::make_pair(v, w);
                // Finding error chains from the path table of decoding graph 
                uint ec = path_table[v_w].path.size() - 1;
                // count the number of error chains. Each cell of array i in ec_all_hw
                // contains the number of error chains with length i. 
                t[k]++;
                //ec_all_hw[k][ec]++;
                if (hw > 10) {
                    t_hhw[k]++;
                    ec_only_hhw[k][ec]++;
                }

            }

            
            }
            shots -= shots_this_round;
        }
        if(n_decoder_logcial_error[k]){
                std::cout << "Processor " << world_rank << " found "<< n_decoder_logcial_error[k] << " in bucket " << k << "(prob_k = " << simulator.prob_k[k] << ") "<< std::endl;
            }
        
        // mwpm_ler += poisson_probability(simulator.return_probab_sum(), k)*((fp_t)n_decoder_logcial_error[k]/shots_per_k);
        // p_mwpm_ler += poisson_probability(simulator.return_probab_sum(), k)*((fp_t)n_predecoder_decoder_logical_error[k]/shots_per_k);
    }
    MPI_Reduce(&n_predecoder_decoder_logical_error, &n_predecoder_decoder_logical_error_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    // sum = 0;
    // MPI_Reduce(&n_remained_hhw, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    // n_remained_hhw_acc = sum;

    MPI_Reduce(&n_decoder_logcial_error, &n_decoder_logcial_error_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&n_predecoder_astrea_logical_error, &n_predecoder_astrea_logical_error_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&n_astrea_logical_error, &n_astrea_logical_error_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&hhw_er0, &hhw_er0_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&lhw_er0, &lhw_er0_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&hhw_er, &hhw_er_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&lhw_er, &lhw_er_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&hw_errors, &hw_errors_acc, n_buckets*2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        
    MPI_Reduce(&hw_errors_pred, &hw_errors_pred_acc, n_buckets*2*n_buckets*2, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    //MPI_Reduce(&ec_all_hw, &ec_all_hw_acc, n_buckets*2*(n_buckets+1), MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&ec_only_hhw, &ec_only_hhw_acc, n_buckets*2*(n_buckets+1), MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    //MPI_Reduce(&hhws, &hhws_acc, n_buckets*2*(n_buckets+1), MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    //MPI_Reduce(&new_hws, &new_hws_acc, n_buckets*2*(n_buckets+1), MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&t, &t_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_hhw, &t_hhw_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&n_hhw, &n_hhw_acc, n_buckets+1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);



    


    if (world_rank == 0) {
        uint64_t total_hhw_er0 = 0;
        uint64_t total_lhw_er0 = 0;
        uint64_t total_hhw_er = 0;
        uint64_t total_lhw_er = 0;
        uint64_t astrea_er = 0;
        uint64_t pred_astrea_er = 0;
        fp_t avg_hhw_er_fraction = 0.0;
        fp_t avg_lhw_er_fraction = 0.0;

        fp_t avg_hhw_er0_fraction = 0.0;
        fp_t avg_lhw_er0_fraction = 0.0;
        for ( uint k = 1; k< max_effective_k+1;  k++){
            // mwpm_ler += binomial_probability(k, simulator.number_of_events, ((fp_t) simulator.return_probab_sum()/simulator.number_of_events)) *((fp_t)n_decoder_logcial_error_acc[k]/shots_per_k);
            // p_mwpm_ler += binomial_probability(k, simulator.number_of_events, ((fp_t)simulator.return_probab_sum()/simulator.number_of_events)) *((fp_t)n_predecoder_decoder_logical_error_acc[k]/shots_per_k);
            // astrea_ler += binomial_probability(k, simulator.number_of_events, ((fp_t)simulator.return_probab_sum()/simulator.number_of_events)) *((fp_t)n_astrea_logical_error_acc[k]/shots_per_k);
            // p_astrea_ler += binomial_probability(k, simulator.number_of_events, ((fp_t)simulator.return_probab_sum()/simulator.number_of_events)) *((fp_t)n_predecoder_astrea_logical_error_acc[k]/shots_per_k);

            mwpm_ler +=  simulator.prob_k[k]*((fp_t)n_decoder_logcial_error_acc[k]/shots_per_k);
            p_mwpm_ler += simulator.prob_k[k]*((fp_t)n_predecoder_decoder_logical_error_acc[k]/shots_per_k);
            astrea_ler +=  simulator.prob_k[k]*((fp_t)n_astrea_logical_error_acc[k]/shots_per_k);
            p_astrea_ler += simulator.prob_k[k]*((fp_t)n_predecoder_astrea_logical_error_acc[k]/shots_per_k);

            total_hhw_er0+=hhw_er0_acc[k];
            total_hhw_er+=hhw_er_acc[k];
            total_lhw_er0+=lhw_er0_acc[k];
            total_lhw_er+=lhw_er_acc[k];
            astrea_er += n_astrea_logical_error_acc[k];
            pred_astrea_er += n_predecoder_astrea_logical_error_acc[k];

            if(n_predecoder_astrea_logical_error_acc[k])
                avg_hhw_er_fraction += ((fp_t)hhw_er_acc[k]/n_predecoder_astrea_logical_error_acc[k]);
            if(n_astrea_logical_error_acc[k])
                avg_hhw_er0_fraction += ((fp_t)hhw_er0_acc[k]/n_astrea_logical_error_acc[k]);
            if(n_predecoder_astrea_logical_error_acc[k])
                avg_lhw_er_fraction += ((fp_t)lhw_er_acc[k]/n_predecoder_astrea_logical_error_acc[k]);
            if(n_astrea_logical_error_acc[k])
                avg_lhw_er0_fraction += ((fp_t)lhw_er0_acc[k]/n_astrea_logical_error_acc[k]);


        }
        avg_hhw_er_fraction /= n_buckets;
        avg_hhw_er0_fraction /= n_buckets;
        avg_lhw_er_fraction /= n_buckets;
        avg_lhw_er0_fraction /= n_buckets;
        out << "Distance: "<< distance << " - #Buckets: " << n_buckets <<  " - Shots Per Bucket: " << shots_per_k <<" - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;

        out << "MWPM Decoder LER: " << mwpm_ler << std::endl;
        out << "Predecoder + MWPM Decoder LER: " << p_mwpm_ler << std::endl;
        out << "Astrea Decoder LER: " << astrea_ler << "(h: " <<avg_hhw_er0_fraction << ", l: "<< avg_lhw_er0_fraction<<")"<< std::endl;
        out << "Predecoder + Astrea Decoder LER: " << p_astrea_ler << "(h: " <<avg_hhw_er_fraction <<", l: "<< avg_lhw_er_fraction<<")" << std::endl;
        out << "***********************************************************" <<std::endl;
        out << " k | # of Err in MWPM - impact on LER | # of Err in Pre+MWPM - impact on LER || # of Err in Asrea - impact on LER | # of Err in Pre+Asrea - impact on LER| " << std::endl;

        for ( uint i = floor(distance/2); i< max_effective_k+1;  i++){
            out << i << " | "<< n_decoder_logcial_error_acc[i] << " - " << (simulator.prob_k[i]*((fp_t)n_decoder_logcial_error_acc[i]/shots_per_k))
            << " | " << n_predecoder_decoder_logical_error_acc[i] << " - " << simulator.prob_k[i]*((fp_t)n_predecoder_decoder_logical_error_acc[i]/shots_per_k) 
            << " | "<< n_astrea_logical_error_acc[i] << "(hhw: "<< hhw_er0_acc[i] <<", lhw: " << lhw_er0_acc[i]<<") - " << simulator.prob_k[i]*((fp_t)n_astrea_logical_error_acc[i]/shots_per_k)
            << " | " << n_predecoder_astrea_logical_error_acc[i] << "(hhw: "<< hhw_er_acc[i] <<", lhw: " << lhw_er_acc[i]<<") - " << simulator.prob_k[i]*((fp_t)n_predecoder_astrea_logical_error_acc[i]/shots_per_k) << std::endl;
        }

        out << "________________________________________________________________" <<std::endl;
        out << "HW | MWPM Error | Pred+MWPM Error " << std::endl;

        for(uint i = 0; i< n_buckets*2; i++){
            if(hw_errors_acc[i] !=0 || std::accumulate(hw_errors_pred_acc[i].begin(), hw_errors_pred_acc[i].begin() + n_buckets*2, 0) !=0)
                out << i << " | "<< hw_errors_acc[i] << " | " << std::accumulate(hw_errors_pred_acc[i].begin(), hw_errors_pred_acc[i].begin() + n_buckets*2, 0) << std::endl;
        }

        //out << "###############################MORE DETAILS#####################################" <<std::endl;
        // out << "HW | MWPM Error " << std::endl;
        // out << "    New HW | Pred+MWPM Error " << std::endl;
        // for(uint i = 0; i< n_buckets*2; i++){
        //     if(hw_errors_acc[i] !=0 || std::accumulate(hw_errors_pred_acc[i].begin(), hw_errors_pred_acc[i].begin() + n_buckets*2, 0) !=0){
        //         out << i << " | "<< hw_errors_acc[i] << std::endl;
            
        //         for(uint j = 0; j<n_buckets*2; j++){
        //             if(hw_errors_pred_acc[i][j] != 0)
        //             out << "\t" << j<<" | " << hw_errors_pred_acc[i][j]<< std:: endl;
        //         }
        //     }
        // }
        std::array<uint64_t, n_buckets+1> total_ec;
        std::array<uint64_t, n_buckets+1> total_ec_hhw;
    
        for(uint i = 0; i<n_buckets+1; i++){

            // total_ec[i] = std::accumulate(ec_all_hw_acc[i].begin(), ec_all_hw_acc[i].begin() + n_buckets*2, 0);
            total_ec_hhw[i] = std::accumulate(ec_only_hhw_acc[i].begin(), ec_only_hhw_acc[i].begin() + n_buckets*2, 0);
        }

        // out << "________________________________________________________________" <<std::endl;
        // out << "EC | Number of that EC  | Frequency" << std::endl;
        // for(uint k = 0; k < n_buckets; k++){
        //     if(total_ec[k] == 0)
        //         continue;
        //     out << " K = " << k << " | Probability of Bucket = " << simulator.prob_k[k] << " | Total HHW in Bucket = " << total_ec[k] << std::endl;

        //     for(uint i = 0; i< n_buckets*2; i++){
        //         if(ec_all_hw_acc[k][i] !=0){
        //             out << i << " | "<< ec_all_hw_acc[k][i] << " | " << ((double)ec_all_hw_acc[k][i]/total_ec[k])*simulator.prob_k[k] << std::endl;
        //         }
        //     }
        // }

        out << "________________________________________________________________" <<std::endl;
        out << "EC for HHW | Number of that EC in HHW Syndromes| Frequency" << std::endl;
        for(uint k = 0; k < n_buckets; k++){
            if(total_ec_hhw[k] == 0)
                continue;
            out << " K = " << k << " | Probability of Bucket = " << simulator.prob_k[k] << " | Total HHW in Bucket = " << total_ec_hhw[k] <<std::endl;

            for(uint i = 0; i< n_buckets*2; i++){
                if(ec_only_hhw_acc[k][i] !=0){
                    out << i << " | "<< ec_only_hhw_acc[k][i] << " | " << ((double)ec_only_hhw_acc[k][i]/total_ec_hhw[k])*simulator.prob_k[k] << std::endl;
                }
            }
        }


        // out << "________________________________________________________________" <<std::endl;
        // out << "Info about HWs before and after predecoding "<< std::endl;
        // for(uint k = 0; k < n_buckets; k++){
        //     if(t_acc[k] == 0)
        //         continue;
        //     out << " K = " << k << " | Probability of Bucket = " << simulator.prob_k[k] << std::endl;
        //     out << "HHW | Frequency" << std::endl;

        //     for(uint i = 0; i< n_buckets*2; i++){
        //         if(hhws_acc[k][i] !=0){
        //             out << i << " | "<< ((double)hhws_acc[k][i]/n_hhw_acc[k])*simulator.prob_k[k]  << std::endl;
        //         }
        //     }

        //     out << "New HW | Frequency" << std::endl;

        //     for(uint i = 0; i< n_buckets*2; i++){
        //         if(new_hws_acc[k][i] !=0){
        //             out << i << " | "<< ((double)new_hws_acc[k][i]/n_hhw_acc[k])*simulator.prob_k[k]  << std::endl;
        //         }
        //     }
        // }
        std::cout << "Distance = " << distance << " LER = " << mwpm_ler << std::endl;
        

    }

}


void bb_ler_calculation_AstreaG(std::ofstream& out, std::string syndromes_directory, uint distance, 
fp_t physcial_error, fp_t meas_er, bool print_logs, bool adding_boundry){
    uint64_t shots_per_k = 12'000'000;
    uint64_t SHOTS_PER_BATCH = 1'000'000;

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    
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
    const uint n_buckets = 20;
    bbsim::BucketBasedSim simulator(circ,n_buckets);
    ;

    qpd::Predecoder predecoder(circ);

    qrc::AstreaParams astreaG_param= {};
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -1*log10(0.03*pow((physcial_error/0.0054), (((double)(distance+1)/2)))); // -1*log10(0.01*1e-9);
    qrc::Astrea decoder(circ,
        n_detectors_per_round,
        weight_filter_cutoff,
        astreaG_param);
        
    std::mt19937_64 rng(world_rank);

    
    bool predecoder_decoder_error = false;
    std::array<uint64_t, n_buckets+1> n_predecoder_decoder_logical_error;
    std::array<uint64_t, n_buckets+1> n_decoder_logcial_error;
    n_predecoder_decoder_logical_error.fill(0);
    n_decoder_logcial_error.fill(0);
    fp_t mwpm_ler = 0;
    fp_t p_mwpm_ler = 0;
    fp_t astrea_ler = 0;
    fp_t p_astrea_ler = 0;
    uint n_observables = circ.count_observables();
    uint n_detectors = circ.count_detectors();

    std::array<uint64_t, n_buckets+1> n_predecoder_decoder_logical_error_acc;
    std::array<uint64_t, n_buckets+1> n_decoder_logcial_error_acc;
    n_predecoder_decoder_logical_error_acc.fill(0);
    n_decoder_logcial_error_acc.fill(0);

    std::vector<uint8_t> syndrome;

    for(uint k = 1; k < n_buckets+1; k++){
        uint64_t shots = shots_per_k / world_size;
        if (world_rank == world_size - 1) {
            shots += shots_per_k % world_size;
        }
        if(world_rank == 0 ){
                std::cout << "--k: " << k <<std::endl;//qpd::print_syndrome(syndrome);
        }
         
        while(shots > 0){
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            std::vector<std::vector<uint8_t>> buffer = simulator.create_syndromes( k, shots_this_round, rng, print_logs);


           // if(world_rank == 0 ){
                 // std::cout << "k: " << k<< ":"<< n_decoder_logcial_error[k]<<", "<<buffer.size()<< " " <<binomial_probability(k, simulator.number_of_events, ((fp_t)simulator.return_probab_sum()/simulator.number_of_events))<<std::endl;//qpd::print_syndrome(syndrome);
		
            //}
            for(auto syndrome : buffer){
            //syndrome = buffer[s];
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
            // std::cout << "11" << std::endl;
            
            // predecoder.print_subgraph(syndrome);
            // number_of_hhw[hw]++;
            // std::cout << "12" << std::endl;
            //Predecoder + MWPM decoder
            auto prematch_results = predecoder.matching_ec_length_one(syndrome, adding_boundry);
            // std::cout << "1, "  << prematch_results.pre_matches.size() << std::endl;
            auto predecoder_results = predecoder.predecode_error(prematch_results);
            // std::cout << "3" << std::endl;
            auto predecoder_decoder_results = decoder.decode_error(predecoder_results.post_syndrome);
            uint new_hw = std::accumulate(predecoder_results.post_syndrome.begin(), predecoder_results.post_syndrome.begin() + n_detectors, 0);
            // hw_after_prematch[hw] ++;
            std::vector<uint8_t> final_correction = predecoder_decoder_results.correction;
            // std::cout << "4" << " n_observables" << n_observables << " " << circ.count_observables() << " n_detectors" << n_detectors<< 
            // "  "  << circ.count_detectors() << "-- "<< syndrome.size()<< " CSIZE"<<final_correction.size() << std::endl;
            for (uint i = 0; i < n_observables; i++) {
                final_correction[i] ^= predecoder_results.correction[i];
            }
            //std::cout << "5" << std::endl;
            predecoder_decoder_error = qrc::is_logical_error(final_correction,syndrome,n_detectors, n_observables);
            if(predecoder_decoder_error){
                n_predecoder_decoder_logical_error[k]++;
                // std::cout << " k : " <<k <<n_predecoder_decoder_logical_error[k]<< std::endl;
                // std::cout << "_______"<<std::endl;
                // std::cout << "Error" << hw << ", "<< new_hw<< std::endl;
                // qpd::print_syndrome(syndrome);
                // std::cout << "Size: "<< final_correction.size()<< ", xx: "<< ((int) final_correction[0])<<std::endl;
                // std::cout << "_______"<<std::endl;
            }
            
        
            // MWPM Decoder
            // std::cout << "6" << std::endl;
            auto decoder_results = decoder.decode_error(syndrome);
            //std::cout << "7" << std::endl;
            // decoder_error = qrc::is_logical_error(decoder_results.correction,syndrome,n_detectors, n_observables);
            //std::cout << "8" << std::endl;
            if(decoder_results.is_logical_error){
                n_decoder_logcial_error[k]++;
            }
            
            // if(predecoder_decoder_error){
            //     uint ex_hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
            //     std::cout << "Could not decode hw of " <<ex_hw<<","<< hw << std::endl;
            //     uint number_of_degree_1 = predecoder.count_degree_one();
            //     predecoder.print_shot_info(syndrome, prematch_results, 
            //                     predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
                
            // }
            // // std::cout << "Just before condition " << std::endl;
            // if(decoder_error==false && predecoder_decoder_error==true){
            //     std::cout << "((()))"<< std::endl;
            //     uint number_of_degree_1 = predecoder.count_degree_one();
            //     predecoder.print_shot_info(syndrome, prematch_results, 
            //                     predecoder_decoder_results, decoder_results, number_of_degree_1, false, false);
            //     //std::cout << "end of the if " << std::endl;
            // }
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
            shots -= shots_this_round;
        }
        
        // mwpm_ler += poisson_probability(simulator.return_probab_sum(), k)*((fp_t)n_decoder_logcial_error[k]/shots_per_k);
        // p_mwpm_ler += poisson_probability(simulator.return_probab_sum(), k)*((fp_t)n_predecoder_decoder_logical_error[k]/shots_per_k);
    }

    // for ( uint i = 0; i< arr_N;  i++){
    //     uint64_t sum = 0;
    //     MPI_Reduce(&hw_after_prematch[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    //     hw_after_prematch_acc[i] = sum;

    //     sum = 0;
    //     MPI_Reduce(&number_of_hhw[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    //     number_of_hhw_acc[i] = sum;
    // }
    uint64_t sum = 0;
    for ( uint i = 1; i< n_buckets+1;  i++){
        sum = 0;
        MPI_Reduce(&n_predecoder_decoder_logical_error[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        n_predecoder_decoder_logical_error_acc[i] = sum;
        // sum = 0;
        // MPI_Reduce(&n_remained_hhw, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        // n_remained_hhw_acc = sum;
        sum = 0;
        MPI_Reduce(&n_decoder_logcial_error[i], &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        n_decoder_logcial_error_acc[i] = sum;

    }
    
    // sum = 0;
    // MPI_Reduce(&n_better_with_predec, &sum, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    // n_better_with_predec_acc = sum;
    //std::cout << "Before gathter" << std:: endl;
    // MPI_Gather(output_str.c_str(), output_str.size() + 1, MPI_CHAR, 
    //        output_strings.data(), output_str.size() + 1, MPI_CHAR, 0, MPI_COMM_WORLD);
    //std::cout << "after" << std:: endl;

    if (world_rank == 0) {
        for ( uint k = 1; k< n_buckets+1;  k++){
            mwpm_ler += binomial_probability(k, simulator.number_of_events, ((fp_t) simulator.return_probab_sum()/simulator.number_of_events)) *((fp_t)n_decoder_logcial_error_acc[k]/shots_per_k);
            p_mwpm_ler += binomial_probability(k, simulator.number_of_events, ((fp_t)simulator.return_probab_sum()/simulator.number_of_events)) *((fp_t)n_predecoder_decoder_logical_error_acc[k]/shots_per_k);
        }
        out << "Distance: "<< distance << " - #Buckets: " << n_buckets <<  " - Shots Per Bucket: " << shots_per_k <<" - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;

        out << "AstreaG Decoder LER: " << mwpm_ler << std::endl;
        out << "Predecoder + AstreaG Decoder LER: " << p_mwpm_ler << std::endl;
        out << "***********************************************************" <<std::endl;
        for ( uint i = 1; i< n_buckets+1;  i++){
            out << " k | # of Err in AstreaG | # of Err in Pre+AstreaG  " << std::endl;
            out << i << " | "<< n_decoder_logcial_error_acc[i] << " | "<< n_predecoder_decoder_logical_error_acc[i] << std::endl;
        }


    }

}

void test_b_statistical_ler(uint64_t shots_per_k, uint distance, fp_t physcial_error, fp_t meas_er){

    // int world_size, world_rank;
    // MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    // MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);  

    
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
    const uint min_fault = floor(distance/2);
    qrc::MWPMDecoder* decoder =  new qrc::MWPMDecoder(circ);

    std::random_device rd;
    std::mt19937_64 rng(rd());//(world_rank);
    //if(world_rank == 0)

    b_statistical_ler(decoder, shots_per_k, rng, true, min_fault, distance);


}

void test_b_statistical_ler_AstreaG(uint64_t shots_per_k, uint distance, fp_t physcial_error, fp_t meas_er){ 

    
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
    const uint min_fault = 6;

    qrc::AstreaParams astreaG_param= {};
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    astreaG_param.use_mld = false;

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
   

    qrc::MWPMDecoder* decoder =  new  qrc::Astrea(circ,
                                        n_detectors_per_round,
                                        weight_filter_cutoff,
                                        astreaG_param);
    std::random_device rd;
    std::mt19937_64 rng(rd());
    
    b_statistical_ler(decoder, shots_per_k, rng, true, min_fault, distance);


}

void bb_ler_calculation(std::string decoder_name, std::ofstream& out, uint64_t shots_per_k, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, uint64_t SHOTS_PER_BATCH){
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_rank);

    std::vector<uint64_t> n_decoder_logcial_error_acc(max_k-min_k, 0);
    std::vector<uint64_t> n_decoder_logcial_error(max_k-min_k, 0);

    fp_t mwpm_ler = 0;

    bool found_le = false;
    uint64_t max_effective_k = max_k;

    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));

    if(world_rank == 0){
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Shots per bucket = " << shots_per_k << std::endl;
    }

    for(uint k = 0; k < max_k-min_k; k++){
        if(simulator.prob_k[k+min_k] < expected_ler){
            max_effective_k = k - 1;
            break;
        }

        uint64_t shots = shots_per_k / world_size;
        if (world_rank == world_size - 1) {
            shots += shots_per_k % world_size;
        }
        if(world_rank == 0 ){
                std::cout << "Started bucket " << k+min_k<< " - probability  = " << simulator.prob_k[k+min_k];
                if(0 < k)
                    std::cout << " Previous bucket faults for process "<<world_rank<< " = " << n_decoder_logcial_error[k-1];
                std::cout<<std::endl;
                
                out << "Started bucket " << k+min_k << " - probability  = " << simulator.prob_k[k+min_k];
                if(0 < k)
                    out << " Previous bucket faults for process "<<world_rank<< " = " <<  n_decoder_logcial_error[k-1];
                out<<std::endl;

        }
        std::vector<std::vector<uint8_t>> buffer;
        qrc::DecoderShotResult decoder_results;

        while(shots > 0){
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            buffer.clear();
            buffer = simulator.create_syndromes( k+min_k, shots_this_round, rng, false);

            for(auto& syndrome : buffer){
                decoder_results = decoder->decode_error(syndrome);
            
                if(decoder_results.is_logical_error){
                    n_decoder_logcial_error[k]++;
                }

            }
            shots -= shots_this_round;
        }

    }

    MPI_Reduce(&n_decoder_logcial_error[0], &n_decoder_logcial_error_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    mwpm_ler = 0;
    if(world_rank == 0){
        for(uint k = 0; k < max_k-min_k; k++){
            mwpm_ler +=  simulator.prob_k[k+min_k]*((fp_t)n_decoder_logcial_error_acc[k]/shots_per_k);
        }
    

        out << "Distance: "<< distance << " - #Buckets Max: " << max_k << " Min: "<< min_k <<  " - Shots Per Bucket: " << shots_per_k <<" - Physical ER: " << physcial_error << std::endl;
        out << decoder_name << " Decoder LER: " << mwpm_ler << std::endl;
        out << "________________________________________________________________" <<std::endl;
        out << " k (p) | # of Err in MWPM - impact on LER " << std::endl;

        for ( uint i = 0; i< max_effective_k+1;  i++){
            if(n_decoder_logcial_error_acc[i]) {
                out << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "<< n_decoder_logcial_error_acc[i] << " - " << (simulator.prob_k[i+min_k]*((fp_t)n_decoder_logcial_error_acc[i]/shots_per_k))<< std::endl;
            }
        }

        std::cout << "Distance: "<< distance << " - #Buckets Max: " << max_k << " Min: "<< min_k <<  " - Shots Per Bucket: " << shots_per_k <<" - Physical ER: " << physcial_error << std::endl;
        std::cout << decoder_name << " Decoder LER: " << mwpm_ler << std::endl;

    }

    delete decoder;     

}

void printing_samples_of_syndromes_and_matchings(uint64_t shots_per_k, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool print_ec_length){
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
    std::mt19937_64 rng;
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::MWPMDecoder decoder(circ);
    qpd::Predecoder predecoder(circ);

    auto path_table = qrc::compute_path_table(simulator.decoding_graph);
    bool had_long_ec = false;
    uint64_t total_sampling = 0;
    std::vector<qpd::SingleMatchingInfo> matchings;
    for(uint k = 0; k < max_k-min_k; k++){
        std::vector<std::vector<uint8_t>> buffer;
        qrc::DecoderShotResult decoder_results;
        std::map<uint, uint> printedPairs;
        total_sampling = 0;
        uint64_t s = 0;
        qpd::SingleMatchingInfo m_info = {};
        std::cout << "************** Number of Errors: " <<k+min_k <<" *****"<<std::endl;
        while(s < shots_per_k){
            total_sampling++;
            buffer = simulator.create_syndromes( k+min_k, 1, rng, false);
            printedPairs.clear();
            matchings.clear();
            //for(auto& syndrome : buffer){
                decoder_results = decoder.decode_error(buffer[0]);
                
                for(const auto& m : decoder_results.matching){
                    m_info = {};
                    uint first = m.first;
                    uint second = m.second;
                    if(printedPairs.count(first) == 0 && printedPairs.count(second) == 0){
                        printedPairs[first] = second;
                        m_info.first = first;
                        m_info.second = second;

                        auto v = simulator.decoding_graph.get_vertex(m.first);
                        auto w = simulator.decoding_graph.get_vertex(m.second);
                        auto v_w = std::make_pair(v, w);
                        // Finding error chains from the path table of decoding graph 
                        uint ec = path_table[v_w].path.size() - 1;
                        m_info.length = ec;
                        if(1 < ec){
                            had_long_ec = true;
                        }
                        m_info.chain_probability = pow(10,-1*path_table[v_w].distance);   
                        matchings.push_back(m_info);
                    }
                    
                }
                if(had_long_ec){
                    had_long_ec = false;
                    s++;
                }
                else{
                    continue;
                }
                std::cout << k+min_k << " errors - ";
                if(!print_ec_length){
                    qpd::print_syndrome(buffer[0]);
                }
                //std::cout << "HW = " << std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0) << std::endl;

                if(decoder_results.is_logical_error){
                    std::cout << "LOGICAL ERROR" << std::endl;
                }
                if(print_ec_length){
                    std::cout <<std::endl;
                    predecoder.update_adjacency_matrix_and_detector_list(buffer[0], 1, true);
                }
                std::cout << "___________MATCHING__________" << std::endl;
                for(auto m  : matchings){
                    std::cout << "(" << m.first << ", " << m.second << "), ";
                    if(print_ec_length){
                        std::cout << "(" << m.length << ": " << m.chain_probability <<") ";
                        predecoder.is_matched_parity_isolated(m);
                        if(!m.isolated){
                            std::cout <<"*- ";
                        }
                        else{
                            std::cout << "- ";
                        }
                    }
                }
                
                std::cout << std::endl<< "_____________________________" << std::endl;

            //}
        }
        if(print_ec_length){
            std::cout << "Total samples = \"" << total_sampling << "\" to have "<< s << " long ec." << std::endl;
        }
        
    }

}

void
error_chains_distribution(uint64_t shots_per_batch, bool use_mpi, uint n_faults, uint distance,  fp_t physcial_error, fp_t meas_er
                        , uint64_t min_k, uint64_t max_k) {
    int world_rank = 0, world_size = 1;
    if (use_mpi) {
        MPI_Comm_size(MPI_COMM_WORLD, &world_size);
        MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    }
    // setting circuit and mwpm decoder
    const stim::Circuit circuit = qrc::build_circuit(
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
    const uint min_fault = floor(distance/2);
    qrc::MWPMDecoder decoder(circuit);

    if(world_rank == 0)
        std::cout << "Distance = " << distance << " - Shots per k = " << shots_per_batch << " - Decoder: " << decoder.name() << std::endl;
    
    
    std::mt19937_64 rng(world_rank);
    //
    ///// For gathering data about EC
    std::vector<uint64_t> ec_distribution(100);
    std::fill(ec_distribution.begin(), ec_distribution.end(), 0);

    std::vector<std::vector<uint64_t>> syndromes_ecs;
    ////

    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circuit);
    auto path_table = qrc::compute_path_table(decoding_graph);
    

    const uint n_detectors = circuit.count_detectors();
    const uint n_observables = circuit.count_observables();
    const uint n_results = n_detectors + n_observables;

    uint64_t local_shots = shots_per_batch / world_size;
    if (world_rank == 0) {
        local_shots += shots_per_batch % world_size;
    }
    stim::simd_bit_table result_table(local_shots, n_results);

#define DELTA(a, b) ((a)-(b))/(a)
#define B_STAT_LER_USE_POLY

#ifdef B_STAT_LER_USE_POLY
    std::array<fp_t, 100'000> log_poly;
    log_poly.fill(1);

    uint eventno = 0;
#else
    fp_t mean_flip_prob = 0.0;
#endif
    std::vector<fp_t> edge_probs;
    std::vector<qrc::DecodingGraph::Edge*> edge_list;
    for (uint d = 0; d < n_detectors; d++) {
        auto v = decoding_graph.get_vertex(d);
        for (auto w : decoding_graph.adjacency_list(v)) {
            if(w->detector < d){
                continue;
            }
            auto edge = decoding_graph.get_edge(v, w);
            edge_list.push_back(edge);
            edge_probs.push_back(edge->error_probability);
#ifdef B_STAT_LER_USE_POLY
            fp_t ep = edge->error_probability;
            if (eventno == 0) {
                log_poly[0] = log(1-ep);
                log_poly[1] = log(ep);
            } else {
                std::array<fp_t, 100'000> log_pa, log_pb;
                log_pa.fill(1);
                log_pb.fill(1);
                for (uint i = 0; i <= eventno; i++) {
                    log_pa[i] = log_poly[i] + log(1-ep);
                    log_pb[i+1] = log_poly[i] + log(ep);
                }
                for (uint i = 0; i < log_poly.size(); i++) {
                    if (log_pa[i] == 1 && log_pb[i] == 1) {
                        log_poly[i] = 1;
                    } else if (log_pa[i] == 1) {
                        log_poly[i] = log_pb[i];
                    } else if (log_pb[i] == 1) {
                        log_poly[i] = log_pa[i];
                    } else {
                        log_poly[i] = log(pow(M_E, log_pa[i]) + pow(M_E, log_pb[i]));
                    }
                }
            }
            eventno++;
#else
            mean_flip_prob += 1.0/edge->error_probability;
#endif
        }
    }


#ifndef B_STAT_LER_USE_POLY
    mean_flip_prob = edge_list.size() / mean_flip_prob;
#endif
    std::discrete_distribution<> edge_dist(edge_probs.begin(), edge_probs.end());

    for (uint64_t s = 0; s < local_shots; s++) {
        for (uint i = 0; i < n_faults-1; i++) {
            // Add a random error.
            auto edge = edge_list[edge_dist(rng)];
            uint d1 = edge->detectors.first;
            uint d2 = edge->detectors.second;
            if (d1 != BOUNDARY_INDEX) {
                result_table[s][d1] ^= 1;
            }
            if (d2 != BOUNDARY_INDEX) {
                result_table[s][d2] ^= 1;
            }
            for (uint obs : edge->frames) {
                if (obs >= 0) {
                    result_table[s][n_detectors+obs] ^= 1;
                }
            }
        }
    }

    fp_t prev_logical_error_rate = 0.0;
    
    
    while (n_faults<max_k) {
        uint64_t local_errors = 0;

        for (uint64_t s = 0; s < local_shots; s++) {
            // Add a random error.
            auto edge = edge_list[edge_dist(rng)];
            uint d1 = edge->detectors.first;
            uint d2 = edge->detectors.second;
            if (d1 != BOUNDARY_INDEX) {
                result_table[s][d1] ^= 1;
            }
            if (d2 != BOUNDARY_INDEX) {
                result_table[s][d2] ^= 1;
            }
            for (uint obs : edge->frames) {
                if (obs >= 0) {
                    result_table[s][n_detectors+obs] ^= 1;
                }
            }

            auto syndrome = qrc::_to_vector(result_table[s], n_detectors, n_observables);
            auto res = decoder.decode_error(syndrome);
            for (auto pair : res.matching) {
                auto v = decoding_graph.get_vertex(pair.first);
                auto w = decoding_graph.get_vertex(pair.second);
                auto v_w = std::make_pair(v, w);
                // Finding error chains from the path table of decoding graph 
                uint ec = path_table[v_w].path.size() - 1;
                // count the number of error chains. Each cell of array i in ec_all_hw
                // contains the number of error chains with length i. 
                ec_distribution[ec]++;
            }
            syndromes_ecs.push_back(ec_distribution);
            std::fill(ec_distribution.begin(), ec_distribution.end(), 0);
            
        }

        n_faults++;
    }
    if (world_rank == 0) {
    

    }

    return;


}

void test_groups( uint distance, fp_t physcial_error, fp_t meas_er){
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

    std::vector<qrc::DecodingGraph::Vertex*> vertices = predecoder.decoding_graph.vertices();
    std::vector<qrc::DecodingGraph::Vertex*> a,b;
    for(uint i = 0; i<10; i++){
        a.push_back(vertices[i]);
        std::cout << vertices[i]->detector << " ";
    }
    std::cout << std::endl;
    for(uint i = 6; i<20; i++){
        b.push_back(vertices[i]);
        std::cout << vertices[i]->detector << " ";
    }
    std::cout << std::endl;
    qpd::CountResult c = qpd::countMembersInAB(a,b);
    std::cout << c.countANotInB << "  " << c.countBNotInA << " "<< c.countBothAB;


}
void test_finding_fast_group(uint64_t shots_per_k, uint distance, 
                            fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, std::string decoder_name){
    /*
    In this fucntion, we want to test our method that group vertices for fast matching.
    In this function, by A we mean group of vertices that are matched to their adjacent vertex by MWPM.
    And by               B we mean group of vertices that our method predict that will be matched with
                         their adjacent vertices (Fast group).
    */
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }

    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);

    auto path_table = qrc::compute_path_table(predecoder->decoding_graph);

    std::vector<fp_t> sum_IR(max_k-min_k, 0);
    std::vector<fp_t> sum_AR(max_k-min_k, 0);
    std::vector<uint64_t> times_IR_is_1(max_k-min_k, 0);
    std::vector<uint64_t> times_AR_is_1(max_k-min_k, 0);

    std::vector<fp_t> sum_IR_acc(max_k-min_k, 0);
    std::vector<fp_t> sum_AR_acc(max_k-min_k, 0);
    std::vector<uint64_t> times_IR_is_1_acc(max_k-min_k, 0);
    std::vector<uint64_t> times_AR_is_1_acc(max_k-min_k, 0);

    std::vector<uint64_t> sum_A_size(max_k-min_k, 0);
    std::vector<uint64_t> sum_B_size(max_k-min_k, 0);

    std::vector<uint64_t> sum_A_size_acc(max_k-min_k, 0);
    std::vector<uint64_t> sum_B_size_acc(max_k-min_k, 0);

    std::vector<uint> sum_remained_HW_A(max_k-min_k, 0);
    std::vector<uint> sum_remained_HW_B(max_k-min_k, 0);

    std::vector<uint> sum_remained_HW_A_acc(max_k-min_k, 0);
    std::vector<uint> sum_remained_HW_B_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_4(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_4_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_5(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_5_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_6(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_6_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_7(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_7_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_8(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_8_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_9(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_9_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_minusAB(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_minusAB_acc(max_k-min_k, 0);

    std::vector<uint64_t> sum_countBNotInA(max_k-min_k, 0);
    std::vector<uint64_t> sum_countBNotInA_acc(max_k-min_k, 0);

    std::vector<uint64_t> times_more_than_1_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_1_incorrect_pair_acc(max_k-min_k, 0);

    std::vector<uint64_t> times_more_than_2_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_2_incorrect_pair_acc(max_k-min_k, 0);

    std::vector<uint64_t> times_more_than_3_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_3_incorrect_pair_acc(max_k-min_k, 0);

    std::vector<uint64_t> times_more_than_4_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_4_incorrect_pair_acc(max_k-min_k, 0);
    
    std::vector<uint64_t> times_more_than_5_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_5_incorrect_pair_acc(max_k-min_k, 0);


    // By bi, we mean boundary node involved
    std::vector<uint64_t> ec_2(max_k-min_k, 0);
    std::vector<uint64_t> ec_2_acc(max_k-min_k, 0);

    std::vector<uint64_t> ec_3(max_k-min_k, 0);
    std::vector<uint64_t> ec_3_acc(max_k-min_k, 0);

    std::vector<uint64_t> ec_4(max_k-min_k, 0);
    std::vector<uint64_t> ec_4_acc(max_k-min_k, 0);

    std::vector<uint64_t> ec_5(max_k-min_k, 0);
    std::vector<uint64_t> ec_5_acc(max_k-min_k, 0);

    std::vector<uint64_t> ec_6(max_k-min_k, 0);
    std::vector<uint64_t> ec_6_acc(max_k-min_k, 0);

    std::vector<uint64_t> ec_7(max_k-min_k, 0);
    std::vector<uint64_t> ec_7_acc(max_k-min_k, 0);

    // By bi, we mean boundary node involved
    std::vector<uint64_t> bi_ec_2(max_k-min_k, 0);
    std::vector<uint64_t> bi_ec_2_acc(max_k-min_k, 0);

    std::vector<uint64_t> bi_ec_3(max_k-min_k, 0);
    std::vector<uint64_t> bi_ec_3_acc(max_k-min_k, 0);

    std::vector<uint64_t> bi_ec_4(max_k-min_k, 0);
    std::vector<uint64_t> bi_ec_4_acc(max_k-min_k, 0);

    std::vector<uint64_t> bi_ec_5(max_k-min_k, 0);
    std::vector<uint64_t> bi_ec_5_acc(max_k-min_k, 0);

    std::vector<uint64_t> bi_ec_6(max_k-min_k, 0);
    std::vector<uint64_t> bi_ec_6_acc(max_k-min_k, 0);

    std::vector<uint64_t> bi_ec_7(max_k-min_k, 0);
    std::vector<uint64_t> bi_ec_7_acc(max_k-min_k, 0);
    
    



    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));

    if(world_rank == 0){
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Shots per bucket = " << shots_per_k << std::endl;
    }

    for(uint k = 0; k < max_k-min_k; k++){

        uint64_t shots = shots_per_k / world_size;
        if (world_rank == world_size - 1) {
            shots += shots_per_k % world_size;
        }
        if(world_rank == 0 ){
                std::cout << "Started bucket " << k+min_k<< " - probability  = " << simulator.prob_k[k+min_k] << std::endl;

        }
        std::vector<std::vector<uint8_t>> buffer;
        std::vector<qrc::DecodingGraph::Vertex*> A; // Actual Fast Group
        std::vector<qrc::DecodingGraph::Vertex*> B; // Predicted Fast Group
        std::vector<qrc::DecodingGraph::Vertex*> C; // Boundary connected group
        qrc::DecoderShotResult decoder_results;
        uint hw = 0;
        bool not_found_ec_2, not_found_ec_3, not_found_ec_4, not_found_ec_5, not_found_ec_6, not_found_ec_7;
        bool not_found_bi_ec_2, not_found_bi_ec_3, not_found_bi_ec_4, not_found_bi_ec_5, not_found_bi_ec_6, not_found_bi_ec_7;
        bool found_boundary_in_chain;
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            buffer.clear();
            buffer = simulator.create_syndromes( k+min_k, shots_this_round, rng, false);

            for(auto& syndrome : buffer){
                A.clear();
                B.clear();
                C.clear();

                // Need these booleans for only counting wheather a ec with
                // a specific length exist or not, rather than couting all of them
                not_found_bi_ec_2 = true; not_found_bi_ec_3 = true; not_found_bi_ec_4 = true; not_found_bi_ec_5 = true;
                not_found_bi_ec_6 = true; not_found_bi_ec_7 = true;

                not_found_ec_2 = true; not_found_ec_3 = true; not_found_ec_4 = true; not_found_ec_5 = true;
                not_found_ec_6 = true; not_found_ec_7 = true;

                found_boundary_in_chain = false;

                hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);


                //Creating B group
                B = predecoder->get_fast_matching_group(syndrome);

                decoder_results = decoder->decode_error(syndrome);

                //Iterating over matching result to find how many of them has ec of length 1
                
                for(const auto& matching: decoder_results.matching){
                    qrc::DecodingGraph::Vertex* v1 = predecoder->decoding_graph.get_vertex(matching.first);
                    qrc::DecodingGraph::Vertex* v2 = predecoder->decoding_graph.get_vertex(matching.second);

                    if(predecoder->decoding_graph.get_edge(v1, v2) != nullptr){
                        //Creating A group
                        qpd::addUniqueVertex(A, v1);
                        qpd::addUniqueVertex(A, v2);
                    }
                    else{
                        auto v12 = std::make_pair(v1, v2);
                        // Finding error chains from the path table of decoding graph 
                        ///Path table and the path!
                        for(const auto& p :path_table[v12].path){
                            if(p->detector == BOUNDARY_INDEX){
                                qpd::addUniqueVertex(C, v1);
                                qpd::addUniqueVertex(C, v2);
                                found_boundary_in_chain = true;

                                if(not_found_bi_ec_2 && (path_table[v12].path.size() == 3)){
                                    bi_ec_2[k]++;
                                    not_found_bi_ec_2 = false;
                                }
                                else if(not_found_bi_ec_3 && (path_table[v12].path.size() == 4)){
                                    bi_ec_3[k]++;
                                    not_found_bi_ec_3 = false;
                                }
                                else if(not_found_bi_ec_4 && (path_table[v12].path.size() == 5)){
                                    bi_ec_4[k]++;
                                    not_found_bi_ec_4 = false;
                                }
                                else if(not_found_bi_ec_5 && (path_table[v12].path.size() == 6)){
                                    bi_ec_5[k]++;
                                    not_found_bi_ec_5 = false;
                                }
                                else if(not_found_bi_ec_6 && (path_table[v12].path.size() == 7)){
                                    bi_ec_6[k]++;
                                    not_found_bi_ec_6 = false;
                                }
                                else if(not_found_bi_ec_7 && (path_table[v12].path.size() == 8)){
                                    bi_ec_7[k]++;
                                    not_found_bi_ec_7 = false;
                                }
                            }
                        }
                        if(!found_boundary_in_chain){
                            
                            if(not_found_ec_2 && (path_table[v12].path.size() == 3)){
                                ec_2[k]++;
                                not_found_ec_2 = false;
                            }
                            else if(not_found_ec_3 && (path_table[v12].path.size() == 4)){
                                ec_3[k]++;
                                not_found_ec_3 = false;
                            }
                            else if(not_found_ec_4 && (path_table[v12].path.size() == 5)){
                                ec_4[k]++;
                                not_found_ec_4 = false;
                            }
                            else if(not_found_ec_5 && (path_table[v12].path.size() == 6)){
                                ec_5[k]++;
                                not_found_ec_5 = false;
                            }
                            else if(not_found_ec_6 && (path_table[v12].path.size() == 7)){
                                ec_6[k]++;
                                not_found_ec_6 = false;
                            }
                            else if(not_found_ec_7 && (path_table[v12].path.size() == 8)){
                                ec_7[k]++;
                                not_found_ec_7 = false;
                            }
                        
                        }
                    }
                }

                qpd::CountResult c = qpd::countMembersInAB(A,B);

                fp_t IR = ((double)c.countBothAB/(c.countBothAB+c.countANotInB));
                if((c.countBothAB + c.countANotInB) != A.size()){
                    std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. A - B + AandB != A";
                }
                sum_IR[k] += IR;
                if(IR == 1){
                    times_IR_is_1[k]++;
                }

                fp_t AR = B.size()!=0 ?((double)c.countBothAB/(c.countBothAB+c.countBNotInA)):0;

                if((c.countBothAB + c.countBNotInA) != B.size()){
                    std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. B - A + AandB != B";
                }
                sum_AR[k] += AR;
                if(AR == 1){
                    times_AR_is_1[k]++;
                }
                sum_A_size[k] += A.size();
                sum_B_size[k] += B.size();

                if( c.countBothAB+c.countANotInB < hw){
                    sum_remained_HW_A[k] += hw-(c.countBothAB+c.countANotInB);
                    if( 10<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar[k]++;
                    }
                    if( 9<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar_9[k]++;
                    }
                    if( 8<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar_8[k]++;
                    }
                    if( 7<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar_7[k]++;
                    }
                    if( 6<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar_6[k]++;
                    }
                    if( 5<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar_5[k]++;
                    }
                    if( 4<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar_4[k]++;
                    }
                }

                if( c.countBothAB+c.countBNotInA < hw){
                    sum_remained_HW_B[k] += hw-(c.countBothAB+c.countBNotInA);
                    
                }
                if( 10<(hw-c.countBothAB)){
                        above_bar_minusAB[k]++;
                    }

                sum_countBNotInA[k] += c.countBNotInA;

                if(2 < c.countBNotInA){
                    times_more_than_1_incorrect_pair[k]++;
                }
                if(4 < c.countBNotInA){
                    times_more_than_2_incorrect_pair[k]++;
                }
                if(6 < c.countBNotInA){
                    times_more_than_3_incorrect_pair[k]++;
                }
                if(8 < c.countBNotInA){
                    times_more_than_4_incorrect_pair[k]++;
                }
                if(10 < c.countBNotInA){
                    times_more_than_5_incorrect_pair[k]++;
                }

            }
            shots -= shots_this_round;
        }

    }
    
    MPI_Reduce(&sum_IR[0], &sum_IR_acc[0], max_k-min_k, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_AR[0], &sum_AR_acc[0], max_k-min_k, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&times_AR_is_1[0], &times_AR_is_1_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_IR_is_1[0], &times_IR_is_1_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
     
    MPI_Reduce(&sum_A_size[0], &sum_A_size_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_B_size[0], &sum_B_size_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&sum_remained_HW_A[0], &sum_remained_HW_A_acc[0], max_k-min_k, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_remained_HW_B[0], &sum_remained_HW_B_acc[0], max_k-min_k, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&above_bar[0], &above_bar_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&above_bar_4[0], &above_bar_4_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&above_bar_5[0], &above_bar_5_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&above_bar_6[0], &above_bar_6_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&above_bar_7[0], &above_bar_7_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&above_bar_8[0], &above_bar_8_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&above_bar_9[0], &above_bar_9_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&above_bar_minusAB[0], &above_bar_minusAB_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&sum_countBNotInA[0], &sum_countBNotInA_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&times_more_than_1_incorrect_pair[0], &times_more_than_1_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_more_than_2_incorrect_pair[0], &times_more_than_2_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_more_than_3_incorrect_pair[0], &times_more_than_3_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_more_than_4_incorrect_pair[0], &times_more_than_4_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_more_than_5_incorrect_pair[0], &times_more_than_5_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&bi_ec_2[0], &bi_ec_2_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&bi_ec_3[0], &bi_ec_3_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&bi_ec_4[0], &bi_ec_4_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&bi_ec_5[0], &bi_ec_5_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&bi_ec_6[0], &bi_ec_6_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&bi_ec_7[0], &bi_ec_7_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&ec_2[0], &ec_2_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ec_3[0], &ec_3_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ec_4[0], &ec_4_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ec_5[0], &ec_5_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ec_6[0], &ec_6_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ec_7[0], &ec_7_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);


    if(world_rank == 0){ 
        std::string printing_shots = std::to_string(shots_per_k);
        if(1'000<=shots_per_k){
            printing_shots = number_of_shots_to_string(shots_per_k);
        }

        std::cout << "Distance: "<< distance << " - #Buckets Max: " << max_k << " Min: "<< min_k <<  " - Shots Per Bucket: " << shots_per_k <<" - Physical ER: " << physcial_error << std::endl;
        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | ";
        std::cout << std::setw(20) << " P(IR = 1) | " ;
        std::cout << std::setw(20) << " P(AR = 1) | ";
        std::cout << std::setw(1) << " E[IR] | ";
        std::cout << std::setw(16) << " E[AR] | ";
        std::cout << std::setw(1) << " E[|A|], E[HW-|A|] | "; 
        std::cout << std::setw(1) << " E[|B|], E[HW-|B|] |" << std::endl;

        for ( uint i = 0; i< max_k-min_k;  i++){
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(10) <<((double)times_IR_is_1_acc[i]/shots_per_k) << "(" << times_IR_is_1_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)times_AR_is_1_acc[i]/shots_per_k) << "(" << times_AR_is_1_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(14) <<sum_IR_acc[i]/shots_per_k << " | ";
            std::cout << std::setw(14) <<sum_AR_acc[i]/shots_per_k << " | ";
            std::cout << std::setw(5) <<sum_A_size_acc[i]/shots_per_k << ", " << sum_remained_HW_A_acc[i]/shots_per_k << " | ";
            std::cout << std::setw(5) <<sum_B_size_acc[i]/shots_per_k << ", " << sum_remained_HW_B_acc[i]/shots_per_k << " | " << std::endl;
        
        }
        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | "; 
        std::cout <<  std::setw(16) << " P(4<|HW-|A|) | "; 
        std::cout <<  std::setw(16) << " P(5<|HW-|A|) | "; 
        std::cout <<  std::setw(16) << " P(6<|HW-|A|) | "; 
        std::cout <<  std::setw(16) << " P(7<|HW-|A|) | "; 
        std::cout <<  std::setw(16) << " P(8<|HW-|A|) | "; 
        std::cout <<  std::setw(16) << " P(9<|HW-|A|) | "; 
        std::cout <<  std::setw(16) << " P(10<|HW-|A|) | "; 
        std::cout << std::setw(14) << " P(10<|HW-|AB||) | "<< std::endl;
        for ( uint i = 0; i< max_k-min_k;  i++){
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(10) <<((double)above_bar_4_acc[i]/shots_per_k) << "(" << above_bar_4_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_5_acc[i]/shots_per_k) << "(" << above_bar_5_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_6_acc[i]/shots_per_k) << "(" << above_bar_6_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_7_acc[i]/shots_per_k) << "(" << above_bar_7_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_8_acc[i]/shots_per_k) << "(" << above_bar_8_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_9_acc[i]/shots_per_k) << "(" << above_bar_9_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_acc[i]/shots_per_k) << "(" << above_bar_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_minusAB_acc[i]/shots_per_k) << "(" << above_bar_minusAB_acc[i] << "/" << printing_shots <<") | "
            << std::endl;
            
        }

        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | "; 
        std::cout << std::setw(6) << " E[|B-A|] | "; 
        std::cout << std::setw(6) << " P(2<|B-A|) |";
        std::cout << std::setw(6) << " P(4<|B-A|) |";
        std::cout << std::setw(6) << " P(6<|B-A|) |";
        std::cout << std::setw(6) << " P(8<|B-A|) |";
        std::cout << std::setw(6) << " P(10<|B-A|) |" << std::endl;

        for ( uint i = 0; i< max_k-min_k;  i++){
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(6) <<((double)sum_countBNotInA_acc[i]/shots_per_k) << " | ";
            std::cout << std::setw(6) <<((double)times_more_than_1_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_1_incorrect_pair_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)times_more_than_2_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_2_incorrect_pair_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)times_more_than_3_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_3_incorrect_pair_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)times_more_than_4_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_4_incorrect_pair_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)times_more_than_5_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_5_incorrect_pair_acc[i] << "/" << printing_shots <<") | "

            << std::endl;
            
        }

        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | "; 
        std::cout << std::setw(6) << " P(|ec|=2 & bi) | "; 
        std::cout << std::setw(6) << " P(|ec|=3 & bi) |";
        std::cout << std::setw(6) << " P(|ec|=4 & bi) |";
        std::cout << std::setw(6) << " P(|ec|=5 & bi) |";
        std::cout << std::setw(6) << " P(|ec|=6 & bi |";
        std::cout << std::setw(6) << " P(|ec|=7 & bi) |" << std::endl;
        // The above ones means probability of existance of a error chain of length x and 
        // The error chain has a boundary in one matching set.

        for ( uint i = 0; i< max_k-min_k;  i++){
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(6) <<((double)bi_ec_2_acc[i]/shots_per_k) << "(" << bi_ec_2_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)bi_ec_3_acc[i]/shots_per_k) << "(" << bi_ec_3_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)bi_ec_4_acc[i]/shots_per_k) << "(" << bi_ec_4_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)bi_ec_5_acc[i]/shots_per_k) << "(" << bi_ec_5_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)bi_ec_6_acc[i]/shots_per_k) << "(" << bi_ec_6_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)bi_ec_7_acc[i]/shots_per_k) << "(" << bi_ec_7_acc[i] << "/" << printing_shots <<") | "

            << std::endl;
            
        }
        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | "; 
        std::cout << std::setw(6) << " P(|ec|=2 & !bi) | "; 
        std::cout << std::setw(6) << " P(|ec|=3 & !bi) |";
        std::cout << std::setw(6) << " P(|ec|=4 & !bi) |";
        std::cout << std::setw(6) << " P(|ec|=5 & !bi) |";
        std::cout << std::setw(6) << " P(|ec|=6 & !bi) |";
        std::cout << std::setw(6) << " P(|ec|=7 & !bi) |" << std::endl;
        // The above ones means probability of existance of a error chain of length x and 
        // The error chain does not have a boundary in one matching set.

        for ( uint i = 0; i< max_k-min_k;  i++){
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(6) <<((double)ec_2_acc[i]/shots_per_k) << "(" << ec_2_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)ec_3_acc[i]/shots_per_k) << "(" << ec_3_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)ec_4_acc[i]/shots_per_k) << "(" << ec_4_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)ec_5_acc[i]/shots_per_k) << "(" << ec_5_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)ec_6_acc[i]/shots_per_k) << "(" << ec_6_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)ec_7_acc[i]/shots_per_k) << "(" << ec_7_acc[i] << "/" << printing_shots <<") | "

            << std::endl;
            
        }

    }
    
}

void printing_incorrect_predictions(uint64_t shots_per_k, uint distance, 
                            fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k
                            , std::string decoder_name){
    /*
    In this fucntion, we want to test our method that group vertices for fast matching.
    In this function, by A we mean group of vertices that are matched to their adjacent vertex by MWPM.
    And by               B we mean group of vertices that our method predict that will be matched with
                         their adjacent vertices (Fast group).
    */
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);

    std::vector<fp_t> sum_IR(max_k-min_k, 0);
    std::vector<fp_t> sum_AR(max_k-min_k, 0);
    std::vector<uint64_t> times_IR_is_1(max_k-min_k, 0);
    std::vector<uint64_t> times_AR_is_1(max_k-min_k, 0);

    std::vector<fp_t> sum_IR_acc(max_k-min_k, 0);
    std::vector<fp_t> sum_AR_acc(max_k-min_k, 0);
    std::vector<uint64_t> times_IR_is_1_acc(max_k-min_k, 0);
    std::vector<uint64_t> times_AR_is_1_acc(max_k-min_k, 0);

    std::vector<uint64_t> sum_A_size(max_k-min_k, 0);
    std::vector<uint64_t> sum_B_size(max_k-min_k, 0);

    std::vector<uint64_t> sum_A_size_acc(max_k-min_k, 0);
    std::vector<uint64_t> sum_B_size_acc(max_k-min_k, 0);

    std::vector<uint> sum_remained_HW_A(max_k-min_k, 0);
    std::vector<uint> sum_remained_HW_B(max_k-min_k, 0);

    std::vector<uint> sum_remained_HW_A_acc(max_k-min_k, 0);
    std::vector<uint> sum_remained_HW_B_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_acc(max_k-min_k, 0);

    std::vector<uint64_t> above_bar_minusAB(max_k-min_k, 0);
    std::vector<uint64_t> above_bar_minusAB_acc(max_k-min_k, 0);

    std::vector<uint64_t> sum_countBNotInA(max_k-min_k, 0);
    std::vector<uint64_t> sum_countBNotInA_acc(max_k-min_k, 0);

    std::vector<uint64_t> times_more_than_1_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_1_incorrect_pair_acc(max_k-min_k, 0);

    std::vector<uint64_t> times_more_than_2_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_2_incorrect_pair_acc(max_k-min_k, 0);

    std::vector<uint64_t> times_more_than_3_incorrect_pair(max_k-min_k, 0);
    std::vector<uint64_t> times_more_than_3_incorrect_pair_acc(max_k-min_k, 0);
    



    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));

    if(world_rank == 0){
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Shots per bucket = " << shots_per_k << std::endl;
    }

    for(uint k = 0; k < max_k-min_k; k++){

        uint64_t shots = shots_per_k / world_size;
        if (world_rank == world_size - 1) {
            shots += shots_per_k % world_size;
        }
        std::cout << "************** Number of Errors: " <<k+min_k <<" - P = " << simulator.prob_k[k+min_k] <<" *****"<<std::endl;

        std::vector<std::vector<uint8_t>> buffer;
        std::vector<qrc::DecodingGraph::Vertex*> A;
        std::vector<qrc::DecodingGraph::Vertex*> B;
        qrc::DecoderShotResult decoder_results;
        std::map<uint, uint> printedPairs;
        uint hw = 0;

        while(shots > 0){
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);

            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                A.clear();
                B.clear();

                hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);


                //Creating B group
                //B = predecoder->get_fast_matching_group(syndrome);
                B = predecoder->get_fast_matching_group_v2(syndrome);

                decoder_results = decoder->decode_error(syndrome);

                //Iterating over matching result to find how many of them has ec of length 1
                for(const auto& matching: decoder_results.matching){
                    qrc::DecodingGraph::Vertex* v1 = predecoder->decoding_graph.get_vertex(matching.first);
                    qrc::DecodingGraph::Vertex* v2 = predecoder->decoding_graph.get_vertex(matching.second);

                    if(predecoder->decoding_graph.get_edge(v1, v2) != nullptr){
                        //Creating A group
                        qpd::addUniqueVertex(A, v1);
                        qpd::addUniqueVertex(A, v2);
                    }
                    
                }

                qpd::CountResult c = qpd::countMembersInAB(A,B);

                fp_t IR = ((double)c.countBothAB/(c.countBothAB+c.countANotInB));
                if((c.countBothAB + c.countANotInB) != A.size()){
                    std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. A - B + AandB != A";
                }
                sum_IR[k] += IR;
                if(IR == 1){
                    times_IR_is_1[k]++;
                }

                fp_t AR = B.size()!=0 ?((double)c.countBothAB/(c.countBothAB+c.countBNotInA)):0;

                if((c.countBothAB + c.countBNotInA) != B.size()){
                    std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. B - A + AandB != B";
                }
                sum_AR[k] += AR;
                if(AR == 1){
                    times_AR_is_1[k]++;
                }
                sum_A_size[k] += A.size();
                sum_B_size[k] += B.size();

                if( c.countBothAB+c.countANotInB < hw){
                    sum_remained_HW_A[k] += hw-(c.countBothAB+c.countANotInB);
                    if( 10<(hw-(c.countBothAB+c.countANotInB))){
                        above_bar[k]++;
                    }
                }

                if( c.countBothAB+c.countBNotInA < hw){
                    sum_remained_HW_B[k] += hw-(c.countBothAB+c.countBNotInA);
                    
                }
                if( 10<(hw-c.countBothAB)){
                        above_bar_minusAB[k]++;
                    }

                sum_countBNotInA[k] += c.countBNotInA;

                if(2 < c.countBNotInA){
                    times_more_than_1_incorrect_pair[k]++;
                }
                if(4 < c.countBNotInA){
                    times_more_than_2_incorrect_pair[k]++;
                }
                if(6 < c.countBNotInA){
                    times_more_than_3_incorrect_pair[k]++;
                }


                // Checking if c.countBNotInA != 0 -> detecting false positive
                // Checking if c.countANotInB != 0 -> detecting false Negative

                
                if(c.countANotInB != 0){
                    printedPairs.clear();
                    std::cout << k+min_k << " errors - ";
                    qpd::print_syndrome(syndrome);
                    if(decoder_results.is_logical_error){
                        std::cout << "LOGICAL ERROR" << std::endl;
                    }
                    std::cout << "___________MATCHING__________" << std::endl;
                    for(const auto& m: decoder_results.matching){
                        uint first = m.first;
                        uint second = m.second;
                        if(printedPairs.count(first) == 0 && printedPairs.count(second) == 0){
                            printedPairs[first] = second;
                            std::cout << "(" << m.first << ", " << m.second << "), ";   
                        }   
                    }
                    std::cout << std::endl<< "_____________________________" << std::endl;

                    std::cout << std::endl;
                    std::cout << "Actual Fast Group (A): ";
                    for(const auto& a: A){
                        std::cout<< a->detector << ",";
                    }
                    std::cout << std::endl;
                    std::cout << "Predicted Fast Group (B): ";
                    for(const auto& b: B){
                        std::cout<< b->detector << ",";
                    }
                    std::cout << std::endl;
                    std::cout << "False Negative (A-B): ";
                    for(const auto& a: A){
                        auto itA = std::find_if(B.begin(), B.end(), [&](qrc::DecodingGraph::Vertex* vertexB) {
                            return vertexB->detector == a->detector;
                        });

                        if(itA == B.end())
                            std::cout<< a->detector << ",";
                    }
                    std::cout << std::endl;
                    std::cout << "False Positive (B-A): ";
                    for(const auto& b: B){
                        auto itB = std::find_if(A.begin(), A.end(), [&](qrc::DecodingGraph::Vertex* vertexA) {
                            return vertexA->detector == b->detector;
                        });

                        if(itB == A.end())
                            std::cout<< b->detector << ",";
                    }
                    std::cout << std::endl;
                    std::cout <<" |A| = " << c.countBothAB+c.countANotInB << 
                                " |B| = " << c.countBothAB+c.countBNotInA <<
                                " |AB| = " << c.countBothAB <<
                                " |A-B| = " << c.countANotInB <<
                                " |B-A| = " << c.countBNotInA << std::endl;
                    std::cout << std::endl<< "_____________________________" << std::endl;



                }

            }
            shots -= shots_this_round;
        }

    }
    
    MPI_Reduce(&sum_IR[0], &sum_IR_acc[0], max_k-min_k, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_AR[0], &sum_AR_acc[0], max_k-min_k, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(&times_AR_is_1[0], &times_AR_is_1_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_IR_is_1[0], &times_IR_is_1_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
     
    MPI_Reduce(&sum_A_size[0], &sum_A_size_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_B_size[0], &sum_B_size_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&sum_remained_HW_A[0], &sum_remained_HW_A_acc[0], max_k-min_k, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sum_remained_HW_B[0], &sum_remained_HW_B_acc[0], max_k-min_k, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&above_bar[0], &above_bar_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&above_bar_minusAB[0], &above_bar_minusAB_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&sum_countBNotInA[0], &sum_countBNotInA_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&times_more_than_1_incorrect_pair[0], &times_more_than_1_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_more_than_2_incorrect_pair[0], &times_more_than_2_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&times_more_than_3_incorrect_pair[0], &times_more_than_3_incorrect_pair_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    if(world_rank == 0){ 
        std::string printing_shots = std::to_string(shots_per_k);
        if(1'000<=shots_per_k){
            printing_shots = number_of_shots_to_string(shots_per_k);
        }

        std::cout << "Distance: "<< distance << " - #Buckets Max: " << max_k << " Min: "<< min_k <<  " - Shots Per Bucket: " << shots_per_k <<" - Physical ER: " << physcial_error << std::endl;
        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | ";
        std::cout << std::setw(20) << " P(IR = 1) | " ;
        std::cout << std::setw(20) << " P(AR = 1) | ";
        std::cout << std::setw(1) << " E[IR] | ";
        std::cout << std::setw(16) << " E[AR] | ";
        std::cout << std::setw(1) << " E[|A|], E[HW-|A|] | "; 
        std::cout << std::setw(1) << " E[|B|], E[HW-|B|] |" << std::endl;

        for ( uint i = 0; i< max_k-min_k;  i++){
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(10) <<((double)times_IR_is_1_acc[i]/shots_per_k) << "(" << times_IR_is_1_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)times_AR_is_1_acc[i]/shots_per_k) << "(" << times_AR_is_1_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(14) <<sum_IR_acc[i]/shots_per_k << " | ";
            std::cout << std::setw(14) <<sum_AR_acc[i]/shots_per_k << " | ";
            std::cout << std::setw(5) <<sum_A_size_acc[i]/shots_per_k << ", " << sum_remained_HW_A_acc[i]/shots_per_k << " | ";
            std::cout << std::setw(5) <<sum_B_size_acc[i]/shots_per_k << ", " << sum_remained_HW_B_acc[i]/shots_per_k << " | " << std::endl;
        
        }
        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | "; 
        std::cout <<  std::setw(16) << " P(10<|HW-|A|) | "; 
        std::cout << std::setw(14) << " P(10<|HW-|AB||) | ";
        std::cout << std::setw(6) << " E[|B-A|] | "; 
        std::cout << std::setw(6) << " P(2<|B-A|) |";
        std::cout << std::setw(6) << " P(4<|B-A|) |";
        std::cout << std::setw(6) << " P(6<|B-A|) |" << std::endl;

        for ( uint i = 0; i< max_k-min_k;  i++){
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(10) <<((double)above_bar_acc[i]/shots_per_k) << "(" << above_bar_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)above_bar_minusAB_acc[i]/shots_per_k) << "(" << above_bar_minusAB_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)sum_countBNotInA_acc[i]/shots_per_k) << " | ";
            std::cout << std::setw(6) <<((double)times_more_than_1_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_1_incorrect_pair_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)times_more_than_2_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_2_incorrect_pair_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(6) <<((double)times_more_than_3_incorrect_pair_acc[i]/shots_per_k) << "(" << times_more_than_3_incorrect_pair_acc[i] << "/" << printing_shots <<") | "

            << std::endl;
            
        }

    }
    
}

void test_priority_groups(uint64_t max_shot, uint distance, fp_t physcial_error, 
            fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name){

/*
    In this fucntion, we want to test our method that group vertices for fast matching.
    In this function, by A we mean group of vertices that are matched to their adjacent vertex by MWPM.
    And by               B we mean group of vertices that our method predict that will be matched with
                         their adjacent vertices (Fast group).
    */
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);
    

    auto path_table = qrc::compute_path_table(predecoder->decoding_graph);

    std::vector<std::vector<uint64_t>> times_AR_is_1(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_AR_is_1_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> times_CR_is_1(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_CR_is_1_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    std::vector<std::vector<uint64_t>> times_AR2_is_1(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_AR2_is_1_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> times_CR2_is_1(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_CR2_is_1_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    std::vector<std::vector<fp_t>> AR_avg(6, std::vector<fp_t>(max_k - min_k, 0));
    std::vector<std::vector<fp_t>> AR_avg_acc(6, std::vector<fp_t>(max_k - min_k, 0));

    std::vector<std::vector<fp_t>> CR_avg(6, std::vector<fp_t>(max_k - min_k, 0));
    std::vector<std::vector<fp_t>> CR_avg_acc(6, std::vector<fp_t>(max_k - min_k, 0));

    std::vector<std::vector<fp_t>> AR2_avg(6, std::vector<fp_t>(max_k - min_k, 0));
    std::vector<std::vector<fp_t>> AR2_avg_acc(6, std::vector<fp_t>(max_k - min_k, 0));

    std::vector<std::vector<fp_t>> CR2_avg(6, std::vector<fp_t>(max_k - min_k, 0));
    std::vector<std::vector<fp_t>> CR2_avg_acc(6, std::vector<fp_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> PG_high(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> PG_high_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> PG_10(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> PG_10_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    std::vector<std::vector<uint64_t>> PG_2(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> PG_2_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    std::vector<std::vector<uint64_t>> PG_0(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> PG_0_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> PG_4(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> PG_4_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> PG_6(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> PG_6_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> PG_8(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> PG_8_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    std::vector<std::vector<uint64_t>> times_1_incorrect_pair(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_1_incorrect_pair_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> times_2_incorrect_pair(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_2_incorrect_pair_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> times_3_incorrect_pair(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_3_incorrect_pair_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> times_4_incorrect_pair(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> times_4_incorrect_pair_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    std::vector<std::vector<uint64_t>> more_than_4_incorrect_pair(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> more_than_4_incorrect_pair_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    // for three incorrect pairs
    std::vector<std::vector<uint64_t>> inc_pair_3_hw_6(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_3_hw_6_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> inc_pair_3_hw_8(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_3_hw_8_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> inc_pair_3_hw_10(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_3_hw_10_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> inc_pair_3_hhw(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_3_hhw_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    // for two incorrect pairs
    std::vector<std::vector<uint64_t>> inc_pair_2_hw_8(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_2_hw_8_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> inc_pair_2_hw_10(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_2_hw_10_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> inc_pair_2_hhw(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_2_hhw_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    // for one incorrect pairs
    std::vector<std::vector<uint64_t>> inc_pair_1_hw_10(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_1_hw_10_acc(6, std::vector<uint64_t>(max_k - min_k, 0));

    std::vector<std::vector<uint64_t>> inc_pair_1_hhw(6, std::vector<uint64_t>(max_k - min_k, 0));
    std::vector<std::vector<uint64_t>> inc_pair_1_hhw_acc(6, std::vector<uint64_t>(max_k - min_k, 0));
    
    
    //For seeing remained HW
    std::vector<uint64_t> hw_r_high(max_k - min_k, 0);
    std::vector<uint64_t> hw_r_high_acc(max_k - min_k, 0);

    std::vector<uint64_t> hw_r_10(max_k - min_k, 0);
    std::vector<uint64_t> hw_r_10_acc(max_k - min_k, 0);
    
    std::vector<uint64_t> hw_r_2(max_k - min_k, 0);
    std::vector<uint64_t> hw_r_2_acc(max_k - min_k, 0);
    
    std::vector<uint64_t> hw_r_0(max_k - min_k, 0);
    std::vector<uint64_t> hw_r_0_acc(max_k - min_k, 0);

    std::vector<uint64_t> hw_r_4(max_k - min_k, 0);
    std::vector<uint64_t> hw_r_4_acc(max_k - min_k, 0);

    std::vector<uint64_t> hw_r_6(max_k - min_k, 0);
    std::vector<uint64_t> hw_r_6_acc(max_k - min_k, 0);

    std::vector<uint64_t> hw_r_8(max_k - min_k, 0);
    std::vector<uint64_t> hw_r_8_acc(max_k - min_k, 0);
    ////



    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, max_shot, max_shot);
    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }

    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        if(world_rank == 0 ){
                std::cout << "Started bucket " << k+min_k<< " - probability  = " << simulator.prob_k[k+min_k] << " - #iterations: " << simulator.shots_per_k[k+min_k]<< std::endl;

        }
        
        std::vector<qrc::DecodingGraph::Vertex*> A; // Actual Fast Group
        std::vector<qrc::DecodingGraph::Edge*> A_edge; // Actual Fast Group as edges
        std::vector<qrc::DecodingGraph::Vertex*> B; // Predicted Fast Group
        std::vector<qrc::DecodingGraph::Vertex*> C; // Boundary connected group
        qrc::DecoderShotResult decoder_results;
        uint hw = 0;
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);

            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);

                hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
                //Creating B group
                predecoder->set_the_priorities(syndrome);

                decoder_results = decoder->decode_error(syndrome);
                for(const auto& matching: decoder_results.matching){
                    qrc::DecodingGraph::Vertex* v1 = predecoder->decoding_graph.get_vertex(matching.first);
                    qrc::DecodingGraph::Vertex* v2 = predecoder->decoding_graph.get_vertex(matching.second);
                    qrc::DecodingGraph::Edge* edge = predecoder->decoding_graph.get_edge(v1, v2);
                    if(edge != nullptr){
                        //Creating A group
                        qpd::addUniqueVertex(A, v1);
                        qpd::addUniqueVertex(A, v2);
                        qpd::addUniqueEdge(A_edge, edge);
                        
                    }
                }
                if( predecoder->number_of_edges < hw){
                    if( 10<(hw-predecoder->number_of_edges)){
                        hw_r_high[k]++;
                    }
                    else if( 10==(hw-predecoder->number_of_edges) || 9==(hw-predecoder->number_of_edges)){
                        hw_r_10[k]++;
                    }
                    else if( 8==(hw-predecoder->number_of_edges) || 7==(hw-predecoder->number_of_edges)){
                        hw_r_8[k]++;
                    }
                    else if( 6==(hw-predecoder->number_of_edges) || 5==(hw-predecoder->number_of_edges)){
                        hw_r_6[k]++;
                    }
                    else if( 4==(hw-predecoder->number_of_edges) || 3==(hw-predecoder->number_of_edges)){
                        hw_r_4[k]++;
                    }
                    else if( 2==(hw-predecoder->number_of_edges) || 1==(hw-predecoder->number_of_edges)){
                        hw_r_2[k]++;
                    }
                }
                else if(predecoder->number_of_edges == hw){
                    hw_r_0[k]++;
                }
                //std::cout << world_rank << ":&&&" <<hw_r_0[k]<<std::endl;

                for(uint p=0;p<6;p++){
                    auto B = predecoder->pf_groups[p];
                    qpd::CountResult c2 = qpd::countMembersInAB(A,B);
                    qpd::CountResult c = {0,0,0};
                    c = qpd::countMembersInAB_edges(A_edge,predecoder->pf_pairs[p]);
                    
                    fp_t CR2 = A.size() !=0 ? ((double)c2.countBothAB/(c2.countBothAB+c2.countANotInB)):0;

                    fp_t CR = A_edge.size() !=0 ? ((double)c.countBothAB/(c.countBothAB+c.countANotInB)):0;
                    // if((c.countBothAB + c.countANotInB) != A.size()){
                    //     std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. A - B + AandB != A";
                    // }
                    if(1 < CR){
                        std::cout << "QQQ";
                    }

                    fp_t AR2 = B.size()!=0 ?((double)c2.countBothAB/(c2.countBothAB+c2.countBNotInA)):0;
                    fp_t AR = predecoder->pf_pairs[p].size()!=0 ?((double)c.countBothAB/(c.countBothAB+c.countBNotInA)):0;

                    if(1 < AR){
                        std::cout << "XXX";
                    }
                    
                    // if((c.countBothAB + c.countBNotInA) != B.size()){
                    //     std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. B - A + AandB != B";
                    // }
                    
                    if(AR == 1){
                        times_AR_is_1[p][k]++;
                    }
                    AR_avg[p][k] += AR;
                    if( CR == 1){
                        times_CR_is_1[p][k]++;
                    }
                    CR_avg[p][k] += CR;

                    if(AR2 == 1){
                        times_AR2_is_1[p][k]++;
                    }
                    AR2_avg[p][k] += AR2;
                    if( CR2 == 1){
                        times_CR2_is_1[p][k]++;
                    }
                    CR2_avg[p][k] += CR2;
                    

                    if( 10<B.size()){
                        PG_high[p][k]++;
                    }
                    else if( 10==B.size() || 9==B.size()){
                        PG_10[p][k]++;
                    }
                    else if( 8==B.size()|| 7==B.size()){
                        PG_8[p][k]++;
                    }
                    else if( 6==B.size() || 5==B.size()){
                        PG_6[p][k]++;
                    }
                    else if( 4==B.size()|| 3==B.size()){
                        PG_4[p][k]++;
                    }
                    else if( 2==B.size()|| 1==B.size()){
                        PG_2[p][k]++;
                    }

                    uint predecoding_size = ((uint)predecoder->pf_pairs[p].size())*2;

                    if(1 == c.countBNotInA){
                        times_1_incorrect_pair[p][k]++;
                        if((hw-predecoding_size)==9 || (hw-predecoding_size)==10){
                            inc_pair_1_hw_10[p][k]++;
                        }
                        else if(10<(hw-predecoding_size)){
                            inc_pair_1_hhw[p][k]++;
                        }
                    }
                    else if(2 == c.countBNotInA){
                        times_2_incorrect_pair[p][k]++;
                        if((hw-predecoding_size)==7 || (hw-predecoding_size)==8){
                            inc_pair_2_hw_8[p][k]++;
                        }
                        else if((hw-predecoding_size)==9 || (hw-predecoding_size)==10){
                            inc_pair_2_hw_10[p][k]++;
                        }
                        else if(10<(hw-predecoding_size)){
                            inc_pair_2_hhw[p][k]++;
                        }
                    }
                    else if(3 == c.countBNotInA){
                        times_3_incorrect_pair[p][k]++;
                        if((hw-predecoding_size)==5 || (hw-predecoding_size)==6){
                            inc_pair_3_hw_6[p][k]++;
                        }
                        else if((hw-predecoding_size)==7 || (hw-predecoding_size)==8){
                            inc_pair_3_hw_8[p][k]++;
                        }
                        else if((hw-predecoding_size)==9 || (hw-predecoding_size)==10){
                            inc_pair_3_hw_10[p][k]++;
                        }
                        else if(10<(hw-predecoding_size)){
                            inc_pair_3_hhw[p][k]++;
                        }

                    }
                    else if(4 == c.countBNotInA){
                        times_4_incorrect_pair[p][k]++;
                    }
                    else if(4 < c.countBNotInA){
                        more_than_4_incorrect_pair[p][k]++;
                    }
                }   
            }
            shots -= shots_this_round;
        }

    }
     
    
    // Flatten the 2D vectors into 1D vectors

    std::vector<uint64_t> flattened_times_AR_is_1;
    std::vector<uint64_t> flattened_times_CR_is_1;
    std::vector<fp_t> flattened_AR_avg;
    std::vector<fp_t> flattened_CR_avg;

    std::vector<uint64_t> flattened_times_AR2_is_1;
    std::vector<uint64_t> flattened_times_CR2_is_1;
    std::vector<fp_t> flattened_AR2_avg;
    std::vector<fp_t> flattened_CR2_avg;

    std::vector<uint64_t> flattened_PG_high;
    std::vector<uint64_t> flattened_PG_10;
    std::vector<uint64_t> flattened_PG_2;
    std::vector<uint64_t> flattened_PG_0;
    std::vector<uint64_t> flattened_PG_4;
    std::vector<uint64_t> flattened_PG_6;
    std::vector<uint64_t> flattened_PG_8;
    std::vector<uint64_t> flattened_times_1_incorrect_pair;
    std::vector<uint64_t> flattened_times_2_incorrect_pair;
    std::vector<uint64_t> flattened_times_3_incorrect_pair;
    std::vector<uint64_t> flattened_times_4_incorrect_pair;
    std::vector<uint64_t> flattened_more_than_4_incorrect_pair;

    std::vector<uint64_t> flattened_inc_pair_3_hw_6;
    std::vector<uint64_t> flattened_inc_pair_3_hw_8;
    std::vector<uint64_t> flattened_inc_pair_3_hw_10;
    std::vector<uint64_t> flattened_inc_pair_3_hhw;
    std::vector<uint64_t> flattened_inc_pair_2_hw_8;
    std::vector<uint64_t> flattened_inc_pair_2_hw_10;
    std::vector<uint64_t> flattened_inc_pair_2_hhw;
    std::vector<uint64_t> flattened_inc_pair_1_hw_10;
    std::vector<uint64_t> flattened_inc_pair_1_hhw;

    for (const auto& row : times_AR_is_1) {
        flattened_times_AR_is_1.insert(flattened_times_AR_is_1.end(), row.begin(), row.end());
    }

    for (const auto& row : times_CR_is_1) {
        flattened_times_CR_is_1.insert(flattened_times_CR_is_1.end(), row.begin(), row.end());
    }

    for (const auto& row : AR_avg) {
        flattened_AR_avg.insert(flattened_AR_avg.end(), row.begin(), row.end());
    }

    for (const auto& row : CR_avg) {
        flattened_CR_avg.insert(flattened_CR_avg.end(), row.begin(), row.end());
    }

    for (const auto& row : times_AR2_is_1) {
        flattened_times_AR2_is_1.insert(flattened_times_AR2_is_1.end(), row.begin(), row.end());
    }

    for (const auto& row : times_CR2_is_1) {
        flattened_times_CR2_is_1.insert(flattened_times_CR2_is_1.end(), row.begin(), row.end());
    }

    for (const auto& row : AR2_avg) {
        flattened_AR2_avg.insert(flattened_AR2_avg.end(), row.begin(), row.end());
    }

    for (const auto& row : CR2_avg) {
        flattened_CR2_avg.insert(flattened_CR2_avg.end(), row.begin(), row.end());
    }

    for (const auto& row : PG_high) {
        flattened_PG_high.insert(flattened_PG_high.end(), row.begin(), row.end());
    }

    for (const auto& row : PG_10) {
        flattened_PG_10.insert(flattened_PG_10.end(), row.begin(), row.end());
    }

    for (const auto& row : PG_2) {
        flattened_PG_2.insert(flattened_PG_2.end(), row.begin(), row.end());
    }

    for (const auto& row : PG_4) {
        flattened_PG_4.insert(flattened_PG_4.end(), row.begin(), row.end());
    }

    for (const auto& row : PG_6) {
        flattened_PG_6.insert(flattened_PG_6.end(), row.begin(), row.end());
    }

    for (const auto& row : PG_8) {
        flattened_PG_8.insert(flattened_PG_8.end(), row.begin(), row.end());
    }

    for (const auto& row : PG_0) {
        flattened_PG_0.insert(flattened_PG_0.end(), row.begin(), row.end());
    }

    for (const auto& row : times_1_incorrect_pair) {
        flattened_times_1_incorrect_pair.insert(flattened_times_1_incorrect_pair.end(), row.begin(), row.end());
    }

    for (const auto& row : times_2_incorrect_pair) {
        flattened_times_2_incorrect_pair.insert(flattened_times_2_incorrect_pair.end(), row.begin(), row.end());
    }

    for (const auto& row : times_3_incorrect_pair) {
        flattened_times_3_incorrect_pair.insert(flattened_times_3_incorrect_pair.end(), row.begin(), row.end());
    }

    for (const auto& row : times_4_incorrect_pair) {
        flattened_times_4_incorrect_pair.insert(flattened_times_4_incorrect_pair.end(), row.begin(), row.end());
    }

    for (const auto& row : more_than_4_incorrect_pair) {
        flattened_more_than_4_incorrect_pair.insert(flattened_more_than_4_incorrect_pair.end(), row.begin(), row.end());
    }

    // for incorrect&no space experiment
    for (const auto& row : inc_pair_3_hw_6) {
        flattened_inc_pair_3_hw_6.insert(flattened_inc_pair_3_hw_6.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_3_hw_8) {
        flattened_inc_pair_3_hw_8.insert(flattened_inc_pair_3_hw_8.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_3_hw_10) {
        flattened_inc_pair_3_hw_10.insert(flattened_inc_pair_3_hw_10.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_3_hhw) {
        flattened_inc_pair_3_hhw.insert(flattened_inc_pair_3_hhw.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_2_hw_8) {
        flattened_inc_pair_2_hw_8.insert(flattened_inc_pair_2_hw_8.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_2_hw_10) {
        flattened_inc_pair_2_hw_10.insert(flattened_inc_pair_2_hw_10.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_2_hhw) {
        flattened_inc_pair_2_hhw.insert(flattened_inc_pair_2_hhw.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_1_hw_10) {
        flattened_inc_pair_1_hw_10.insert(flattened_inc_pair_1_hw_10.end(), row.begin(), row.end());
    }
    for (const auto& row : inc_pair_1_hhw) {
        flattened_inc_pair_1_hhw.insert(flattened_inc_pair_1_hhw.end(), row.begin(), row.end());
    }

    // Perform the reduction for the remaining flattened vectors

    std::vector<uint64_t> flattened_times_AR_is_1_acc(flattened_times_AR_is_1.size());
    std::vector<uint64_t> flattened_times_CR_is_1_acc(flattened_times_CR_is_1.size());
    std::vector<fp_t> flattened_AR_avg_acc(flattened_AR_avg.size());
    std::vector<fp_t> flattened_CR_avg_acc(flattened_CR_avg.size());
    
    std::vector<uint64_t> flattened_times_AR2_is_1_acc(flattened_times_AR2_is_1.size());
    std::vector<uint64_t> flattened_times_CR2_is_1_acc(flattened_times_CR2_is_1.size());
    std::vector<fp_t> flattened_AR2_avg_acc(flattened_AR2_avg.size());
    std::vector<fp_t> flattened_CR2_avg_acc(flattened_CR2_avg.size());
    
    std::vector<uint64_t> flattened_PG_high_acc(flattened_PG_high.size());
    std::vector<uint64_t> flattened_PG_10_acc(flattened_PG_10.size());
    std::vector<uint64_t> flattened_PG_2_acc(flattened_PG_2.size());
    std::vector<uint64_t> flattened_PG_0_acc(flattened_PG_0.size());
    std::vector<uint64_t> flattened_PG_4_acc(flattened_PG_4.size());
    std::vector<uint64_t> flattened_PG_6_acc(flattened_PG_6.size());
    std::vector<uint64_t> flattened_PG_8_acc(flattened_PG_8.size());
    std::vector<uint64_t> flattened_times_1_incorrect_pair_acc(flattened_times_1_incorrect_pair.size());
    std::vector<uint64_t> flattened_times_2_incorrect_pair_acc(flattened_times_2_incorrect_pair.size());
    std::vector<uint64_t> flattened_times_3_incorrect_pair_acc(flattened_times_3_incorrect_pair.size());
    std::vector<uint64_t> flattened_times_4_incorrect_pair_acc(flattened_times_4_incorrect_pair.size());
    std::vector<uint64_t> flattened_more_than_4_incorrect_pair_acc(flattened_more_than_4_incorrect_pair.size());

    std::vector<uint64_t> flattened_inc_pair_3_hw_6_acc(flattened_inc_pair_3_hw_6.size());
    std::vector<uint64_t> flattened_inc_pair_3_hw_8_acc(flattened_inc_pair_3_hw_8.size());
    std::vector<uint64_t> flattened_inc_pair_3_hw_10_acc(flattened_inc_pair_3_hw_10.size());
    std::vector<uint64_t> flattened_inc_pair_3_hhw_acc(flattened_inc_pair_3_hhw.size());
    std::vector<uint64_t> flattened_inc_pair_2_hw_8_acc(flattened_inc_pair_2_hw_8.size());
    std::vector<uint64_t> flattened_inc_pair_2_hw_10_acc(flattened_inc_pair_2_hw_10.size());
    std::vector<uint64_t> flattened_inc_pair_2_hhw_acc(flattened_inc_pair_2_hhw.size());
    std::vector<uint64_t> flattened_inc_pair_1_hw_10_acc(flattened_inc_pair_1_hw_10.size());
    std::vector<uint64_t> flattened_inc_pair_1_hhw_acc(flattened_inc_pair_1_hhw.size());

    // Perform the reduction for the remaining flattened vectors

    MPI_Reduce(flattened_times_AR_is_1.data(), flattened_times_AR_is_1_acc.data(), flattened_times_AR_is_1.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(flattened_times_CR_is_1.data(), flattened_times_CR_is_1_acc.data(), flattened_times_CR_is_1.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(flattened_AR_avg.data(), flattened_AR_avg_acc.data(), flattened_AR_avg.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(flattened_CR_avg.data(), flattened_CR_avg_acc.data(), flattened_CR_avg.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_times_AR2_is_1.data(), flattened_times_AR2_is_1_acc.data(), flattened_times_AR2_is_1.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(flattened_times_CR2_is_1.data(), flattened_times_CR2_is_1_acc.data(), flattened_times_CR2_is_1.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(flattened_AR2_avg.data(), flattened_AR2_avg_acc.data(), flattened_AR2_avg.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(flattened_CR2_avg.data(), flattened_CR2_avg_acc.data(), flattened_CR2_avg.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_PG_high.data(), flattened_PG_high_acc.data(), flattened_PG_high.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_PG_10.data(), flattened_PG_10_acc.data(), flattened_PG_10.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(flattened_PG_0.data(), flattened_PG_0_acc.data(), flattened_PG_0.size(),
                MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_PG_2.data(), flattened_PG_2_acc.data(), flattened_PG_2.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_PG_4.data(), flattened_PG_4_acc.data(), flattened_PG_4.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_PG_6.data(), flattened_PG_6_acc.data(), flattened_PG_6.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_PG_8.data(), flattened_PG_8_acc.data(), flattened_PG_8.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(flattened_times_1_incorrect_pair.data(), flattened_times_1_incorrect_pair_acc.data(), flattened_times_1_incorrect_pair.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_times_2_incorrect_pair.data(), flattened_times_2_incorrect_pair_acc.data(), flattened_times_2_incorrect_pair.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_times_3_incorrect_pair.data(), flattened_times_3_incorrect_pair_acc.data(), flattened_times_3_incorrect_pair.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_times_4_incorrect_pair.data(), flattened_times_4_incorrect_pair_acc.data(), flattened_times_4_incorrect_pair.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_more_than_4_incorrect_pair.data(), flattened_more_than_4_incorrect_pair_acc.data(), flattened_more_than_4_incorrect_pair.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    //
    MPI_Reduce(flattened_inc_pair_3_hw_6.data(), flattened_inc_pair_3_hw_6_acc.data(), flattened_inc_pair_3_hw_6.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_inc_pair_3_hw_8.data(), flattened_inc_pair_3_hw_8_acc.data(), flattened_inc_pair_3_hw_8.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
 
    MPI_Reduce(flattened_inc_pair_3_hw_10.data(), flattened_inc_pair_3_hw_10_acc.data(), flattened_inc_pair_3_hw_10.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(flattened_inc_pair_3_hhw.data(), flattened_inc_pair_3_hhw_acc.data(), flattened_inc_pair_3_hhw.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(flattened_inc_pair_2_hw_8.data(), flattened_inc_pair_2_hw_8_acc.data(), flattened_inc_pair_2_hw_8.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_inc_pair_2_hw_10.data(), flattened_inc_pair_2_hw_10_acc.data(), flattened_inc_pair_2_hw_10.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
 
    MPI_Reduce(flattened_inc_pair_2_hhw.data(), flattened_inc_pair_2_hhw_acc.data(), flattened_inc_pair_2_hhw.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(flattened_inc_pair_1_hw_10.data(), flattened_inc_pair_1_hw_10_acc.data(), flattened_inc_pair_1_hw_10.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Reduce(flattened_inc_pair_1_hhw.data(), flattened_inc_pair_1_hhw_acc.data(), flattened_inc_pair_1_hhw.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // Reshape the remaining flattened vectors in the same way
    for (int i = 0; i < times_AR_is_1_acc.size(); i++) {
        times_AR_is_1_acc[i] = std::vector<uint64_t>(flattened_times_AR_is_1_acc.begin() + i * (max_k - min_k),
                                                 flattened_times_AR_is_1_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < times_CR_is_1_acc.size(); i++) {
        times_CR_is_1_acc[i] = std::vector<uint64_t>(flattened_times_CR_is_1_acc.begin() + i * (max_k - min_k),
                                                 flattened_times_CR_is_1_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < AR_avg_acc.size(); i++) {
        AR_avg_acc[i] = std::vector<fp_t>(flattened_AR_avg_acc.begin() + i * (max_k - min_k),
                                                 flattened_AR_avg_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < CR_avg_acc.size(); i++) {
        CR_avg_acc[i] = std::vector<fp_t>(flattened_CR_avg_acc.begin() + i * (max_k - min_k),
                                                 flattened_CR_avg_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < times_AR2_is_1_acc.size(); i++) {
        times_AR2_is_1_acc[i] = std::vector<uint64_t>(flattened_times_AR2_is_1_acc.begin() + i * (max_k - min_k),
                                                 flattened_times_AR2_is_1_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < times_CR2_is_1_acc.size(); i++) {
        times_CR2_is_1_acc[i] = std::vector<uint64_t>(flattened_times_CR2_is_1_acc.begin() + i * (max_k - min_k),
                                                 flattened_times_CR2_is_1_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < AR2_avg_acc.size(); i++) {
        AR2_avg_acc[i] = std::vector<fp_t>(flattened_AR2_avg_acc.begin() + i * (max_k - min_k),
                                                 flattened_AR2_avg_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < CR2_avg_acc.size(); i++) {
        CR2_avg_acc[i] = std::vector<fp_t>(flattened_CR2_avg_acc.begin() + i * (max_k - min_k),
                                                 flattened_CR2_avg_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < PG_high_acc.size(); i++) {
        PG_high_acc[i] = std::vector<uint64_t>(flattened_PG_high_acc.begin() + i * (max_k - min_k),
                                            flattened_PG_high_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < PG_10_acc.size(); i++) {
        PG_10_acc[i] = std::vector<uint64_t>(flattened_PG_10_acc.begin() + i * (max_k - min_k),
                                            flattened_PG_10_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < PG_2_acc.size(); i++) {
        PG_2_acc[i] = std::vector<uint64_t>(flattened_PG_2_acc.begin() + i * (max_k - min_k),
                                            flattened_PG_2_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < PG_0_acc.size(); i++) {
        PG_0_acc[i] = std::vector<uint64_t>(flattened_PG_0_acc.begin() + i * (max_k - min_k),
                                            flattened_PG_0_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < PG_4_acc.size(); i++) {
        PG_4_acc[i] = std::vector<uint64_t>(flattened_PG_4_acc.begin() + i * (max_k - min_k),
                                            flattened_PG_4_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < PG_6_acc.size(); i++) {
        PG_6_acc[i] = std::vector<uint64_t>(flattened_PG_6_acc.begin() + i * (max_k - min_k),
                                            flattened_PG_6_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < PG_8_acc.size(); i++) {
        PG_8_acc[i] = std::vector<uint64_t>(flattened_PG_8_acc.begin() + i * (max_k - min_k),
                                            flattened_PG_8_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < times_1_incorrect_pair_acc.size(); i++) {
        times_1_incorrect_pair_acc[i] = std::vector<uint64_t>(flattened_times_1_incorrect_pair_acc.begin() + i * (max_k - min_k),
                                            flattened_times_1_incorrect_pair_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < times_2_incorrect_pair_acc.size(); i++) {
        times_2_incorrect_pair_acc[i] = std::vector<uint64_t>(flattened_times_2_incorrect_pair_acc.begin() + i * (max_k - min_k),
                                            flattened_times_2_incorrect_pair_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < times_3_incorrect_pair_acc.size(); i++) {
        times_3_incorrect_pair_acc[i] = std::vector<uint64_t>(flattened_times_3_incorrect_pair_acc.begin() + i * (max_k - min_k),
                                            flattened_times_3_incorrect_pair_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < times_4_incorrect_pair_acc.size(); i++) {
        times_4_incorrect_pair_acc[i] = std::vector<uint64_t>(flattened_times_4_incorrect_pair_acc.begin() + i * (max_k - min_k),
                                            flattened_times_4_incorrect_pair_acc.begin() + (i + 1) * (max_k - min_k));
    }

    for (int i = 0; i < more_than_4_incorrect_pair_acc.size(); i++) {
        more_than_4_incorrect_pair_acc[i] = std::vector<uint64_t>(flattened_more_than_4_incorrect_pair_acc.begin() + i * (max_k - min_k),
                                            flattened_more_than_4_incorrect_pair_acc.begin() + (i + 1) * (max_k - min_k));
    }
    //
    for (int i = 0; i < inc_pair_3_hw_6_acc.size(); i++) {
        inc_pair_3_hw_6_acc[i] = std::vector<uint64_t>(flattened_inc_pair_3_hw_6_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_3_hw_6_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_3_hw_8_acc.size(); i++) {
        inc_pair_3_hw_8_acc[i] = std::vector<uint64_t>(flattened_inc_pair_3_hw_8_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_3_hw_8_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_3_hw_10_acc.size(); i++) {
        inc_pair_3_hw_10_acc[i] = std::vector<uint64_t>(flattened_inc_pair_3_hw_10_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_3_hw_10_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_3_hhw_acc.size(); i++) {
        inc_pair_3_hhw_acc[i] = std::vector<uint64_t>(flattened_inc_pair_3_hhw_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_3_hhw_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_2_hw_8_acc.size(); i++) {
        inc_pair_2_hw_8_acc[i] = std::vector<uint64_t>(flattened_inc_pair_2_hw_8_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_2_hw_8_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_2_hw_10_acc.size(); i++) {
        inc_pair_2_hw_10_acc[i] = std::vector<uint64_t>(flattened_inc_pair_2_hw_10_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_2_hw_10_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_2_hhw_acc.size(); i++) {
        inc_pair_2_hhw_acc[i] = std::vector<uint64_t>(flattened_inc_pair_2_hhw_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_2_hhw_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_1_hw_10_acc.size(); i++) {
        inc_pair_1_hw_10_acc[i] = std::vector<uint64_t>(flattened_inc_pair_1_hw_10_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_1_hw_10_acc.begin() + (i + 1) * (max_k - min_k));
    }
    for (int i = 0; i < inc_pair_1_hhw_acc.size(); i++) {
        inc_pair_1_hhw_acc[i] = std::vector<uint64_t>(flattened_inc_pair_1_hhw_acc.begin() + i * (max_k - min_k),
                                            flattened_inc_pair_1_hhw_acc.begin() + (i + 1) * (max_k - min_k));
    }
    

    MPI_Reduce(&hw_r_0[0], &hw_r_0_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hw_r_2[0], &hw_r_2_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hw_r_4[0], &hw_r_4_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hw_r_6[0], &hw_r_6_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hw_r_8[0], &hw_r_8_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hw_r_10[0], &hw_r_10_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&hw_r_high[0], &hw_r_high_acc[0], max_k-min_k, MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);



    if(world_rank == 0){ 
        std::cout << "________________________________________________________________" <<std::endl;
        std::cout << std::setw(18) <<  " k (p) | "; 
        std::cout <<  std::setw(16) << " P(HW-|B| = 0) | "; 
        std::cout <<  std::setw(16) << " P(HW-|B| = 2) | "; 
        std::cout <<  std::setw(16) << " P(HW-|B| = 4) | "; 
        std::cout <<  std::setw(16) << " P(HW-|B| = 6) | "; 
        std::cout <<  std::setw(16) << " P(HW-|B| = 8) | "; 
        std::cout <<  std::setw(16) << " P(HW-|B| = 10) | "; 
        std::cout <<  std::setw(16) << " P(10 < HW-|B|) | " << std::endl;
        for ( uint i = 0; i< max_k-min_k;  i++){
            std::string printing_shots = std::to_string(simulator.shots_per_k[i+min_k]);
            if(1'000<=simulator.shots_per_k[i+min_k]){
                printing_shots = number_of_shots_to_string(simulator.shots_per_k[i+min_k]);
            }
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(10) <<((double)hw_r_0_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << hw_r_0_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)hw_r_2_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << hw_r_2_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)hw_r_4_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << hw_r_4_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)hw_r_6_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << hw_r_6_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)hw_r_8_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << hw_r_8_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)hw_r_10_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << hw_r_10_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)hw_r_high_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << hw_r_high_acc[i] << "/" << printing_shots <<") | "
            << std::endl;
            
        }
        

        
        for(int p = 0;p < 1;p++){
            std::cout << "######################PRIORITY GROUP " << p << "######################" << std::endl;
            std::cout << "Distance: "<< distance << " - #Buckets Max: " << max_k << " Min: "<< min_k <<" - Physical ER: " << physcial_error << std::endl;
            std::cout << "________________________________________________________________" <<std::endl;
            std::cout << std::setw(18) <<  " k (p) | ";
            std::cout << std::setw(20) << " P(AR = 1) | ";
            std::cout << std::setw(20) << " P(CR = 1) | ";
            std::cout << std::setw(20) << " P(AR2 = 1) | ";
            std::cout << std::setw(20) << " P(CR2 = 1) | "<< std::endl;
            for ( uint i = 0; i< max_k-min_k;  i++){
                std::string printing_shots = std::to_string(simulator.shots_per_k[i+min_k]);
                if(1'000<=simulator.shots_per_k[i+min_k]){
                    printing_shots = number_of_shots_to_string(simulator.shots_per_k[i+min_k]);
                }
                std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
                std::cout << std::setw(10) <<((double)times_AR_is_1_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_AR_is_1_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)times_CR_is_1_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_CR_is_1_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)times_AR2_is_1_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_AR2_is_1_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)times_CR2_is_1_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_CR2_is_1_acc[p][i] << "/" << printing_shots <<") | "<< std::endl;
               
               
        
            }
            std::cout << "________________________________________________________________" <<std::endl;
            std::cout << std::setw(18) <<  " k (p) | ";
            std::cout << std::setw(20) << " AVG AR | ";
            std::cout << std::setw(20) << " AVG CR | ";
            std::cout << std::setw(20) << " AVG AR2 | ";
            std::cout << std::setw(20) << " AVG CR2 | "<< std::endl;
            
            for ( uint i = 0; i< max_k-min_k;  i++){
                std::string printing_shots = std::to_string(simulator.shots_per_k[i+min_k]);
                if(1'000<=simulator.shots_per_k[i+min_k]){
                    printing_shots = number_of_shots_to_string(simulator.shots_per_k[i+min_k]);
                }
                std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
                std::cout << std::setw(10) <<((double)AR_avg_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << AR_avg_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)CR_avg_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << CR_avg_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)AR2_avg_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << AR2_avg_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)CR2_avg_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << CR2_avg_acc[p][i] << "/" << printing_shots <<") | "<< std::endl;
               
               
        
            }
            std::cout << "________________________________________________________________" <<std::endl;            
            std::cout << std::setw(18) <<  " k (p) | "; 
            std::cout <<  std::setw(21) << " P(|PG| = 2) | "; 
            std::cout <<  std::setw(21) << " P(|PG| = 4) | "; 
            std::cout <<  std::setw(21) << " P(|PG| = 6) | "; 
            std::cout <<  std::setw(21) << " P(|PG| = 8) | "; 
            std::cout <<  std::setw(21) << " P(|PG| = 10) | "; 
            std::cout <<  std::setw(21) << " P(10 < |PG|) | " << std::endl;
            for ( uint i = 0; i< max_k-min_k;  i++){
                std::string printing_shots = std::to_string(simulator.shots_per_k[i+min_k]);
                if(1'000<=simulator.shots_per_k[i+min_k]){
                    printing_shots = number_of_shots_to_string(simulator.shots_per_k[i+min_k]);
                }
                std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
                std::cout << std::setw(10) <<((double)PG_2_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << PG_2_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)PG_4_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << PG_4_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)PG_6_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << PG_6_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)PG_8_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << PG_8_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)PG_10_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << PG_10_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(10) <<((double)PG_high_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << PG_high_acc[p][i] << "/" << printing_shots <<") | "
                << std::endl;
                
            }

            std::cout << "________________________________________________________________" <<std::endl;
            std::cout << std::setw(18) <<  " k (p) | "; 
            std::cout << std::setw(6) << " P(|InCorr. Pairs|=1) | "; 
            std::cout << std::setw(6) << " P(|InCorr. Pairs|=2) |";
            std::cout << std::setw(6) << " P(|InCorr. Pairs|=3) |";
            std::cout << std::setw(6) << " P(|InCorr. Pairs|=4) |";
            std::cout << std::setw(6) << " P(4<|InCorr. Pairs|)|" << std::endl;

            for ( uint i = 0; i< max_k-min_k;  i++){
                std::string printing_shots = std::to_string(simulator.shots_per_k[i+min_k]);
                if(1'000<=simulator.shots_per_k[i+min_k]){
                    printing_shots = number_of_shots_to_string(simulator.shots_per_k[i+min_k]);
                }
                std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
                std::cout << std::setw(14) <<((double)times_1_incorrect_pair_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_1_incorrect_pair_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)times_2_incorrect_pair_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_2_incorrect_pair_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)times_3_incorrect_pair_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_3_incorrect_pair_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)times_4_incorrect_pair_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << times_4_incorrect_pair_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)more_than_4_incorrect_pair_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << more_than_4_incorrect_pair_acc[p][i] << "/" << printing_shots <<") | "
                << std::endl;
                
            }

            std::cout << "________________________________________________________________" <<std::endl;
            std::cout << std::setw(18) <<  " k (p) | "; 
            std::cout << std::setw(12) << " P(IP=1,RHW=10) | "; 
            std::cout << std::setw(12) << " P(IP=1,10<RHW) | ";
            std::cout << std::setw(12) << " P(IP=2,RHW=8) | ";
            std::cout << std::setw(12) << " P(IP=2,RHW=10) | ";
            std::cout << std::setw(12) << " P(IP=2,10<RHW) | ";
            std::cout << std::setw(12) << " P(IP=3,RHW=6) | ";
            std::cout << std::setw(12) << " P(IP=3,RHW=8) | ";
            std::cout << std::setw(12) << " P(IP=3,RHW=10) | ";
            std::cout << std::setw(12) << " P(IP=3,10<RHW) | " << std::endl;

            for ( uint i = 0; i< max_k-min_k;  i++){
                std::string printing_shots = std::to_string(simulator.shots_per_k[i+min_k]);
                if(1'000<=simulator.shots_per_k[i+min_k]){
                    printing_shots = number_of_shots_to_string(simulator.shots_per_k[i+min_k]);
                }
                std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
                std::cout << std::setw(14) <<((double)inc_pair_1_hw_10_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_1_hw_10_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_1_hhw_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_1_hhw_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_2_hw_8_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_2_hw_8_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_2_hw_10_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_2_hw_10_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_2_hhw_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_2_hhw_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_3_hw_6_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_3_hw_6_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_3_hw_8_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_3_hw_8_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_3_hw_10_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_3_hw_10_acc[p][i] << "/" << printing_shots <<") | ";
                std::cout << std::setw(14) <<((double)inc_pair_3_hhw_acc[p][i]/simulator.shots_per_k[i+min_k]) << "(" << inc_pair_3_hhw_acc[p][i] << "/" << printing_shots <<") | "
                << std::endl;
                
            }
        }
        
        
    
    }
    

}

void test_size_of_predecoding(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, std::string decoder_name){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);
    

    auto path_table = qrc::compute_path_table(predecoder->decoding_graph);

    std::vector<std::vector<uint64_t>> difference(max_k - min_k, std::vector<uint64_t>(100, 0));
    std::vector<std::vector<uint64_t>> difference_acc(max_k - min_k, std::vector<uint64_t>(100, 0));

    
    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot, true, max_shot,max_shot);
 
    if(world_rank == 0){
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }

    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        if(world_rank == 0 ){
                std::cout << "Started bucket " << k+min_k<< " - probability  = " << simulator.prob_k[k+min_k] << " - #iterations: " << simulator.shots_per_k[k+min_k]<< std::endl;

        }
        
        std::vector<qrc::DecodingGraph::Vertex*> A; // Actual Fast Group
        std::vector<qrc::DecodingGraph::Edge*> A_edge; // Actual Fast Group as edges
        std::vector<qrc::DecodingGraph::Vertex*> B; // Predicted Fast Group
        std::vector<qrc::DecodingGraph::Vertex*> C; // Boundary connected group
        qrc::DecoderShotResult decoder_results;
        uint hw = 0;
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);

            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);

                B = predecoder->get_fast_matching_group(syndrome);

                //Creating B group
                predecoder->set_the_priorities(syndrome);
                int size_diff = B.size()-predecoder->predecoding_vertices.size();
                if(size_diff<0){
                    std::cout << "ERROR! THIS MESSAGE SHOULD NOT BE PRINTED. Fast Group is smaller than the predecoding group" <<std::endl;
                }
                difference[k][size_diff]++;

            }
            shots -= shots_this_round;
        }
    }

        // Flatten the 2D vector into a 1D vector
    std::vector<uint64_t> flattened_difference;
    for (const auto& row : difference) {
        flattened_difference.insert(flattened_difference.end(), row.begin(), row.end());
    }

    // Perform the reduction on the flattened 1D vector
    std::vector<uint64_t> flattened_difference_acc(flattened_difference.size());
    MPI_Reduce(flattened_difference.data(), flattened_difference_acc.data(), flattened_difference.size(),
            MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    // Reshape the flattened 1D vector back into a 2D vector
    for (int i = 0; i < difference_acc.size(); i++) {
        difference_acc[i] = std::vector<uint64_t>(flattened_difference_acc.begin() + i * 100,
                                                flattened_difference_acc.begin() + (i + 1) * 100);
    }

    if(world_rank == 0){
        for(uint k =0;k<max_k-min_k;k++){
            std::cout << "K = " << k+min_k << " p " <<simulator.prob_k[k+min_k] <<"___________"<< std::endl;
            for(uint j =0;j<100;j++){
                if(difference_acc[k][j]){
                    std::cout << j <<": " <<difference_acc[k][j]<<std::endl;
                }
            }
        }
    }
}

void bb_ler_calc_promatch(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name, bool generalized, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string syndrome_clock_challenging_file, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // std::string syndrome_file = "../NFDecoder/data/challengingSyndromes_d13/challenging_d13.txt";
    if(world_rank == 0){
        std::cout << "round#" <<round_n <<" Promatch-A + " << decoder_name <<std::endl;
    }



    // for saving syndromes
    uint64_t round = 1;
    uint64_t round_t = 1;
    if (world_rank == 0 && save_syndromes){
        std::cout << syndrome_folder_name<< std::endl;
        // To check if we are reading correct files
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
    
    }

    if (world_rank == 0 && save_syndromes){
        std::cout << syndrome_clock_challenging_file << std::endl;
        // To check if we are reading correct files
        if (opendir(syndrome_clock_challenging_file.c_str()) == NULL) {
            if (mkdir(syndrome_clock_challenging_file.c_str(), 0777) == 0) {
                std::cout << "Directory created successfully\n";
            } else {
                std::cout << "Failed to create directory\n";
            }
        } 
        else {
            std::cout << "Directory already exists\n";
        }
    
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
    // if(distance == 13){
    //     max_k = 17;
    // } 
    // else if(distance == 11){
    //     max_k = 12;
    // }
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;
    qrc::benchmark::StatisticalResult statres;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_size*round_n+world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ, generalized);
    predecoder->set_decoder(decoder);

    //////
    std::vector<vector<vector<uint16_t>>> saved_hhw_syndromes;
    saved_hhw_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits;
    //////
    /// For saving time challenging syndromes
    std::vector<vector<vector<uint16_t>>> saved_timechallenging_syndromes;
    saved_timechallenging_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits_t;
    ////

    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, hshots_replc, lshots_rplac);

    uint64_t max_k_ = max_k;
    fp_t sum_normalized_prob = 0;
    for(uint k = 0; k < max_k-min_k; k++){
        sum_normalized_prob += simulator.prob_k[k+min_k];

        // To eliminate the buckets that their probability is lower than LER and accordinly has zero number of shots
        if(simulator.shots_per_k[k+min_k] == 0){
            max_k_ --;
        }
    }
    max_k = max_k_;

    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << " Physical Err = "<< physcial_error<<std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }
    fp_t prev_logical_error_rate = 0.0;

    
    std::array<fp_t,6> normalized_stages;
    normalized_stages.fill(0);

    std::array<fp_t,6> weighted_stages;
    weighted_stages.fill(0);

    ///////Clock Simulation
    uint64_t local_sum_cycle = 0;

    uint64_t local_max_cycle = 0;

    uint64_t local_max_round = 0;

    uint64_t local_sum_cycle_parallel = 0;

    uint64_t local_max_cycle_parallel = 0;


    fp_t final_average_cycle = 0;
    fp_t final_average_cycle_normalized = 0;
    fp_t final_average_cycle_parallel = 0;
    fp_t final_average_cycle_parallel_normalized = 0;

    uint64_t local_number_above_threshold_parallel = 0;
    uint64_t local_number_above_threshold = 0;

    uint64_t local_n_Astrea10_used = 0;
    uint64_t local_n_Astrea8_used = 0;
    uint64_t local_n_Astrea6_used = 0;

    fp_t total_Astrea10_used = 0;
    fp_t total_Astrea8_used = 0;
    fp_t total_Astrea6_used = 0;

    fp_t time_constrained_ler = 0;
    uint64_t local_number_above_threshold_parallel_prior = 0;


    //////Clock Simulation


    //////
    // std::vector<std::vector<uint16_t>> vecss = qpd::read_vector_of_vectors(syndrome_file);
    // std::cout << "Size of incorrect pairs = " << vecss.size() << std::endl;
    //////
    // uint k=0;
    // for(auto compressed_synd : vecss){
    //     k++;
    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t local_errors = 0;
        
        std::array<uint64_t, 6> local_stages;
        local_stages.fill(0);

        uint64_t local_times_lhw = 0;

        local_sum_cycle_parallel = 0;
        local_sum_cycle = 0;

        local_n_Astrea10_used = 0;
        local_n_Astrea8_used = 0;
        local_n_Astrea6_used = 0;


        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                uint hw = std::accumulate(syndrome.begin(), syndrome.begin()+predecoder->n_detectors,0);
                if(hw <= MAX_BF10_HW){
                    local_times_lhw++;
                }
            /////////
            // std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);

            /////////
                qrc::DecoderShotResult res;

                res = predecoder->adaptively_decode_error(syndrome);
                // res = predecoder->ensemble_decode_error(syndrome, 0);

                // if(res.matching.size() == 0){
                //     std::cout << "res.matching.size() == 0 " << std::endl;
                //     flipped_bits = qpd::syndrome_compressed(syndrome);
                //     continue;

                // }
                
                local_errors += res.is_logical_error;
                for(uint stg = 0; stg <6; stg++){
                    local_stages[stg] += predecoder->reached_stage[stg];
                }
                
                if(res.is_logical_error){
                    //std::cout << "Error" << k <<"-";
                    flipped_bits = qpd::syndrome_compressed(syndrome);
                    saved_hhw_syndromes[world_rank].push_back(flipped_bits);
    
                }

                if(predecoder->total_cycles > MAX_BF6_CYCLE){
                    local_number_above_threshold++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                if(predecoder->total_cycles_parallel_updating > MAX_BF6_CYCLE){
                    local_number_above_threshold_parallel++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                

                //Setting clock related varaibles
                if (local_max_cycle < predecoder->total_cycles){
                    local_max_cycle = predecoder->total_cycles;
                }
                if(local_max_cycle_parallel < predecoder->total_cycles_parallel_updating){
                    local_max_cycle_parallel = predecoder->total_cycles_parallel_updating;
                }
                local_sum_cycle += predecoder->total_cycles;
                local_sum_cycle_parallel += predecoder->total_cycles_parallel_updating;

                if(local_max_round < predecoder->number_of_rounds){
                    local_max_round = predecoder->number_of_rounds;
                }
                if(predecoder->MAX_BF_HW == MAX_BF10_HW && hw > MAX_BF10_HW){
                    local_n_Astrea10_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF8_HW && hw > MAX_BF10_HW){
                    local_n_Astrea8_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF6_HW && hw > MAX_BF10_HW){
                    local_n_Astrea6_used++;
                }
            }
            shots -= shots_this_round;
        }
        uint64_t fault_errors = 0;
        MPI_Allreduce(&local_errors, &fault_errors, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t times_lhw = 0;
        MPI_Allreduce(&local_times_lhw, &times_lhw, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t, 6> stages;
        stages.fill(0);
        MPI_Allreduce(&local_stages, &stages, 6, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle = 0;
        MPI_Allreduce(&local_max_cycle, &max_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle = 0;
        MPI_Allreduce(&local_sum_cycle, &sum_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle_parallel = 0;
        MPI_Allreduce(&local_max_cycle_parallel, &max_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle_parallel = 0;
        MPI_Allreduce(&local_sum_cycle_parallel, &sum_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold = 0;
        MPI_Allreduce(&local_number_above_threshold, &number_above_threshold, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold_parallel = 0;
        MPI_Allreduce(&local_number_above_threshold_parallel, &number_above_threshold_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_round = 0;
        MPI_Allreduce(&local_max_round, &max_round, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        
        uint64_t n_Astrea10_used = 0;
        MPI_Allreduce(&local_n_Astrea10_used, &n_Astrea10_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea8_used = 0;
        MPI_Allreduce(&local_n_Astrea8_used, &n_Astrea8_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea6_used = 0;
        MPI_Allreduce(&local_n_Astrea6_used, &n_Astrea6_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }

        if (fault_errors > 0 || ((number_above_threshold_parallel -  local_number_above_threshold_parallel_prior) > 0)) {
            fp_t failure_rate = ((fp_t)fault_errors/simulator.shots_per_k[k+min_k]);
        
            prev_logical_error_rate = statres.logical_error_rate;
            statres.logical_error_rate += (failure_rate*simulator.prob_k[k+min_k]);
            statres.n_logical_errors += fault_errors;
            fp_t time_failure_rate = (fp_t)(number_above_threshold_parallel -  local_number_above_threshold_parallel_prior)/simulator.shots_per_k[k+min_k];
            time_constrained_ler += (failure_rate*simulator.prob_k[k+min_k]) + (time_failure_rate*simulator.prob_k[k+min_k]);
            local_number_above_threshold_parallel_prior = number_above_threshold_parallel;

            if (world_rank == 0) {
                std::cout << "\tFailure rate = " << failure_rate << std::endl
                        << "\tNumber of errors = " << statres.n_logical_errors << std::endl;
                std::cout << "\tLogical error rate = " << statres.logical_error_rate << std::endl;
                std::cout << "\t\t Logical error rate timed = " << time_constrained_ler << std::endl;
            }
        }
        std::array<fp_t,6> stages_rate;
        stages_rate.fill(0);
        for(uint stg = 0; stg < 6; stg++){
            stages_rate[stg] = ((fp_t)stages[stg]/(simulator.shots_per_k[k+min_k]-times_lhw));
            normalized_stages[stg] += stages_rate[stg]*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
            weighted_stages[stg] += stages_rate[stg]*simulator.prob_k[k+min_k];
 
        }

        // std::cout << world_rank <<" sum_cycle"  <<sum_cycle<< " sum_cycle_parallel: " <<sum_cycle_parallel<<std::endl;

        final_average_cycle += ((fp_t)sum_cycle/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];        
        final_average_cycle_normalized += ((fp_t)sum_cycle/(simulator.shots_per_k[k+min_k]))*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        final_average_cycle_parallel += ((fp_t)sum_cycle_parallel/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];
        final_average_cycle_parallel_normalized += ((fp_t)sum_cycle_parallel/(simulator.shots_per_k[k+min_k]))*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        fp_t percentage_astrea10 = ((fp_t)n_Astrea10_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea10_used +=percentage_astrea10*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
        
        fp_t percentage_astrea8 = ((fp_t)n_Astrea8_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea8_used +=percentage_astrea8*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
        
        fp_t percentage_astrea6 = ((fp_t)n_Astrea6_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea6_used +=percentage_astrea6*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);


        // std::cout << "world_rank: "<<world_rank  << " final_average_cycle: " << final_average_cycle <<
        // " final_average_cycle_normalized: " << final_average_cycle_normalized << " final_average_cycle_parallel: " << final_average_cycle_parallel << 
        // " final_average_cycle_parallel_normalized: " << final_average_cycle_parallel_normalized << std::endl;
        
        if (world_rank == 0) {
            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\t\tS1: " << stages_rate[0] << " S2_1: " << stages_rate[1] << " S2_2: " << stages_rate[2] << " S3: " << stages_rate[3] << " S4_1: " << stages_rate[4] 
                            << " S4_2: " << stages_rate[5] << std::endl;
            std::cout << "\tNormalized Stage Rates:" << std::endl;
            std::cout << "\t\tS1: " << normalized_stages[0] << " S2_1: " << normalized_stages[1] << " S2_2: " << normalized_stages[2] << " S3: " << normalized_stages[3] << " S4_1: " << normalized_stages[4] 
                            << " S4_2: " << normalized_stages[5] << std::endl;
            std::cout << "\tStage Rates:" << std::endl;
            std::cout << "\t\tS1: " << weighted_stages[0] << " S2_1: " << weighted_stages[1] << " S2_2: " << weighted_stages[2] << " S3: " << weighted_stages[3] << " S4_1: " << weighted_stages[4] 
                            << " S4_2: " << weighted_stages[5] << std::endl;

            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\tNormalized:" << std::endl;
            std::cout << "\t\t Avg Cycles (UP) = " << final_average_cycle_parallel_normalized  << "(" << final_average_cycle_normalized <<")"<<std::endl;
            std::cout << "\tNot Normalized:" << std::endl;
            std::cout << "\t\t Avg Cycles(UP) = " << final_average_cycle_parallel << "(" << final_average_cycle <<")"<<std::endl;
            std::cout << "\tMAX Cycles (UP) = " << max_cycle_parallel  << "(" << max_cycle << ") Number above threshold(UP):" << number_above_threshold_parallel << "(" << number_above_threshold << ")"<< std::endl;
            std::cout << "\tMAX Round = " << max_round  << std::endl;
            std::cout << "\t ----------ASTREA PERCENTAGES-------"<< std::endl;
            std::cout << "\t This bucket" << std::endl;
            std::cout << "\t\tAstrea10: " << percentage_astrea10 << " Astrea8: " << percentage_astrea8 << " Astrea6: " << percentage_astrea6 << std::endl;
            std::cout << "\t Total(Normalized)" << std::endl;
            std::cout << "\t\tAstrea10: " << total_Astrea10_used << " Astrea8: " << total_Astrea8_used << " Astrea6: " << total_Astrea6_used << std::endl;

            std::cout << "_______________________________________________"<< std::endl;


        }

        if(int(saved_hhw_syndromes[world_rank].size()) != 0){
            std::string file_name = syndrome_folder_name + syndrome_file_name("sgen",distance,physcial_error,
            meas_er, max_shot, world_rank, round);
            qpd::write_vector_of_vectors(saved_hhw_syndromes[world_rank], file_name);
            saved_hhw_syndromes[world_rank].clear();
            round++;
        }

        if(int(saved_timechallenging_syndromes[world_rank].size()) != 0){
            std::string file_name = syndrome_clock_challenging_file + syndrome_file_name("t_sgen",distance,physcial_error,
            meas_er, max_shot, world_rank, round);
            qpd::write_vector_of_vectors(saved_timechallenging_syndromes[world_rank], file_name);
            saved_timechallenging_syndromes[world_rank].clear();
            round_t++;
        }
    }

}

void bb_decoder_ler_calc(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string decoder_name, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    uint round = 0;
    if(world_rank == 0){
        std::cout << "round#" <<round_n <<" " << decoder_name <<std::endl;
    }

    if(distance == 11){
        min_k = 7;
        max_k = 10;
    }
    else if(distance == 13){
        min_k = 9;
        max_k = 12;
    }

    // for saving syndromes
    if (world_rank == 0 && save_syndromes){
        std::cout << syndrome_folder_name<< std::endl;
        // To check if we are reading correct files
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

    if(world_rank == 0){
        std::cout << "Circuit created, ";
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::benchmark::StatisticalResult statres;
    qrc::Decoder* decoder;
    if(world_rank == 0){
        std::cout << "Simulator created, ";
    }

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    
    else if(decoder_name.find("ASTREAG") != std::string::npos){
        qrc::AstreaParams astreaG_param= {};
        astreaG_param.bfu_fetch_width = 2;
        astreaG_param.bfu_priority_queue_size = 8;
        astreaG_param.main_clock_frequency = 250e6;
        astreaG_param.bfu_compute_stages = 2;
        astreaG_param.n_registers = 2000;
        astreaG_param.use_mld = false;

        uint n_detectors_per_round = (distance*distance - 1)/2;
        uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(0.01*0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    

        decoder =  new  qrc::Astrea(circ,
                                            n_detectors_per_round,
                                            weight_filter_cutoff,
                                            astreaG_param);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    if(world_rank == 0){
        std::cout << "Decoder " << decoder_name << " created, ";
    }
    std::mt19937_64 rng(world_size*round_n+world_rank);

    std::vector<vector<vector<uint16_t>>> saved_hhw_syndromes;
    saved_hhw_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits;

    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot, true, hshots_replc, lshots_rplac);
    
    uint64_t max_k_ = max_k;
    for(uint k = 0; k < max_k-min_k; k++){
        if(simulator.shots_per_k[k+min_k] == 0){
            max_k_ --;
        }
    }
    max_k = max_k_;

    if(world_rank == 0){
        std::cout << "Simulator set the number of shots. " << std::endl;
    }
    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot  << " Physical Err = "<< physcial_error<<std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }
    fp_t prev_logical_error_rate = 0.0;
    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t local_errors = 0;
        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                // std::cout << "!!1" << std::endl;
                auto res = decoder->decode_error(syndrome);
                
                local_errors += res.is_logical_error;
                if(res.is_logical_error){
                    flipped_bits = qpd::syndrome_compressed(syndrome);
                    saved_hhw_syndromes[world_rank].push_back(flipped_bits);

                }
                
            }
            shots -= shots_this_round;
        }
        uint64_t fault_errors = 0;
        MPI_Allreduce(&local_errors, &fault_errors, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }

        if (fault_errors > 0) {
            fp_t failure_rate = ((fp_t)fault_errors/simulator.shots_per_k[k+min_k]);
        
            prev_logical_error_rate = statres.logical_error_rate;
            statres.logical_error_rate += (failure_rate*simulator.prob_k[k+min_k]);
            statres.n_logical_errors += fault_errors;

            if (world_rank == 0) {
                std::cout << "\tFailure rate = " << failure_rate << std::endl
                        << "\tNumber of errors = " << statres.n_logical_errors << std::endl;
                std::cout << "\tLogical error rate = " << statres.logical_error_rate << std::endl;
            }
        }
        if (world_rank == 0) {
            std::cout << "_______________________________________________"<< std::endl;
        }

        if(int(saved_hhw_syndromes[world_rank].size()) != 0){
            std::string file_name = syndrome_folder_name + syndrome_file_name("sgen",distance,physcial_error,
            meas_er, max_shot, world_rank, round);
            qpd::write_vector_of_vectors(saved_hhw_syndromes[world_rank], file_name);
            saved_hhw_syndromes[world_rank].clear();
            round++;
        }
    }

}

void bb_predecoder_ler_calc(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string decoder_name, std::string predecoder_name, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    uint round = 0;
    uint64_t local_clique_predecoded = 0;
    uint64_t clique_predecoded = 0;

    uint64_t local_hhw = 0;
    uint64_t hhw = 0;
    if(world_rank == 0){
        std::cout << "round#" <<round_n << " " << predecoder_name << "+" << decoder_name <<std::endl;
    }

    // for saving syndromes
    if (world_rank == 0 && save_syndromes){
        std::cout << syndrome_folder_name<< std::endl;
        // To check if we are reading correct files
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

    if(world_rank == 0){
        std::cout << "Circuit created, ";
    }
    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::benchmark::StatisticalResult statres;
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);
    qrc::Decoder* decoder;
    
    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    } 
    else if(decoder_name.find("ASTREAG") != std::string::npos){
        qrc::AstreaParams astreaG_param= {};
        astreaG_param.bfu_fetch_width = 2;
        astreaG_param.bfu_priority_queue_size = 8;
        astreaG_param.main_clock_frequency = 250e6;
        astreaG_param.bfu_compute_stages = 2;
        astreaG_param.n_registers = 2000;
        astreaG_param.use_mld = false;
        fp_t threshold_scale = 100;
        uint n_detectors_per_round = (distance*distance - 1)/2;
        uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(threshold_scale*0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    

        decoder =  new  qrc::Astrea(circ,       
                                    n_detectors_per_round,
                                    weight_filter_cutoff,
                                    astreaG_param);
    }

    std::function<qpd::PredecoderDecoderShotResult(const std::vector<uint8_t>&)> decode_func;

    if (predecoder_name == "CLIQUE") {
        decode_func = [predecoder](const std::vector<uint8_t>& synd) -> qpd::PredecoderDecoderShotResult { return predecoder->clique_decode_error(synd); };
    } else if (predecoder_name == "SMITH") {
        decode_func = [predecoder](const std::vector<uint8_t>& synd) -> qpd::PredecoderDecoderShotResult  { return predecoder->smith_decode_error(synd); };
    }
    
    

    predecoder->set_decoder(decoder);
    if(world_rank == 0){
        std::cout << "Simulator created, ";
    }
    std::mt19937_64 rng(world_size*round_n+world_rank);

    std::vector<vector<vector<uint16_t>>> saved_hhw_syndromes;
    saved_hhw_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits;

    simulator.set_shots_per_k(expected_ler, max_shot, true, hshots_replc, lshots_rplac);

    uint64_t max_k_ = max_k;
    fp_t sum_normalized_prob = 0;
    for(uint k = 0; k < max_k-min_k; k++){
        sum_normalized_prob += simulator.prob_k[k+min_k];
        if(simulator.shots_per_k[k+min_k] == 0){
            max_k_ --;
        }
    }
    max_k = max_k_;

    if(world_rank == 0){
        std::cout << "Simulator set the number of shots. " << std::endl;
    }
    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }
    fp_t prev_logical_error_rate = 0.0;
    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t local_errors = 0;
        local_clique_predecoded = 0;
        local_hhw = 0;
        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
   
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                auto res = decode_func(syndrome);
               
                
                local_errors += res.is_logical_error;
                // if(res.is_logical_error){
                //     flipped_bits = qpd::syndrome_compressed(syndrome);
                //     // saved_hhw_syndromes[world_rank].push_back(flipped_bits);

                // }
                if(res.pre_matching.size() != 0){
                    local_clique_predecoded++;
                }
                if(res.dec_matching.size()*2 > 10){
                    local_hhw++;
                }
                
            }
            shots -= shots_this_round;
        }
        uint64_t fault_errors = 0;
        MPI_Allreduce(&local_errors, &fault_errors, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t predecoded_n = 0;
        MPI_Allreduce(&local_clique_predecoded, &predecoded_n, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t hhw_n = 0;
        MPI_Allreduce(&local_hhw, &hhw_n, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }

        if (fault_errors > 0) {
            fp_t failure_rate = ((fp_t)fault_errors/simulator.shots_per_k[k+min_k]);
        
            prev_logical_error_rate = statres.logical_error_rate;
            statres.logical_error_rate += (failure_rate*simulator.prob_k[k+min_k]);
            statres.n_logical_errors += fault_errors;

            if (world_rank == 0) {
                std::cout << "\tFailure rate = " << failure_rate << std::endl
                        << "\tNumber of errors = " << statres.n_logical_errors << std::endl;
                std::cout << "\tLogical error rate = " << statres.logical_error_rate << std::endl;
            }
        }

        fp_t percentage_predecoded = ((fp_t)predecoded_n/(simulator.shots_per_k[k+min_k]));
        clique_predecoded +=percentage_predecoded*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        fp_t hhw_prob = ((fp_t)hhw_n/(simulator.shots_per_k[k+min_k]));
        hhw +=hhw_prob*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        if (world_rank == 0 && predecoder_name == "CLIQUE") {
            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\t\tThis bucket Clique%: " << percentage_predecoded*100 << std::endl;
            std::cout << "\t\tTotal Clique%: " << clique_predecoded*100 << std::endl;
            std::cout << "_______________________________________________"<< std::endl;
        }

        if (world_rank == 0 && predecoder_name == "SMITH") {
            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\t\tThis Bucket HHW: " << hhw_prob << std::endl;
            std::cout << "\t\tTotal HHW Prob(normalized): " << hhw << std::endl;
            std::cout << "_______________________________________________"<< std::endl;
        }

        // if(int(saved_hhw_syndromes[world_rank].size()) != 0){
        //     std::string file_name = syndrome_folder_name + syndrome_file_name("sgen",distance,physcial_error,
        //     meas_er, max_shot, world_rank, round);
        //     qpd::write_vector_of_vectors(saved_hhw_syndromes[world_rank], file_name);
        //     saved_hhw_syndromes[world_rank].clear();
        //     round++;
        // }
    }

}

void check_the_incorrect_syndrome(uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k,  std::ofstream& out,  bool generalized,
 std::string decoder_name, bool write_file){
    std::string syndrome_file =  "../NFDecoder/promatch_mistakes_d11_all.txt";//data/ChallengingSyndromes/Rtests/P11/promatch_mistakes_d11_all.txt";//All_PromatchAM_d13.txt";//A13/All_AstreaG_d13.txt" challengingSyndromes_d13/All_AstreaG_d13_final.txt";//Promatch_A_d13_all_Multi_Astrea.txt";////Promatch_A_d11_all_Multi_Astrea.txt";//final_promatchA_mistakes.txt";//ler_miscorrection_d13_ASTREAG_r1_b1416.txt";//All_AstreaG_d13.txt";//all_promatch_u1_d13.txt";//ler_miscorrection_d13_ASTREAG_r1_b1416.txt";//PromatchA_d13_r0_all.txt";//PromatchA_d13_r0.txt";;//P//promatchA_first_trial.txt";//promatchA_first_trial.txt";//;

    // std::cout << "Enter a the path to syndromes: ";
    // std::getline(std::cin, syndrome_file);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }

    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ, generalized);
    predecoder->set_decoder(decoder);

    qrc::AstreaParams astreaG_param= {};
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    astreaG_param.use_mld = false;


    vector<vector<uint16_t>> saved_syndromes_promatch_fails_astreag_decodes;
    vector<vector<uint16_t>> saved_syndromes_astreag_fails_promatch_decodes;
    vector<vector<uint16_t>> saved_syndromes_both_fails;

    vector<std::pair<fp_t,fp_t>> weight_promatch_fails_astreag_decodes;
    vector<std::pair<fp_t,fp_t>> weight_astreag_fails_promatch_decodes;
    vector<std::pair<fp_t,fp_t>> weight_both_fails;
    
    std::vector<uint64_t> promatch_accuracy(40,0);
    std::vector<uint64_t> promatch_skipped(40,0);
    std::vector<uint64_t> promatch_incorrect(40,0);

    std::vector<uint64_t> astreag_accuracy(40,0);
    std::vector<uint64_t> astreag_skipped(40,0);
    std::vector<uint64_t> astreag_incorrect(40,0);
    std::vector<fp_t> vector_sparcity_weighted;
    std::vector<fp_t> vector_sparcity_unweighted;
    fp_t sparcity_uw;
    

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(0.01*0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    

    qrc::Decoder*  astreaG=  new  qrc::Astrea(circ,
                                            n_detectors_per_round,
                                            weight_filter_cutoff,
                                            astreaG_param);

    std::vector<std::vector<uint16_t>> vecss = qpd::read_vector_of_vectors(syndrome_file);
    std::cout << "Size of incorrect pairs = " << vecss.size() << " - generalized = "<< predecoder->general<<std::endl;

    auto path_table = qrc::compute_path_table(predecoder->decoding_graph);

    uint promatch_le = 0;
    uint promatch_le_u1 = 0;
    uint promatch_le_u2 = 0;
    uint promatch_le_u0 = 0;
    uint promatch_le_promatch_a = 0;
    uint mwpm_le = 0;
    uint astreag_le = 0;
    uint count_icc = 0;
    uint smith_le = 0;

    uint low_pr_matching_total = 0;
    uint low_pr_matching_synd = 0;
    bool counted_lp = false;
    std::vector<uint> hw_cout(40,0);

    uint count_mistakes = 0;

    qpd::CountResult comp_res;
    for(auto compressed_synd : vecss){
        count_icc++;
        std::cout << "__________________________________________________****** Inc. Syndrome #" << count_icc << std::endl;
        std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);
        uint hw = std::accumulate(syndrome.begin(), syndrome.begin()+predecoder->n_detectors,0);
        std::cout << "HW = " << hw << std::endl;
        for(auto det: compressed_synd){
            if(det < predecoder->n_detectors){
                std::cout << det << " ";
            }
            if(det >= predecoder->n_detectors){
                std::cout << " - " << det;
            }
        }
        std::cout <<std::endl;
        fp_t sw = predecoder->graph_sparcity_calculator(compressed_synd, sparcity_uw, distance);
        vector_sparcity_weighted.push_back(sw);
        vector_sparcity_unweighted.push_back(sparcity_uw);

        auto res_u1 = predecoder->ensemble_decode_error(syndrome,1);
        auto res_u0 = predecoder->ensemble_decode_error(syndrome,0);
        auto res_promatch_a = predecoder->adaptively_decode_error(syndrome);
        
        auto res_mwpm = decoder->decode_error(syndrome);
        auto res_astreag = astreaG->decode_error(syndrome);

        auto res_smith = predecoder->smith_decode_error(syndrome);
        
        //qrc::DecoderShotResult res_astreag;
        mwpm_le += res_mwpm.is_logical_error;
        astreag_le += res_astreag.is_logical_error;
        promatch_le_u1 += res_u1.is_logical_error;
        promatch_le_u0 += res_u0.is_logical_error;
        promatch_le_promatch_a += res_promatch_a.is_logical_error;
        smith_le += res_smith.is_logical_error;

        std::cout << std::endl<< "_____________ MWPM Paths________________" << std::endl;
        predecoder->print_paths_matching(res_mwpm.matching);
        std::cout << std::endl<< "_____________________________" << std::endl;

        std::cout << std::endl<< "_____________ Promatch Paths________________" << std::endl;
        uint promatch_w = predecoder->print_paths_matching(res_promatch_a.matching);
        std::cout << std::endl<< "_____________________________" << std::endl;

        std::cout << std::endl<< "_____________ AstreaG Paths________________" << std::endl;
        uint astreag_w = predecoder->print_paths_matching(res_astreag.matching);
        std::cout << std::endl<< "_____________________________" << std::endl;

        if((promatch_w == astreag_w) && (!res_promatch_a.is_logical_error && res_astreag.is_logical_error)){
            std::cout << "equal edges" << std::endl;
        }


        //This part is for checking if group of edges 
        // that has low error rates are possible to be in 
        // matching.
        if(res_promatch_a.is_logical_error && res_astreag.is_logical_error){
            saved_syndromes_both_fails.push_back(compressed_synd);
            fp_t promatch_w = predecoder->calc_matching_weight(res_promatch_a.matching);
            fp_t astreag_w = predecoder->calc_matching_weight(res_astreag.matching);
            weight_both_fails.push_back(std::pair<fp_t,fp_t>(promatch_w,astreag_w));
        }
        if(res_promatch_a.is_logical_error && (!res_astreag.is_logical_error)){
            saved_syndromes_promatch_fails_astreag_decodes.push_back(compressed_synd);
            fp_t promatch_w = predecoder->calc_matching_weight(res_promatch_a.matching);
            fp_t astreag_w = predecoder->calc_matching_weight(res_astreag.matching);
            weight_promatch_fails_astreag_decodes.push_back(std::pair<fp_t,fp_t>(promatch_w,astreag_w));
            if(promatch_w < astreag_w){
                std::cout <<"\t\tERROR HERE | MWPM ER " << res_mwpm.is_logical_error << std::endl;
            }
        }
        if((!res_promatch_a.is_logical_error )&& (res_astreag.is_logical_error)){
            saved_syndromes_astreag_fails_promatch_decodes.push_back(compressed_synd);
            fp_t promatch_w = predecoder->calc_matching_weight(res_promatch_a.matching);
            fp_t astreag_w = predecoder->calc_matching_weight(res_astreag.matching);
            weight_astreag_fails_promatch_decodes.push_back(std::pair<fp_t,fp_t>(promatch_w,astreag_w));
            if(promatch_w > astreag_w){
                std::cout <<"\t\tERROR HERE | MWPM ER " << res_mwpm.is_logical_error << std::endl;
            }
        }
        counted_lp = false;
        std::set<uint> visited;
        for(const auto& m : res_mwpm.matching){
            if (visited.count(m.first) || visited.count(m.second)) {
            continue;
        }
            qrc::DecodingGraph::Vertex* v1 = predecoder->decoding_graph.get_vertex(m.first);
            qrc::DecodingGraph::Vertex* v2 = predecoder->decoding_graph.get_vertex(m.second);
            qrc::DecodingGraph::Edge* edge = predecoder->decoding_graph.get_edge(v1, v2);
            auto v_w = std::make_pair(v1, v2);
            // Finding error chains from the path table of decoding graph 
            uint ec = path_table[v_w].path.size() - 1;
            if(ec == 1 && pow(10,-1*path_table[v_w].distance) < 0.0001){
                std::cout << "FOUND A DISTANCE 1 OF LOW PROB (" << m.first<< "," << m.second<<")- " << pow(10,-1*path_table[v_w].distance)<< std::endl;
                if(!counted_lp){
                    low_pr_matching_synd++;
                }
                counted_lp = true;
                low_pr_matching_total++;
            }
             visited.insert(m.first);
            visited.insert(m.second);
            
        }
        
        if(res_promatch_a.is_logical_error && !res_mwpm.is_logical_error){
            std::cout << "INCORRECTLY DECODED _____________________________****** Inc. Syndrome #" << count_icc << std::endl; 

        }
        else if(!res_promatch_a.is_logical_error && res_mwpm.is_logical_error){
            std::cout << "EXTRA CORRECTLY DECODED _____________________________****** Inc. Syndrome #" << count_icc << std::endl; 

        }
        else if(res_promatch_a.is_logical_error && res_mwpm.is_logical_error){
            std::cout << "BOTH(M&P) INCORRECTLY DECODED _____________________________****** Inc. Syndrome #" << count_icc << std::endl; 
        }
        else{
            std::cout << "CORRECTLY DECODED _____________________________****** Inc. Syndrome #" << count_icc << std::endl; 
        }

        if(!res_promatch_a.is_logical_error && res_astreag.is_logical_error){
            std::cout << "  ASTREAG FAILED, PROMATCH CORRECTLY DECODED _____________________________****** Inc. Syndrome #" << count_icc << std::endl; 
        }
        std::cout << " ProMatch-A Weight = " << predecoder->calc_matching_weight(res_promatch_a.matching) << " MWPM Weight = " << predecoder->calc_matching_weight(res_mwpm.matching)
        << " AstreaG Weight = " << predecoder->calc_matching_weight(res_astreag.matching)
        << " Promatch-U0 Weight = " << predecoder->calc_matching_weight(res_u0.matching) 
        << " MWPM er " << res_mwpm.is_logical_error
        << " Astreag er " << res_astreag.is_logical_error
        << " Promatch-A er " << res_promatch_a.is_logical_error
        << " Promatch-U0 er " << res_u0.is_logical_error
        << " Smith er " << res_smith.is_logical_error<<std::endl;
            //qpd::print_syndrome(syndrome);
        // write_file = res_promatch_a.is_logical_error;
        for(auto det: compressed_synd){
            if(det < predecoder->n_detectors){
                if(write_file){
                    out << det << " ";
                }
            }
        }
        if(write_file){
            out << std::endl;
        }
        
        if (write_file){
            count_mistakes++;
        }
        std::cout << "___________MWPM MATCHING__________" << std::endl;
        // if(write_file) 
        //         out << "MWPM: "<< std::endl; 
        for(const auto& m: res_mwpm.matching){
            std::cout << "(" << m.first << ", " << m.second << ")-"; 
            if(write_file) 
                out << "(" << m.first << ", " << m.second << ")-"; 
        } 
        if(write_file) 
                out << std::endl;// << "Promatch: "<< std::endl;   
        std::cout << std::endl<< "_____________________________" << std::endl;
        std::cout << "___________Promatch-U0 MATCHING__________" << std::endl;
        for(const auto& m: res_u0.matching){
                std::cout << "(" << m.first << ", " << m.second << ")-";   
        }   
        std::cout << std::endl<< "_____________________________" << std::endl;
        std::cout << "___________Promatch-A MATCHING__________" << std::endl;
        for(const auto& m: res_promatch_a.matching){
                std::cout << "(" << m.first << ", " << m.second << ")-";
                if(write_file) 
                    out << "(" << m.first << ", " << m.second << ")-";  
        }  
        // if(write_file) 
        //         out << std::endl << "AstreaG: " << std::endl;  
        std::cout << std::endl<< "_____________________________" << std::endl;
        std::cout << "___________AstreaG MATCHING__________" << std::endl;
        for(const auto& m: res_astreag.matching){
                std::cout << "(" << m.first << ", " << m.second << ")-";
                // if(write_file)
                //     out << "(" << m.first << ", " << m.second << ")-";   
        }   
        // if(write_file) 
        //         out << std::endl; 
        std::cout << std::endl<< "_____________________________" << std::endl;
        std::cout<<"Adaptive-Predecoding pf_group" << std::endl;
        for(auto e: predecoder->adaptive_predecoding_map){
            std::cout << "(" << e.first << ", " << e.second << ")-"; 
            // if(write_file)  
            //     out << "(" << e.first << ", " << e.second << ")-";
        } 
        if(write_file)  
            out << std::endl << "________"<< std::endl;
        std::cout << std::endl<< "_____________MWPM Diff Promatch________________" << std::endl;
        auto comp_res = qpd::compare_matchings(res_mwpm.matching, res_promatch_a.matching);
        comp_res = qpd::accuracy_coverage(res_mwpm.matching, predecoder->adaptive_predecoding_map);
        promatch_accuracy[res_mwpm.matching.size()/2] += comp_res.countBothAB;
        promatch_incorrect[res_mwpm.matching.size()/2] += comp_res.countBNotInA;
        promatch_skipped[res_mwpm.matching.size()/2] += comp_res.countANotInB;
        std::cout << std::endl<< "_____________MWPM Diff AstreaG________________" << std::endl;
        comp_res = qpd::compare_matchings(res_mwpm.matching, res_astreag.matching);
        comp_res = qpd::accuracy_coverage(res_mwpm.matching, res_astreag.matching);
        astreag_accuracy[res_mwpm.matching.size()/2] += comp_res.countBothAB;
        astreag_incorrect[res_mwpm.matching.size()/2] += comp_res.countBNotInA;
        astreag_skipped[res_mwpm.matching.size()/2] += comp_res.countANotInB;
        hw_cout[res_mwpm.matching.size()/2]++;
        std::cout << std::endl<< "_____________Promatch Diff AstreaG________________" << std::endl;
        comp_res = qpd::compare_matchings(res_promatch_a.matching, res_astreag.matching);
    }

    std::cout << "mwpm #errors = " << mwpm_le << " astreaG #errors = " <<astreag_le << " Promatch-A #errors = " << promatch_le_promatch_a <<" Promatch-U1 #errors = " 
    << promatch_le_u1 << " Promatch-U0 #errors = " << promatch_le_u0 << " Smith #errors = " << smith_le << std::endl; ;
    std::cout << "Number of Synd wtih low prob matching: " << low_pr_matching_synd << " total of matchings with low prob: "<<low_pr_matching_total << std::endl;

    std::cout << "-------Cases that AstreaG fails but Promatch decodes (" << saved_syndromes_astreag_fails_promatch_decodes.size() << ")---" << std::endl;
    for(int i=0; i<saved_syndromes_astreag_fails_promatch_decodes.size(); i++ ){
        for(auto det: saved_syndromes_astreag_fails_promatch_decodes[i]){
            if(det < predecoder->n_detectors){
                std::cout << det << " ";
                // if(write_file){
                //     out << det << " ";
                // }
            }
            if(det >= predecoder->n_detectors){
                std::cout << " - " << det;
            }
        }

        std::cout << std::endl;
        std::cout << "\t(P: " << weight_astreag_fails_promatch_decodes[i].first <<", A: "
        << weight_astreag_fails_promatch_decodes[i].second << ")"<<std::endl;

        if(weight_astreag_fails_promatch_decodes[i].first  > weight_astreag_fails_promatch_decodes[i].second){
                std::cout <<"\t\t*ERROR HERE" << std::endl;
        }
    }
    std::cout << "-------Cases that Promatch fails but AstreaG decodes (" << saved_syndromes_promatch_fails_astreag_decodes.size() << ")---" << std::endl;
    for(int i=0; i<saved_syndromes_promatch_fails_astreag_decodes.size(); i++){
        for(auto det: saved_syndromes_promatch_fails_astreag_decodes[i]){
            if(det < predecoder->n_detectors){
                std::cout << det << " ";
                // if(write_file){
                //     out << det << " ";
                // }
            }
            if(det >= predecoder->n_detectors){
                std::cout << " - " << det;
            }
        }
        std::cout << std::endl;
        std::cout << "\t(P: " << weight_promatch_fails_astreag_decodes[i].first <<", A: "
        << weight_promatch_fails_astreag_decodes[i].second << ")"<<std::endl;

        if(weight_promatch_fails_astreag_decodes[i].first  < weight_promatch_fails_astreag_decodes[i].second){
                std::cout <<"\t\tERROR HERE" << std::endl;
        }
        
    }

    std::cout << "-------Cases that both fails (" << saved_syndromes_both_fails.size() << ")---" << std::endl;
    for(int i=0; i<saved_syndromes_both_fails.size(); i++){
        for(auto det: saved_syndromes_both_fails[i]){
            if(det < predecoder->n_detectors){
                std::cout << det << " ";
                // if(write_file){
                //     out << det << " ";
                // }
            }
            if(det >= predecoder->n_detectors){
                std::cout << " - " << det;
            }
        }
        std::cout << std::endl;
        std::cout << "\t(P: " << weight_both_fails[i].first <<", A: "
        << weight_both_fails[i].second <<")"<< std::endl;

        std::vector<uint8_t> mwpm_synd = qpd::syndrome_decompressed(saved_syndromes_both_fails[i], distance);

        std::cout <<"\tMWPM ER " << decoder->decode_error(mwpm_synd).is_logical_error << std::endl;
    }
    std::cout << "Promatch values --- accuracy" << std::endl;
    for(uint i = 0; i<40;i++){
        if(hw_cout[i]!=0)
            std::cout << i <<": " << ((double)promatch_accuracy[i]/(2*hw_cout[i])) << std::endl;
    }
    std::cout << "Promatch values --- skipped" << std::endl;
    for(uint i = 0; i<40;i++){
        if(hw_cout[i]!=0)
            std::cout << i <<": " << ((double)promatch_skipped[i]/(2*hw_cout[i])) << std::endl;
    }
    std::cout << "Promatch values --- incorrect" << std::endl;
    for(uint i = 0; i<40;i++){
        if(hw_cout[i]!=0)
            std::cout << i <<": " << ((double)promatch_incorrect[i]/(2*hw_cout[i])) << std::endl;
    }
    std::cout << "astreag values --- accuracy" << std::endl;
    for(uint i = 0; i<40;i++){
        if(hw_cout[i]!=0)
            std::cout << i <<": " << ((double)astreag_accuracy[i]/(2*hw_cout[i])) << std::endl;
    }
    std::cout << "astreag values --- skipped" << std::endl;
    for(uint i = 0; i<40;i++){
        if(hw_cout[i]!=0)
            std::cout << i <<": " << ((double)astreag_skipped[i]/(2*hw_cout[i]))<< std::endl;
    }
    std::cout << "astreag values --- incorrect" << std::endl;
    for(uint i = 0; i<40;i++){
        if(hw_cout[i]!=0)
            std::cout << i <<": " << ((double)astreag_incorrect[i]/(2*hw_cout[i])) << std::endl;
    }
    std::cout << count_mistakes << std::endl;
    std::cout << "----------- Unweighted Sparcity -----------" << std::endl;
    for(const auto& s: vector_sparcity_unweighted){
        std::cout << s <<", ";
    }
    double sum = std::accumulate(vector_sparcity_unweighted.begin(), vector_sparcity_unweighted.end(), 0.0);
    double mean = sum / vector_sparcity_unweighted.size();

    // Calculate variance
    double sq_sum = std::inner_product(vector_sparcity_unweighted.begin(), vector_sparcity_unweighted.end(), vector_sparcity_unweighted.begin(), 0.0);
    double variance = sq_sum / vector_sparcity_unweighted.size() - std::pow(mean, 2);

    // Calculate 2.5th and 97.5th percentiles for the 95% percentile interval
    std::sort(vector_sparcity_unweighted.begin(), vector_sparcity_unweighted.end());
    int lowerIndex = static_cast<int>(std::ceil(0.025 * vector_sparcity_unweighted.size()) - 1);
    int upperIndex = static_cast<int>(std::ceil(0.975 * vector_sparcity_unweighted.size()) - 1);
    double percentile2_5 = vector_sparcity_unweighted[lowerIndex];
    double percentile97_5 = vector_sparcity_unweighted[upperIndex];

    // Output results
    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Variance: " << variance << std::endl;
    std::cout << "95% Percentile Interval: [" << percentile2_5 << ", " << percentile97_5 << "]" << std::endl;
            
    std::cout<<"------------------------------" << std::endl;
    std::cout << "----------- Weighted Sparcity -----------" << std::endl;
    for(const auto& s: vector_sparcity_weighted){
        std::cout << s <<", ";
    }
    sum = std::accumulate(vector_sparcity_weighted.begin(), vector_sparcity_weighted.end(), 0.0);
    mean = sum / vector_sparcity_weighted.size();

    // Calculate variance
    sq_sum = std::inner_product(vector_sparcity_weighted.begin(), vector_sparcity_weighted.end(), vector_sparcity_weighted.begin(), 0.0);
    variance = sq_sum / vector_sparcity_weighted.size() - std::pow(mean, 2);

    // Calculate 2.5th and 97.5th percentiles for the 95% percentile interval
    std::sort(vector_sparcity_weighted.begin(), vector_sparcity_weighted.end());
    lowerIndex = static_cast<int>(std::ceil(0.025 * vector_sparcity_weighted.size()) - 1);
    upperIndex = static_cast<int>(std::ceil(0.975 * vector_sparcity_weighted.size()) - 1);
    percentile2_5 = vector_sparcity_weighted[lowerIndex];
    percentile97_5 = vector_sparcity_weighted[upperIndex];

    // Output results
    std::cout << "Mean: " << mean << std::endl;
    std::cout << "Variance: " << variance << std::endl;
    std::cout << "95% Percentile Interval: [" << percentile2_5 << ", " << percentile97_5 << "]" << std::endl;
            



}

void stats_of_incorrect_syndromes(uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool generalized,
 std::string decoder_name){
    std::string syndrome_file = "../NFDecoder/data/challengingSyndromes/RTests/P13/All_PromatchAM_d13";
    // std::cout << "Enter a the path to syndromes: ";
    // std::getline(std::cin, syndrome_file);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }

    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ, generalized);
    predecoder->set_decoder(decoder);


    qrc::AstreaParams astreaG_param= {};
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    astreaG_param.use_mld = false;

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    

    qrc::Decoder*  astreaG=  new  qrc::Astrea(circ,
                                            n_detectors_per_round,
                                            weight_filter_cutoff,
                                            astreaG_param);

    std::vector<std::vector<uint16_t>> vecss = qpd::read_vector_of_vectors(syndrome_file);
    uint total_n = ((uint)vecss.size());
    std::cout << "Size of incorrect pairs = " << total_n << std::endl;
    

    uint promatch_le_u1 = 0;
    uint promatch_le_u0 = 0;
    uint mwpm_le = 0;
    uint astreag_le = 0;
    uint count_icc = 0;

    auto path_table = qrc::compute_path_table(predecoder->decoding_graph);

    uint64_t times_AR_is_1 = 0;
    uint64_t times_CR_is_1 = 0;
    uint64_t times_AR2_is_1 = 0;
    uint64_t times_CR2_is_1 = 0;
    fp_t AR_avg = 0;
    fp_t CR_avg = 0;
    fp_t AR2_avg = 0;
    fp_t CR2_avg = 0;
    uint64_t PG_high = 0;
    uint64_t PG_10 = 0;
    uint64_t PG_2 = 0;
    uint64_t PG_0 = 0;
    uint64_t PG_4 = 0;
    uint64_t PG_6 = 0;
    uint64_t PG_8 = 0;
    uint64_t times_1_incorrect_pair = 0;
    uint64_t times_2_incorrect_pair = 0;
    uint64_t times_3_incorrect_pair = 0;
    uint64_t times_4_incorrect_pair = 0;
    uint64_t more_than_4_incorrect_pair = 0;

    uint64_t hw_r_high = 0;
    uint64_t hw_r_10 = 0;
    uint64_t hw_r_2 = 0;
    uint64_t hw_r_0 = 0;
    uint64_t hw_r_4 = 0;
    uint64_t hw_r_6 = 0;
    uint64_t hw_r_8 = 0;

    // for three incorrect pairs
    uint64_t inc_pair_3_hw_6 = 0;
    uint64_t inc_pair_3_hw_8 = 0;
    uint64_t inc_pair_3_hw_10 = 0;
    uint64_t inc_pair_3_hhw = 0;
    uint64_t inc_pair_2_hw_8 = 0;
    uint64_t inc_pair_2_hw_10 = 0;
    uint64_t inc_pair_2_hhw = 0;
    uint64_t inc_pair_1_hw_10 = 0;
    uint64_t inc_pair_1_hhw = 0;
    std::array<uint64_t, 10> ec_all_hw;
    ec_all_hw.fill(0);

    std::array<uint64_t, 40> hw_numbers;
    hw_numbers.fill(0);

    std::array<fp_t, 10> ec_ratio;
    ec_ratio.fill(0);

    std::array<fp_t, 10> local_ec;
    local_ec.fill(0);
    
    uint64_t sum_hw = 0;
    
    for(auto compressed_synd : vecss){
        count_icc++;
        //std::cout << "__________________________________________________****** Inc. Syndrome #" << count_icc;
        std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);
        auto res_u1 = predecoder->ensemble_decode_error(syndrome,1);
        // auto res_u2 = predecoder->ensemble_decode_error(syndrome,2);
        auto res_u0 = predecoder->ensemble_decode_error(syndrome,0);
        
        auto res2 = decoder->decode_error(syndrome);
        auto res3 = astreaG->decode_error(syndrome);
        mwpm_le += res2.is_logical_error;
        astreag_le += res3.is_logical_error;
        promatch_le_u1 += res_u1.is_logical_error;
        promatch_le_u0 += res_u0.is_logical_error;

        if(res_u0.is_logical_error and !res2.is_logical_error){
            std::cout << "Incorrect U0 s" << count_icc<< endl;
            for(auto xx :  compressed_synd){
                std::cout << xx << " ";
            }
            std::cout <<std::endl;
        }



        std::vector<qrc::DecodingGraph::Vertex*> A; // Actual Fast Group
        std::vector<qrc::DecodingGraph::Edge*> A_edge; // Actual Fast Group as edges
        std::vector<qrc::DecodingGraph::Vertex*> B; // Predicted Fast Group
        
        uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
        sum_hw += hw;
        hw_numbers[hw] ++;
                //Creating B group
                predecoder->set_the_priorities(syndrome);

                auto decoder_results = decoder->decode_error(syndrome);
                uint mn = 0;
                local_ec.fill(0);
                for(const auto& matching: decoder_results.matching){
                    mn++;
                    qrc::DecodingGraph::Vertex* v1 = predecoder->decoding_graph.get_vertex(matching.first);
                    qrc::DecodingGraph::Vertex* v2 = predecoder->decoding_graph.get_vertex(matching.second);
                    qrc::DecodingGraph::Edge* edge = predecoder->decoding_graph.get_edge(v1, v2);
                    if(edge != nullptr){
                        //Creating A group
                        qpd::addUniqueVertex(A, v1);
                        qpd::addUniqueVertex(A, v2);
                        qpd::addUniqueEdge(A_edge, edge);
                        
                    }
                    auto v_w = std::make_pair(v1, v2);
                    // Finding error chains from the path table of decoding graph 
                    uint ec = path_table[v_w].path.size() - 1;
                    if(count_icc == 193){
                        std::cout << ec << " ";
                    }
                    ec_all_hw[ec]++;
                    local_ec[ec]++;   
                }
                if(count_icc == 193){
                        std::cout << endl;
                    }
                for(uint l=0;l<10;l++){
                    ec_ratio[l] += ((double)local_ec[l]/mn);
                }
                if( predecoder->number_of_edges < hw){
                    if( 10<(hw-predecoder->number_of_edges)){
                        hw_r_high++;
                    }
                    else if( 10==(hw-predecoder->number_of_edges) || 9==(hw-predecoder->number_of_edges)){
                        hw_r_10++;
                    }
                    else if( 8==(hw-predecoder->number_of_edges) || 7==(hw-predecoder->number_of_edges)){
                        hw_r_8++;
                    }
                    else if( 6==(hw-predecoder->number_of_edges) || 5==(hw-predecoder->number_of_edges)){
                        hw_r_6++;
                    }
                    else if( 4==(hw-predecoder->number_of_edges) || 3==(hw-predecoder->number_of_edges)){
                        hw_r_4++;
                    }
                    else if( 2==(hw-predecoder->number_of_edges) || 1==(hw-predecoder->number_of_edges)){
                        hw_r_2++;
                    }
                }
                else if(predecoder->number_of_edges == hw){
                    hw_r_0++;
                }
                //std::cout << world_rank << ":&&&" <<hw_r_0[k]<<std::endl;
                
                uint p=0;
                //for(uint p=0;p<6;p++){
                
                    B = predecoder->pf_groups[p];
                    qpd::CountResult c2 = qpd::countMembersInAB(A,B);
                    qpd::CountResult c = {0,0,0};
                    c = qpd::countMembersInAB_edges(A_edge,predecoder->pf_pairs[p]);
                    
                    fp_t CR2 = A.size() !=0 ? ((double)c2.countBothAB/(c2.countBothAB+c2.countANotInB)):0;

                    fp_t CR = A_edge.size() !=0 ? ((double)c.countBothAB/(c.countBothAB+c.countANotInB)):0;
                    // if((c.countBothAB + c.countANotInB) != A.size()){
                    //     std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. A - B + AandB != A";
                    // }
                    

                    fp_t AR2 = B.size()!=0 ?((double)c2.countBothAB/(c2.countBothAB+c2.countBNotInA)):0;
                    fp_t AR = predecoder->pf_pairs[p].size()!=0 ?((double)c.countBothAB/(c.countBothAB+c.countBNotInA)):0;

                    
                    
                    // if((c.countBothAB + c.countBNotInA) != B.size()){
                    //     std::cout<< "THIS MESSAGE SHOULD NOT BE PRINTED. B - A + AandB != B";
                    // }
                    
                    if(AR == 1){
                        times_AR_is_1++;
                    }
                    AR_avg += AR;
                    if( CR == 1){
                        times_CR_is_1++;
                    }
                    CR_avg += CR;

                    if(AR2 == 1){
                        times_AR2_is_1++;
                    }
                    AR2_avg += AR2;
                    if( CR2 == 1){
                        times_CR2_is_1++;
                    }
                    CR2_avg += CR2;
                    

                    if( 10<B.size()){
                        PG_high++;
                    }
                    else if( 10==B.size() || 9==B.size()){
                        PG_10++;
                    }
                    else if( 8==B.size()|| 7==B.size()){
                        PG_8++;
                    }
                    else if( 6==B.size() || 5==B.size()){
                        PG_6++;
                    }
                    else if( 4==B.size()|| 3==B.size()){
                        PG_4++;
                    }
                    else if( 2==B.size()|| 1==B.size()){
                        PG_2++;
                    }

                    uint predecoding_size = ((uint)predecoder->pf_pairs[p].size())*2;

                    if(1 == c.countBNotInA){
                        times_1_incorrect_pair++;
                        if((hw-predecoding_size)==9 || (hw-predecoding_size)==10){
                            inc_pair_1_hw_10++;
                        }
                        else if(10<(hw-predecoding_size)){
                            inc_pair_1_hhw++;
                        }
                    }
                    else if(2 == c.countBNotInA){
                        times_2_incorrect_pair++;
                        if((hw-predecoding_size)==7 || (hw-predecoding_size)==8){
                            inc_pair_2_hw_8++;
                        }
                        else if((hw-predecoding_size)==9 || (hw-predecoding_size)==10){
                            inc_pair_2_hw_10++;
                        }
                        else if(10<(hw-predecoding_size)){
                            inc_pair_2_hhw++;
                        }
                    }
                    else if(3 == c.countBNotInA){
                        times_3_incorrect_pair++;
                        if((hw-predecoding_size)==5 || (hw-predecoding_size)==6){
                            inc_pair_3_hw_6++;
                        }
                        else if((hw-predecoding_size)==7 || (hw-predecoding_size)==8){
                            inc_pair_3_hw_8++;
                        }
                        else if((hw-predecoding_size)==9 || (hw-predecoding_size)==10){
                            inc_pair_3_hw_10++;
                        }
                        else if(10<(hw-predecoding_size)){
                            inc_pair_3_hhw++;
                        }

                    }
                    else if(4 == c.countBNotInA){
                        times_4_incorrect_pair++;
                    }
                    else if(4 < c.countBNotInA){
                        more_than_4_incorrect_pair++;
                    }
                //}

        
    }

  std::cout << "________________________________________________________________" <<std::endl;
        std::cout <<  std::setw(16) << " N(HW-|B| = 0) | "; 
        std::cout <<  std::setw(16) << " N(HW-|B| = 2) | "; 
        std::cout <<  std::setw(16) << " N(HW-|B| = 4) | "; 
        std::cout <<  std::setw(16) << " N(HW-|B| = 6) | "; 
        std::cout <<  std::setw(16) << " N(HW-|B| = 8) | "; 
        std::cout <<  std::setw(16) << " N(HW-|B| = 10) | "; 
        std::cout <<  std::setw(16) << " N(10 < HW-|B|) | " << std::endl;
        std::cout << std::setw(10) <<hw_r_0 << " | ";
        std::cout << std::setw(10) <<hw_r_2  <<" | ";
        std::cout << std::setw(10) <<hw_r_4 <<" | ";
        std::cout << std::setw(10) <<hw_r_6 <<" | ";
        std::cout << std::setw(10) << hw_r_8 <<" | ";
        std::cout << std::setw(10) <<hw_r_10<<" | ";
        std::cout << std::setw(10) << hw_r_high <<" | "
        << std::endl;
            
        

        
        
            std::cout << "######################PRIORITY GROUP " << 0 << "######################" << std::endl;
            std::cout << "Distance: "<< distance << " - #Buckets Max: " << max_k << " Min: "<< min_k <<" - Physical ER: " << physcial_error << std::endl;
            std::cout << "mean hw = " << ((double)sum_hw/total_n);
            std::cout << "________________________________________________________________" <<std::endl;

            std::cout << " P(AR = 1) | ";
            std::cout << " P(CR = 1) | ";
            std::cout << " P(AR2 = 1) | ";
            std::cout << " P(CR2 = 1) | "<< std::endl;
 
            std::cout <<times_AR_is_1 <<" | ";
            std::cout <<times_CR_is_1  <<" | ";
            std::cout <<times_AR2_is_1  <<" | ";
            std::cout <<times_CR2_is_1 <<" | "<< std::endl;

            std::cout << "________________________________________________________________" <<std::endl;
            std::cout  << " AVG AR | ";
            std::cout << " AVG CR | ";
            std::cout << " AVG AR2 | ";
            std::cout  << " AVG CR2 | "<< std::endl;
            
           
            std::cout  <<AR_avg<<" | ";
            std::cout  <<CR_avg <<" | ";
            std::cout  <<AR2_avg <<" | ";
            std::cout  <<CR2_avg <<" | "<< std::endl;

            std::cout  <<((double)AR_avg/total_n)<<" | ";
            std::cout  <<((double)CR_avg/total_n) <<" | ";
            std::cout  <<((double)AR2_avg/total_n) <<" | ";
            std::cout  <<((double)CR2_avg /total_n)<<" | "<< std::endl;

            std::cout << "________________________________________________________________" <<std::endl;            
            std::cout << " N(|PG| = 2) | "; 
            std::cout  << " N(|PG| = 4) | "; 
            std::cout  << " N(|PG| = 6) | "; 
            std::cout << " N(|PG| = 8) | "; 
            std::cout << " N(|PG| = 10) | "; 
            std::cout  << " N(10 < |PG|) | " << std::endl;
            
            std::cout <<PG_2 << " | ";
            std::cout <<PG_4 << " | ";
            std::cout <<PG_6 << " | ";
            std::cout <<PG_8 << " | ";
            std::cout <<PG_10 << " | ";
            std::cout  <<PG_high << " | "
            << std::endl;

            std::cout << "________________________________________________________________" <<std::endl;
            std::cout << " N(|InCorr. Pairs|=1) | "; 
            std::cout << " N(|InCorr. Pairs|=2) |";
            std::cout << " N(|InCorr. Pairs|=3) |";
            std::cout << " N(|InCorr. Pairs|=4) |";
            std::cout << " N(4<|InCorr. Pairs|)|" << std::endl;

        
            std::cout <<times_1_incorrect_pair <<" | ";
            std::cout <<times_2_incorrect_pair <<" | ";
            std::cout <<times_3_incorrect_pair <<" | ";
            std::cout <<times_4_incorrect_pair <<" | ";
            std::cout <<more_than_4_incorrect_pair <<" | "
            << std::endl;
                

            std::cout << "________________________________________________________________" <<std::endl;
            std::cout << " N(IP=1,RHW=10) | "; 
            std::cout << " N(IP=1,10<RHW) | ";
            std::cout << " N(IP=2,RHW=8) | ";
            std::cout << " N(IP=2,RHW=10) | ";
            std::cout << " N(IP=2,10<RHW) | ";
            std::cout << " N(IP=3,RHW=6) | ";
            std::cout << " N(IP=3,RHW=8) | ";
            std::cout << " N(IP=3,RHW=10) | ";
            std::cout << " N(IP=3,10<RHW) | " << std::endl;

            
            std::cout <<inc_pair_1_hw_10 <<" | ";
            std::cout <<inc_pair_1_hhw <<" | ";
            std::cout <<inc_pair_2_hw_8 <<" | ";
            std::cout <<inc_pair_2_hw_10 <<" | ";
            std::cout <<inc_pair_2_hhw <<" | ";
            std::cout <<inc_pair_3_hw_6 <<" | ";
            std::cout <<inc_pair_3_hw_8 <<" | ";
            std::cout <<inc_pair_3_hw_10 <<" | ";
            std::cout <<inc_pair_3_hhw <<" | "
            << std::endl;

            std::cout << "_____________________________ec_all_hw___________________________________" <<std::endl;
            for(uint i = 0; i < 10; i++){
                if (ec_all_hw[i]){
                    std::cout << i << "," << ec_all_hw[i] << "\n";
                }
            }
            std::cout << "________________________________ec_ratio________________________________" <<std::endl;
            for(uint i = 0; i < 10; i++){
                if (ec_ratio[i]){
                    std::cout << i << "," << ec_ratio[i] << " - " << ((double)ec_ratio[i]/total_n) << std::endl;
                }
            }
            std::cout << "____________________________HW____________________________________" <<std::endl;
            for(uint i = 0; i < 10; i++){
                if (ec_ratio[i]){
                    std::cout << i << "," << hw_numbers[i] << " - " << ((double)hw_numbers[i]/total_n) << std::endl;
                }
            }

            std::cout << "MWPM #errors = " << mwpm_le << " AstreaG #errors = " <<astreag_le <<" PROM-U1 #errors = " << promatch_le_u1<<" PROM-U0 #errors = " << promatch_le_u0;

}




void bb_avg_ec1(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // for saving syndromes

    std::vector<fp_t> avg_ec1(max_k - min_k, 0);
    std::vector<fp_t> avg_ec1_acc(max_k - min_k, 0);

    std::vector<uint64_t> n_over_bar(max_k - min_k, 0);
    std::vector<uint64_t> n_over_bar_acc(max_k - min_k, 0);

    std::vector<std::vector<fp_t>> ec_avg(max_k - min_k, std::vector<fp_t>(50, 0));
    std::vector<std::vector<fp_t>> ec_avg_acc(max_k - min_k, std::vector<fp_t>(50, 0));

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::benchmark::StatisticalResult statres;
    qrc::Decoder* decoder;
    

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    auto path_table = qrc::compute_path_table(simulator.decoding_graph);

    std::mt19937_64 rng(world_rank);

    std::vector<vector<vector<uint16_t>>> saved_hhw_syndromes;
    saved_hhw_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits;

    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true);
    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }
    fp_t prev_logical_error_rate = 0.0;
    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t local_errors = 0;
        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                // std::cout << "!!1" << std::endl;
                auto res = decoder->decode_error(syndrome);
                uint n_ec1 = 0;
                for(auto m : res.matching){
                    auto v = simulator.decoding_graph.get_vertex(m.first);
                    auto w = simulator.decoding_graph.get_vertex(m.second);
                    auto v_w = std::make_pair(v, w);
                    uint ec = path_table[v_w].path.size() -1;
                    if(ec == 1){
                        n_ec1 += 1;
                    }
                    ec_avg[k][ec] += ((double)1/res.matching.size());;
                }
                if(10<(res.matching.size()-n_ec1)){
                    n_over_bar[k]++;
                }
                avg_ec1[k] += ((double)n_ec1/res.matching.size());    
            }
            shots -= shots_this_round;
        }
    }

    MPI_Reduce(avg_ec1.data(), avg_ec1_acc.data(), avg_ec1.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(n_over_bar.data(), n_over_bar_acc.data(), n_over_bar.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);

    std::vector<fp_t> flattened_ec_avg;

    for (const auto& row : ec_avg) {
        flattened_ec_avg.insert(flattened_ec_avg.end(), row.begin(), row.end());
    }

    std::vector<fp_t> flattened_ec_avg_acc(flattened_ec_avg.size());

    MPI_Reduce(flattened_ec_avg.data(), flattened_ec_avg_acc.data(), flattened_ec_avg.size(),
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);


    for (int i = 0; i < ec_avg_acc.size(); i++) {
        ec_avg_acc[i] = std::vector<fp_t>(flattened_ec_avg_acc.begin() + i * 50,
                                                 flattened_ec_avg_acc.begin() + (i + 1) * 50);
    }

    if (world_rank == 0) {
    std::cout << "________________________________________________________________" <<std::endl;
    std::cout << std::setw(18) <<  " k (p) | ";  
    std::cout <<  std::setw(16) << " Avg of Percentage of #(|EC|=1) |";
    std::cout <<  std::setw(16) << " Prob. Over 10 remains|"
        << std::endl;
        for ( uint i = 0; i< max_k-min_k;  i++){
            std::string printing_shots = std::to_string(simulator.shots_per_k[i+min_k]);
            if(1'000<=simulator.shots_per_k[i+min_k]){
                printing_shots = number_of_shots_to_string(simulator.shots_per_k[i+min_k]);
            }
            std::cout << i+min_k << " (" << simulator.prob_k[i+min_k]<< ") | "; 
            std::cout << std::setw(10) <<((double)avg_ec1_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << avg_ec1_acc[i] << "/" << printing_shots <<") | ";
            std::cout << std::setw(10) <<((double)n_over_bar_acc[i]/simulator.shots_per_k[i+min_k]) << "(" << n_over_bar_acc[i] << "/" << printing_shots <<") | "

            << std::endl;
            
        } 

        std::cout << "________________________________________________________________" <<std::endl;

        std::vector<fp_t> ec_freq(50,0);
        std::vector<fp_t> ec_freq_no_bucket(50,0);
        fp_t normalize_probs = 0;
        for ( uint i = 0; i< max_k-min_k;  i++){
            normalize_probs += simulator.prob_k[i+min_k];
        }
        for(uint j=0;j < 50;j++){
            for ( uint i = 0; i< max_k-min_k;  i++){
                fp_t fr = ec_avg_acc[i][j];
                if(fr){
                    ec_freq[j] += fr*((double)1/simulator.shots_per_k[i+min_k])*simulator.prob_k[i+min_k];  
                    ec_freq_no_bucket[j] += fr*((double)1/simulator.shots_per_k[i+min_k])*(simulator.prob_k[i+min_k]/normalize_probs);  
                }
            }
        }
        std::cout << "________________________________________________________________" <<std::endl;

        for(uint j=0;j < 50;j++){
            if(ec_freq[j]){
                std::cout << j << ": " << ec_freq[j] << std::endl;
                    
            }
        }
        std::cout << "________________no bucket________________________________________________" <<std::endl;

        for(uint j=0;j < 50;j++){
            if(ec_freq_no_bucket[j]){
                std::cout << j << ": " << ec_freq_no_bucket[j] << std::endl;
                    
            }
        }
    }  
}

void sample_test_degree_priority(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, 
    bool& print_time, std::string decoder_name){


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


    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }

    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);
    predecoder->set_decoder(decoder);

    qpd::Predecoder* predecoder_gen = new qpd::Predecoder(circ, true);
    predecoder_gen->set_decoder(decoder);
    
    std::cout << "Enter the Syndrome: " << std::endl;
    std::string numbers;
    std::getline(std::cin, numbers);

    std::vector<uint16_t> result;
    std::istringstream iss(numbers);
    uint16_t num;

    while (iss >> num) {
        result.push_back(num);
    }

    qrc::AstreaParams astreaG_param= {};
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    astreaG_param.use_mld = false;

    

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(0.01*0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    

    qrc::Decoder*  astreaG=  new  qrc::Astrea(circ,
                                            n_detectors_per_round,
                                            weight_filter_cutoff,
                                            astreaG_param);



    std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(result, distance);
    std::cout << std::endl << "HW = " <<std::accumulate(syndrome.begin(), syndrome.begin()+predecoder->n_detectors,0) << std::endl;

    auto res = predecoder->adaptively_decode_error(syndrome);
    auto res_promatch_g = predecoder_gen->adaptively_decode_error(syndrome);
    auto res_mwpm = decoder->decode_error(syndrome);
    auto res_astreag = astreaG->decode_error(syndrome);
    auto res_smith = predecoder->smith_decode_error(syndrome);
    auto res_clique = predecoder->clique_decode_error(syndrome);

    std::cout << "___________MWPM MATCHING__________" << std::endl;
    for(const auto& m: res_mwpm.matching){
        std::cout << "(" << m.first << ", " << m.second << ")-";   
    }   
    std::cout << std::endl<< "_____________________________" << std::endl;
 
    std::cout << "___________Promatch MATCHING__________" << std::endl;
    for(const auto& m: res.matching){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    }   
    std::cout << std::endl<< "_____________________________" << std::endl;
    std::cout << "___________Promatch Predecoding MATCHING__________" << std::endl;
    for(const auto& m: predecoder->adaptive_predecoding_map){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    }   
    std::cout << std::endl<< "_____________________________" << std::endl;
    std::cout << "___________G-Promatch MATCHING__________" << std::endl;
    for(const auto& m: res_promatch_g.matching){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    }   
    std::cout << std::endl<< "_____________________________" << std::endl;
    std::cout << "___________G-Promatch Predecoding MATCHING__________" << std::endl;
    for(const auto& m: predecoder_gen->adaptive_predecoding_map){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    } 
    std::cout << std::endl << "___________AstreaG MATCHING__________" << std::endl;
    for(const auto& m: res_astreag.matching){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    } 
    std::cout << std::endl << "___________Smith MATCHING__________" << std::endl;
    for(const auto& m: res_smith.pre_matching){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    } 
    std::cout << std::endl << "\t -- Smith Decoder__________" << std::endl;
    for(const auto& m: res_smith.dec_matching){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    } 
    std::cout << std::endl << "\t -- Smith Solution Edges__________" << std::endl;
    for(const auto& m: res_smith.correcting_edge){
            std::cout << "(" << m.first.first << ", " << m.first.second << "): " <<
            m.second << " - ";    
    } 
    std::cout << std::endl << "___________Clique MATCHING__________" << std::endl;
    for(const auto& m: res_clique.pre_matching){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    } 
    std::cout << std::endl << "\t -- Clique Decoder__________" << std::endl;
    for(const auto& m: res_clique.dec_matching){
            std::cout << "(" << m.first << ", " << m.second << ")-";   
    } 
    std::cout << std::endl << "\t -- Clique Solution Edges__________" << std::endl;
    for(const auto& m: res_clique.correcting_edge){
            std::cout << "(" << m.first.first << ", " << m.first.second << "): " <<
            m.second << " - ";    
    }

    std::cout << std::endl<< "_____________________________" << std::endl;
    std::cout << std::endl<< "_____________MWPM Diff Promatch________________" << std::endl;
    auto comp_res = qpd::compare_matchings(res_mwpm.matching, res.matching);
    std::cout << std::endl<< "_____________________________" << std::endl;
    std::cout << std::endl<< "_____________MWPM Diff G-Promatch________________" << std::endl;
    comp_res = qpd::compare_matchings(res_mwpm.matching, res_promatch_g.matching);
    std::cout << std::endl<< "_____________________________" << std::endl;
    std::cout << std::endl<< "_____________MWPM Diff AstreaG" << std::endl;
    comp_res = qpd::compare_matchings(res_mwpm.matching, res_astreag.matching);
    std::cout << std::endl<< "_____________________________" << std::endl;
    
    std::cout << std::endl<< "_____________ MWPM Paths________________" << std::endl;
    predecoder->print_paths_matching(res_mwpm.matching);
    std::cout << std::endl<< "_____________________________" << std::endl;

    std::cout << std::endl<< "_____________ Promatch Paths________________" << std::endl;
    predecoder->print_paths_matching(res.matching);
    std::cout << std::endl<< "_____________________________" << std::endl;

    std::cout << std::endl<< "_____________ AstreaG Paths________________" << std::endl;
    predecoder->print_paths_matching(res_astreag.matching);
    std::cout << std::endl<< "_____________________________" << std::endl;


                    


    std::cout << "MWPM Weight = " << predecoder->calc_matching_weight(res_mwpm.matching) << " LE " << res_mwpm.is_logical_error << std::endl;
    std::cout << "PROMATCH Weight = " << predecoder->calc_matching_weight(res.matching) << " LE " << res.is_logical_error << std::endl;
    std::cout << "G-PROMATCH Weight = " << predecoder->calc_matching_weight(res_promatch_g.matching) << " LE " << res_promatch_g.is_logical_error << std::endl;
    std::cout << "AstreaG Weight = " << predecoder->calc_matching_weight(res_astreag.matching) << " LE " << res_astreag.is_logical_error << std::endl;
    std::cout << "Smith Weight = " << res_smith.weight << " LE " << res_smith.is_logical_error<< std::endl;
    std::cout << "Clique Weight = " << res_clique.weight << " LE " << res_clique.is_logical_error<< std::endl;

    // auto all_matchings = qpd::getAllCompleteMatchings(result);
    // uint all_matchings_number = 0;
    // for (const auto& matches :all_matchings){
    //     all_matchings_number++;
    //     std::cout << all_matchings_number << ": ";
    //     for(const auto& m: matches){
    //         std::cout << "(" << m.first << ", " << m.second << ")-";   
    //     } 
    // }

    std::cout << "DONE!";
}

void degree_distribution(std::ofstream& out, uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::array<fp_t,100> degree_distr;
    degree_distr.fill(0);
    std::array<uint64_t,100> degree_sum;
    degree_sum.fill(0);
    std::array<uint64_t,100> degree_sum_ac;
    degree_sum_ac.fill(0);

    std::array<fp_t,100> degree_distr_2;
    degree_distr_2.fill(0);
    std::array<uint64_t,100> degree_sum_2;
    degree_sum_2.fill(0);
    std::array<uint64_t,100> degree_sum_2_ac;
    degree_sum_2_ac.fill(0);

    std::array<uint64_t,100> single_syndrome;
    single_syndrome.fill(0);


    std::array<fp_t,100> complex_nodes_distr;
    complex_nodes_distr.fill(0);
    std::array<uint64_t,100> complex_nodes_sum;
    complex_nodes_sum.fill(0);
    std::array<uint64_t,100> complex_nodes_sum_ac;
    complex_nodes_sum_ac.fill(0);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }
    else if(distance == 15){
        max_k = 21;
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
  
    std::mt19937_64 rng(world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);


    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, hshots_replc, lshots_rplac);
    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }

    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t local_errors = 0;
        degree_sum.fill(0);
        complex_nodes_sum.fill(0);
        degree_sum_2.fill(0);
        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);

            
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                /////////
                // std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);

                predecoder->update_adjacency_matrix_and_detector_list(syndrome, 1, true);
                auto g = predecoder->get_adjacency_matrix();

                single_syndrome.fill(0);

                for(const auto& det : g.detector_list){
                    auto v = predecoder->decoding_graph.get_vertex(det);
                    if(g.adjacency_matrix.count(v) == 0 && single_syndrome[0] == 0){
                        degree_sum[0]++;
                        single_syndrome[0] == 1;
                    }
                    else{
                        uint degree = g.adjacency_matrix[v].size();
                        if(single_syndrome[degree] == 0){
                            degree_sum[degree]++;
                            single_syndrome[degree] == 1;

                        }
                        
                    }
                }

                uint n_complex_nodes = 0;
                g = predecoder->get_adjacency_matrix();
                for(const auto& det : g.detector_list){
                    auto v = predecoder->decoding_graph.get_vertex(det);
                    if(g.adjacency_matrix.count(v) == 0){
                        n_complex_nodes++;
                    }
                    else{
                        uint degree = g.adjacency_matrix[v].size();
                        if(degree > 1){
                            n_complex_nodes++;
                        }
                    }
                }

                complex_nodes_sum[n_complex_nodes]++;

                predecoder->update_adjacency_matrix_and_detector_list(syndrome, 1, false);
                g = predecoder->get_adjacency_matrix();

                single_syndrome.fill(0);

                for(const auto& det : g.detector_list){
                    auto v = predecoder->decoding_graph.get_vertex(det);
                    if(g.adjacency_matrix.count(v) == 0 && single_syndrome[0] == 0){
                        degree_sum_2[0]++;
                        single_syndrome[0] == 1;
                    }
                    else{
                        uint degree = g.adjacency_matrix[v].size();
                        if(single_syndrome[degree] == 0){
                            degree_sum_2[degree]++;
                            single_syndrome[degree] == 1;

                        }
                        
                    }
                }
                
            }
            shots -= shots_this_round;
        }

        degree_sum_ac.fill(0);
        MPI_Allreduce(&degree_sum, &degree_sum_ac, 100, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        degree_sum_2_ac.fill(0);
        MPI_Allreduce(&degree_sum_2, &degree_sum_2_ac, 100, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);


        complex_nodes_sum_ac.fill(0);
        MPI_Allreduce(&complex_nodes_sum, &complex_nodes_sum_ac, 100, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        
        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << std::endl;
        }
        for(uint l = 0;  l< 100; l++){   
            if (degree_sum_ac[l] ) {
                fp_t bucket_freq = ((fp_t)degree_sum_ac[l]/simulator.shots_per_k[k+min_k]);
                degree_distr[l] += bucket_freq*simulator.prob_k[k+min_k];
            }
        }
        for(uint l = 0;  l< 100; l++){   
            if (degree_sum_2_ac[l] ) {
                fp_t bucket_freq = ((fp_t)degree_sum_2_ac[l]/simulator.shots_per_k[k+min_k]);
                degree_distr_2[l] += bucket_freq*simulator.prob_k[k+min_k];
            }
        }
        for(uint l = 0;  l< 100; l++){   
            if (complex_nodes_sum_ac[l] ) {
                fp_t bucket_freq = ((fp_t)complex_nodes_sum_ac[l]/simulator.shots_per_k[k+min_k]);
                complex_nodes_distr[l] += bucket_freq*simulator.prob_k[k+min_k];
            }
        }
    }

    if(world_rank == 0){
        out << "Distance: "<< distance << " Expected LER = " << expected_ler <<" - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;

        for(uint l = 0; l< 100; l++){
            if(degree_distr[l] !=0){
                out << l << ": " <<degree_distr[l] << std::endl;
            }
        }

        out << "_________________ complex nodes _____________________" << std::endl;

        for(uint l = 0; l< 100; l++){
            if(complex_nodes_distr[l] !=0){
                out << l << ": " <<complex_nodes_distr[l] << std::endl;
            }
        } 

        out << "_________________No boundary:_____________________" << std::endl;

        for(uint l = 0; l< 100; l++){
            if(degree_distr_2[l] !=0){
                out << l << ": " <<degree_distr_2[l] << std::endl;
            }
        }
   
    }

}


void edge_distributions(std::ofstream& out, uint64_t max_shot, uint distance, 
    fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    min_k = 0;
    std::array<fp_t,500> edge_distr;
    edge_distr.fill(0);
    std::array<uint64_t,500> edge_sum;
    edge_sum.fill(0);
    std::array<uint64_t,500> edge_sum_ac;
    edge_sum_ac.fill(0);

    uint64_t syndrome_edges;
    uint64_t max_edges_local = 0;
    uint64_t sum_edges_local = 0;

    fp_t average_edge_total = 0;
    fp_t average_edge_per_k = 0;


    uint64_t incorrect_per_k_local = 0;
    fp_t incorrect_per_k = 0;

    uint64_t accurate_per_k_local = 0;
    fp_t accurate_per_k = 0;

    uint64_t skipped_per_k_local = 0;
    fp_t skipped_per_k = 0;

    std::array<uint64_t,50> incorrect_per_hw_local;
    incorrect_per_hw_local.fill(0);
    std::array<uint64_t,50> accurate_per_hw_local;
    accurate_per_hw_local.fill(0);
    std::array<uint64_t,50> skipped_per_hw_local;
    skipped_per_hw_local.fill(0);
    std::array<uint64_t,50> incorrect_per_hw;
    incorrect_per_hw.fill(0);
    std::array<uint64_t,50> accurate_per_hw;
    accurate_per_hw.fill(0);
    std::array<uint64_t,50> skipped_per_hw;
    skipped_per_hw.fill(0);


    std::array<fp_t, 50> incorrect_per_hw_acc;
    incorrect_per_hw_acc.fill(0);
    std::array<fp_t, 50> accurate_per_hw_acc;;
    accurate_per_hw_acc.fill(0);
    std::array<fp_t,50> skipped_per_hw_acc;;
    skipped_per_hw_acc.fill(0);

    std::array<fp_t, 50> incorrect_per_hw_normalized;
    incorrect_per_hw_normalized.fill(0);
    std::array<fp_t, 50> accurate_per_hw_normalized;;
    accurate_per_hw_normalized.fill(0);
    std::array<fp_t,50> skipped_per_hw_normalized;;
    skipped_per_hw_normalized.fill(0);

    std::array<uint64_t,50> n_hw_local;
    n_hw_local.fill(0);
    std::array<fp_t,50> n_hw_all;
    n_hw_all.fill(0);

    std::array<uint64_t,50> clique_n_hw_local;
    clique_n_hw_local.fill(0);
    std::array<fp_t,50> clique_n_hw_all;
    clique_n_hw_all.fill(0);

    std::array<uint64_t,50> smith_n_hw_local;
    smith_n_hw_local.fill(0);
    std::array<fp_t,50> smith_n_hw_all;
    smith_n_hw_all.fill(0);

    std::array<uint64_t,50> n_hw;
    n_hw.fill(0);

    

    

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }
    else if(distance == 15){
        max_k = 21;
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
  
    std::mt19937_64 rng(world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);
    auto decoder = new qrc::MWPMDecoder(circ);
    predecoder->set_decoder(decoder);
    uint n_detectors = circ.count_detectors();
    
    fp_t sum_normalized_prob = 0;
        for(uint k = 0; k < max_k-min_k; k++){
            sum_normalized_prob += simulator.prob_k[k+min_k];
        }


    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, hshots_replc, lshots_rplac);
    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }

    for(uint k = 0; k < max_k-min_k; k++){
        incorrect_per_k_local = 0;
        accurate_per_k_local = 0;
        skipped_per_k_local = 0;
        sum_edges_local = 0;
        incorrect_per_hw_local.fill(0);
        accurate_per_hw_local.fill(0);
        skipped_per_hw_local.fill(0);
        n_hw_local.fill(0);
        smith_n_hw_local.fill(0);
        clique_n_hw_local.fill(0);

        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);

            
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                uint hw = std::accumulate(syndrome.begin(), syndrome.begin()+n_detectors,0);
                n_hw_local[hw]++;
                /////////
                // std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);
                std::vector<std::pair<uint, uint>> smith_map;
                smith_map.clear();
                auto res = decoder->decode_error(syndrome);


                auto res_clique = predecoder->clique_predecoder(syndrome, false);
                uint clique_hw = std::accumulate(res_clique.post_syndrome.begin(), res_clique.post_syndrome.begin()+n_detectors,0);
                clique_n_hw_local[clique_hw]++;


                auto res_smith = predecoder->smith_predecoder(syndrome, true);
                uint smith_hw = std::accumulate(res_smith.post_syndrome.begin(), res_smith.post_syndrome.begin()+n_detectors,0);
                smith_n_hw_local[smith_hw]++;


                syndrome_edges = 0;
                predecoder->update_adjacency_matrix_and_detector_list(syndrome, 1, true);
                auto g = predecoder->get_adjacency_matrix();

                for(const auto& det : g.detector_list){
                    auto v = predecoder->decoding_graph.get_vertex(det);
                    if(g.adjacency_matrix.count(v) == 0){
                        continue;
                    }
                    else{
                        uint degree = g.adjacency_matrix[v].size();
                        syndrome_edges += degree;
                        
                    }
                }
                syndrome_edges /=2;
                if(max_edges_local < syndrome_edges){
                    max_edges_local = syndrome_edges;
                }
                sum_edges_local += syndrome_edges;


                for(const auto& det : g.detector_list){
                    auto v = predecoder->decoding_graph.get_vertex(det);
                    if(g.adjacency_matrix.count(v) == 0){
                        continue;
                    }
                    auto neighbors = g.adjacency_matrix[v];
                    for(const auto& n : neighbors){
                        smith_map.push_back(std::make_pair(v->detector,n->detector));
                    }
                    
                }
                syndrome_edges /=2;
                if(max_edges_local < syndrome_edges){
                    max_edges_local = syndrome_edges;
                }
                sum_edges_local += syndrome_edges;

                auto compare = qpd::accuracy_coverage(res.matching, smith_map);
                incorrect_per_k_local += compare.countBNotInA;
                accurate_per_k_local += compare.countBothAB;
                skipped_per_k_local += compare.countANotInB;

                incorrect_per_hw_local[hw] += compare.countBNotInA;
                accurate_per_hw_local[hw] += compare.countBothAB;
                skipped_per_hw_local[hw] += compare.countANotInB;


                
            }
            shots -= shots_this_round;
        }

        uint64_t max_edges = 0;
        MPI_Allreduce(&max_edges_local, &max_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

        uint64_t sum_edges= 0;
        MPI_Allreduce(&sum_edges_local, &sum_edges, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t incorrect= 0;
        MPI_Allreduce(&incorrect_per_k_local, &incorrect, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t accurate= 0;
        MPI_Allreduce(&accurate_per_k_local, &accurate, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t skipped= 0;
        MPI_Allreduce(&skipped_per_k_local, &skipped, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        incorrect_per_hw.fill(0);
        MPI_Allreduce(incorrect_per_hw_local.data(), incorrect_per_hw.data(), incorrect_per_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        accurate_per_hw.fill(0);
        MPI_Allreduce(accurate_per_hw_local.data(), accurate_per_hw.data(), accurate_per_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        skipped_per_hw.fill(0);
        MPI_Allreduce(skipped_per_hw_local.data(), skipped_per_hw.data(), skipped_per_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        

        n_hw.fill(0);
        MPI_Allreduce(n_hw_local.data(), n_hw.data(), n_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t,50> clique_n_hw;
        clique_n_hw.fill(0);
        MPI_Allreduce(clique_n_hw_local.data(), clique_n_hw.data(), clique_n_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t,50> smith_n_hw;
        smith_n_hw.fill(0);
        MPI_Allreduce(smith_n_hw_local.data(), smith_n_hw.data(), smith_n_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }

        average_edge_per_k = ((fp_t)sum_edges/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];        
        average_edge_total += ((fp_t)sum_edges/(simulator.shots_per_k[k+min_k]))*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        incorrect_per_k = ((fp_t)incorrect/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];
        accurate_per_k = ((fp_t)accurate/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];
        skipped_per_k = ((fp_t)skipped/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];

        for(uint h=0;h<n_hw.size();h++){
            if(n_hw[h] != 0){
                incorrect_per_hw_acc[h] += ((fp_t)incorrect_per_hw[h]/n_hw[h])*simulator.prob_k[k+min_k];
                incorrect_per_hw_normalized[h] += ((fp_t)incorrect_per_hw[h]/n_hw[h])*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

                accurate_per_hw_acc[h] += ((fp_t)accurate_per_hw[h]/n_hw[h])*simulator.prob_k[k+min_k];
                accurate_per_hw_normalized[h] += ((fp_t)accurate_per_hw[h]/n_hw[h])*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

                skipped_per_hw_acc[h] += ((fp_t)skipped_per_hw[h]/n_hw[h])*simulator.prob_k[k+min_k];
                skipped_per_hw_normalized[h] += ((fp_t)skipped_per_hw[h]/n_hw[h])*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

                n_hw_all[h] += (((fp_t)n_hw[h]/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k]);
                clique_n_hw_all[h] += (((fp_t)clique_n_hw[h]/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k]);
                smith_n_hw_all[h] += (((fp_t)smith_n_hw[h]/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k]);
            }

        }

        if(world_rank == 0 ){
            std::cout << "\t ----------Edge Info-------" << std::endl;
            std::cout << "\t\t Avg #Edges (*prob) = " << ((fp_t)sum_edges/(simulator.shots_per_k[k+min_k]))  <<
                         "(" << average_edge_per_k <<")"<<std::endl;
            std::cout << "\tNormalized:" << std::endl;
            std::cout << "\t\t Avg #Edges = " << average_edge_total << std::endl; 
            
            std::cout << "\tMAX Edge = " << max_edges  << std::endl;
            std::cout << "\t ----------Accuracy and Coverage-------"<< std::endl;
            std::cout << "\t This bucket" << std::endl;
            std::cout << "\t\t accurate(*prob): " <<  ((fp_t)incorrect/(simulator.shots_per_k[k+min_k])) << "(" <<
            accurate_per_k << ")"<< std::endl;
            std::cout << "\t\t incorrect(*prob): " <<  ((fp_t)accurate/(simulator.shots_per_k[k+min_k])) << "(" <<
            incorrect_per_k << ")"<< std::endl;
            std::cout << "\t\t skipped(*prob): " <<  ((fp_t)skipped/(simulator.shots_per_k[k+min_k])) << "(" <<
            skipped_per_k << ")"<< std::endl;

            std::cout << "_______________________________________________"<< std::endl;
        }




    }

    if(world_rank == 0 ){
        std::cout << "  -------------COVERAGE AND ACC----------" << std::endl;
        for(uint h=0;h<n_hw.size();h++){
            if(accurate_per_hw_acc[h]!=0 || incorrect_per_hw_acc[h]!=0 ||  skipped_per_hw_acc[h]!=0)
            std::cout << h << ": " << accurate_per_hw_acc[h] << " - " << incorrect_per_hw_acc[h] << " - " << skipped_per_hw_acc[h] << std::endl;
        }

        std::cout << "  -------------COVERAGE AND ACC (NORMALIZED)----------" << std::endl;
        for(uint h=0;h<n_hw.size();h++){
            if(accurate_per_hw_acc[h]!=0 || incorrect_per_hw_acc[h]!=0 ||  skipped_per_hw_acc[h]!=0)
            std::cout << h << ": " << accurate_per_hw_normalized[h] << " - " << incorrect_per_hw_normalized[h] << " - " << skipped_per_hw_normalized[h] << std::endl;
        }
        std::cout << "  -------------HW Distribution----------" << std::endl;
        for(uint h=0;h<n_hw.size();h++){
            if(n_hw_all[h] != 0)
                std::cout << h << ": " << (n_hw_all[h]) << std::endl;
        }
        std::cout << "  -------------HW After Promatch Distribution----------" << std::endl;
        std::array<fp_t, 11> hw_promatch;
        hw_promatch.fill(0);
        for(uint h=0;h<n_hw_all.size();h++){
            if(n_hw_all[h] != 0){
                if(h>10){
                    fp_t hw_9_10 = 0.995985*n_hw_all[h];
                    fp_t hw_7_8 = 0.003974*n_hw_all[h];
                    fp_t hw_5_6 = 4.1e-05*n_hw_all[h];
                    if(h % 2 == 1){
                        hw_promatch[9] +=hw_9_10;
                        hw_promatch[7] +=hw_7_8;
                        hw_promatch[5] +=hw_9_10;
                    }
                    else{
                        hw_promatch[10] +=hw_9_10;
                        hw_promatch[8] +=hw_7_8;
                        hw_promatch[6] +=hw_9_10;
                    }
                }
                else{
                    hw_promatch[h] += n_hw_all[h];
                }
                
            }
        }

        for(uint h=0;h<hw_promatch.size();h++){
            std::cout << h << ": " << (hw_promatch[h]) << std::endl;
        }

        std::cout << "  -------------HW After Clique Distribution----------" << std::endl;
        std::array<fp_t,50> hw_clique;
        hw_clique.fill(0);
        for(uint h=0;h<n_hw_all.size();h++){
            if(n_hw_all[h] != 0){
                if(h<=10){
                    hw_clique[h] = n_hw_all[h];
                }
                hw_clique[h] += clique_n_hw_all[h];                
            }
        }

        for(uint h=0;h<hw_clique.size();h++){
            if(n_hw_all[h] != 0){
                std::cout << h << ": " << (hw_clique[h]) << std::endl;
            }
            
        }

        std::cout << "  -------------HW After Smith Distribution----------" << std::endl;
        std::array<fp_t, 50> hw_smith;
        hw_smith.fill(0);
        for(uint h=0;h<n_hw_all.size();h++){
            if(n_hw_all[h] != 0){
                if(h<=10){
                    hw_smith[h] = n_hw_all[h];
                }
                hw_smith[h] += smith_n_hw_all[h];                
            }
        }

        for(uint h=0;h<hw_smith.size();h++){
            if(n_hw_all[h] != 0){
                std::cout << h << ": " << (hw_smith[h]) << std::endl;
            }
            
        }
    }

}

void HW_distributions(std::ofstream& out, uint64_t max_shot, uint distance, 
    fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    min_k = 0;

    std::array<uint64_t,50> n_hw_local;
    n_hw_local.fill(0);
    std::array<fp_t,50> n_hw_all;
    n_hw_all.fill(0);

    std::array<uint64_t,50> clique_n_hw_local;
    clique_n_hw_local.fill(0);
    std::array<fp_t,50> clique_n_hw_all;
    clique_n_hw_all.fill(0);

    std::array<uint64_t,50> smith_n_hw_local;
    smith_n_hw_local.fill(0);
    std::array<fp_t,50> smith_n_hw_all;
    smith_n_hw_all.fill(0);

    std::array<uint64_t,50> n_hw;
    n_hw.fill(0);

    

    

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }
    else if(distance == 15){
        max_k = 21;
    }
    bbsim::BucketBasedSim simulator(circ,max_k);
  
    std::mt19937_64 rng(world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ);
    auto decoder = new qrc::MWPMDecoder(circ);
    predecoder->set_decoder(decoder);
    uint n_detectors = circ.count_detectors();
    
    fp_t sum_normalized_prob = 0;
        for(uint k = 0; k < max_k-min_k; k++){
            sum_normalized_prob += simulator.prob_k[k+min_k];
        }


    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, hshots_replc, lshots_rplac);
    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }

    for(uint k = 0; k < max_k-min_k; k++){
        n_hw_local.fill(0);
        smith_n_hw_local.fill(0);
        clique_n_hw_local.fill(0);

        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);

            
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                uint hw = std::accumulate(syndrome.begin(), syndrome.begin()+n_detectors,0);
                n_hw_local[hw]++;
                /////////
                // std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);
                std::vector<std::pair<uint, uint>> smith_map;
                smith_map.clear();

                if(hw > 10){
                    auto res_clique = predecoder->clique_predecoder(syndrome, false);
                    uint clique_hw = std::accumulate(res_clique.post_syndrome.begin(), res_clique.post_syndrome.begin()+n_detectors,0);
                    clique_n_hw_local[clique_hw]++;


                    auto res_smith = predecoder->smith_predecoder(syndrome, true);
                    uint smith_hw = std::accumulate(res_smith.post_syndrome.begin(), res_smith.post_syndrome.begin()+n_detectors,0);
                    smith_n_hw_local[smith_hw]++;
                }
                

                
            }
            shots -= shots_this_round;
        }

        n_hw.fill(0);
        MPI_Allreduce(n_hw_local.data(), n_hw.data(), n_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t,50> clique_n_hw;
        clique_n_hw.fill(0);
        MPI_Allreduce(clique_n_hw_local.data(), clique_n_hw.data(), clique_n_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t,50> smith_n_hw;
        smith_n_hw.fill(0);
        MPI_Allreduce(smith_n_hw_local.data(), smith_n_hw.data(), smith_n_hw_local.size(),
            MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }


        for(uint h=0;h<n_hw.size();h++){
                n_hw_all[h] += (((fp_t)n_hw[h]/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k]);
                clique_n_hw_all[h] += (((fp_t)clique_n_hw[h]/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k]);
                smith_n_hw_all[h] += (((fp_t)smith_n_hw[h]/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k]);

        }

        if(world_rank == 0 ){
            std::cout << "_______________________________________________"<< std::endl;
        }




    }

    if(world_rank == 0 ){
        std::cout << "  -------------HW Distribution----------" << std::endl;
        for(uint h=0;h<n_hw_all.size();h++){
            if(n_hw_all[h] != 0)
                std::cout << h << ": " << (n_hw_all[h]) << std::endl;
        }
        std::cout << "  -------------HW After Promatch Distribution----------" << std::endl;
        std::array<fp_t, 11> hw_promatch;
        hw_promatch.fill(0);
        for(uint h=0;h<n_hw_all.size();h++){
            if(n_hw_all[h] != 0){
                if(h>10){
                    fp_t hw_9_10 = 0.995985*n_hw_all[h];
                    fp_t hw_7_8 = 0.003974*n_hw_all[h];
                    fp_t hw_5_6 = 4.1e-05*n_hw_all[h];
                    if(h % 2 == 1){
                        hw_promatch[9] +=hw_9_10;
                        hw_promatch[7] +=hw_7_8;
                        hw_promatch[5] +=hw_9_10;
                    }
                    else{
                        hw_promatch[10] +=hw_9_10;
                        hw_promatch[8] +=hw_7_8;
                        hw_promatch[6] +=hw_9_10;
                    }
                }
                else{
                    hw_promatch[h] += n_hw_all[h];
                }
                
            }
        }

        for(uint h=0;h<hw_promatch.size();h++){
            std::cout << h << ": " << (hw_promatch[h]) << std::endl;
        }

        std::cout << "  -------------HW After Clique Distribution----------" << std::endl;
        std::array<fp_t,50> hw_clique;
        hw_clique.fill(0);
        for(uint h=0;h<n_hw_all.size();h++){
            if(n_hw_all[h] != 0){
                if(h<=10){
                    hw_clique[h] = n_hw_all[h];
                }
                hw_clique[h] += clique_n_hw_all[h];                
            }
        }

        for(uint h=0;h<hw_clique.size();h++){
            if(n_hw_all[h] != 0){
                std::cout << h << ": " << (hw_clique[h]) << std::endl;
            }
            
        }

        std::cout << "  -------------HW After Smith Distribution----------" << std::endl;
        std::array<fp_t, 50> hw_smith;
        hw_smith.fill(0);
        for(uint h=0;h<n_hw_all.size();h++){
            if(n_hw_all[h] != 0){
                if(h<=10){
                    hw_smith[h] = n_hw_all[h];
                }
                hw_smith[h] += smith_n_hw_all[h];                
            }
        }

        for(uint h=0;h<hw_smith.size();h++){
            if(n_hw_all[h] != 0){
                std::cout << h << ": " << (hw_smith[h]) << std::endl;
            }
            
        }
    }

}


void ler_predecoders(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, fp_t threshold_scale, bool& print_time, std::string decoder_name, bool generalized, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string syndrome_clock_challenging_file, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if(distance == 11){
        min_k = 10;
        max_k = 13;
    }
    else if(distance == 13){
        min_k = 9;
        max_k = 12;
    }

    // std::string syndrome_file = "../NFDecoder/data/challengingSyndromes_d13/challenging_d13.txt";
    if(world_rank == 0){
        std::cout << "round#" <<round_n <<" Promatch-A + " << decoder_name <<std::endl;
    }



    // for saving syndromes
    uint64_t round = 1;
    uint64_t round_t = 1;
    if (world_rank == 0 && save_syndromes){
        std::cout << syndrome_folder_name<< std::endl;
        // To check if we are reading correct files
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

    
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;
    qrc::benchmark::StatisticalResult statres;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_size*round_n+world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ, generalized);
    predecoder->set_decoder(decoder);

    //////
    std::vector<vector<vector<uint16_t>>> saved_hhw_syndromes;
    saved_hhw_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits;
    //////
    /// For saving time challenging syndromes
    std::vector<vector<vector<uint16_t>>> saved_timechallenging_syndromes;
    saved_timechallenging_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits_t;
    ////



    ///// ================== For parallel running with AstreaG ========================

    qrc::AstreaParams astreaG_param= {};
    qrc::Decoder* AstreaG_decoder;
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    astreaG_param.use_mld = false;

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(threshold_scale*0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    if(world_rank == 0){
        std::cout << "threshold_scale:" << threshold_scale << std::endl;
    }


    AstreaG_decoder =  new  qrc::Astrea(circ,
                                        n_detectors_per_round,
                                        weight_filter_cutoff,
                                        astreaG_param);


    ///// ==============================================================================

    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, hshots_replc, lshots_rplac);

    uint64_t max_k_ = max_k;
    fp_t sum_normalized_prob = 0;
    for(uint k = 0; k < max_k-min_k; k++){
        sum_normalized_prob += simulator.prob_k[k+min_k];

        // To eliminate the buckets that their probability is lower than LER and accordinly has zero number of shots
        if(simulator.shots_per_k[k+min_k] == 0){
            max_k_ --;
        }
    }
    max_k = max_k_;

    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << " Physical Err = "<< physcial_error<<std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }
    fp_t prev_logical_error_rate = 0.0;

    
    std::array<fp_t,6> normalized_stages;
    normalized_stages.fill(0);

    std::array<fp_t,6> weighted_stages;
    weighted_stages.fill(0);

    ///////Clock Simulation
    uint64_t local_sum_cycle = 0;

    uint64_t local_max_cycle = 0;

    uint64_t local_max_round = 0;

    uint64_t local_sum_cycle_parallel = 0;

    uint64_t local_max_cycle_parallel = 0;


    fp_t final_average_cycle = 0;
    fp_t final_average_cycle_normalized = 0;
    fp_t final_average_cycle_parallel = 0;
    fp_t final_average_cycle_parallel_normalized = 0;

    uint64_t local_number_above_threshold_parallel = 0;
    uint64_t local_number_above_threshold = 0;

    uint64_t local_n_Astrea10_used = 0;
    uint64_t local_n_Astrea8_used = 0;
    uint64_t local_n_Astrea6_used = 0;

    fp_t total_Astrea10_used = 0;
    fp_t total_Astrea8_used = 0;
    fp_t total_Astrea6_used = 0;

    fp_t time_constrained_ler = 0;

    fp_t ler_PAG = 0;
    fp_t ler_clique = 0;
    fp_t ler_smith = 0;
    fp_t ler_clique_AG = 0;
    fp_t ler_smith_AG = 0;
    fp_t ler_AG = 0;

    uint64_t number_of_logical_err_PAG = 0;
    uint64_t number_of_logical_err_clique = 0;
    uint64_t number_of_logical_err_smith = 0;
    uint64_t number_of_logical_err_clique_AG = 0;
    uint64_t number_of_logical_err_smith_AG = 0;
    uint64_t number_of_logical_err_AG = 0;
    uint64_t local_number_above_threshold_parallel_prior = 0;


    //////Clock Simulation


    //////
    // std::vector<std::vector<uint16_t>> vecss = qpd::read_vector_of_vectors(syndrome_file);
    // std::cout << "Size of incorrect pairs = " << vecss.size() << std::endl;
    //////
    // uint k=0;
    // for(auto compressed_synd : vecss){
    //     k++;
    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t local_errors = 0;
        uint64_t local_errors_smith = 0;
        uint64_t local_errors_clique = 0;
        uint64_t local_errors_PAG = 0;
        uint64_t local_errors_smith_AG = 0;
        uint64_t local_errors_clique_AG = 0;
        uint64_t local_errors_AG = 0;
        
        std::array<uint64_t, 6> local_stages;
        local_stages.fill(0);

        uint64_t local_times_lhw = 0;

        local_sum_cycle_parallel = 0;
        local_sum_cycle = 0;

        local_n_Astrea10_used = 0;
        local_n_Astrea8_used = 0;
        local_n_Astrea6_used = 0;


        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                uint hw = std::accumulate(syndrome.begin(), syndrome.begin()+predecoder->n_detectors,0);
                if(hw <= MAX_BF10_HW){
                    local_times_lhw++;
                }
            /////////
            // std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);

            /////////
                qrc::DecoderShotResult res;
                qpd::PredecoderDecoderShotResult res_smith;
                qpd::PredecoderDecoderShotResult res_clique;

                qrc::DecoderShotResult res_astreag;


                res = predecoder->adaptively_decode_error(syndrome);
                res_smith = predecoder->smith_decode_error(syndrome);
                res_clique = predecoder->clique_decode_error(syndrome);

                res_astreag = AstreaG_decoder->decode_error(syndrome);
                

                fp_t astreag_weight = predecoder->calc_matching_weight(res_astreag.matching);
                fp_t promatch_weight = predecoder->calc_matching_weight(res.matching);
                local_errors+= res.is_logical_error;
                bool PAG_err = (astreag_weight < promatch_weight)? res_astreag.is_logical_error :res.is_logical_error;
                local_errors_AG += res_astreag.is_logical_error;
                local_errors_PAG += PAG_err;
                local_errors_clique += res_clique.is_logical_error;
                local_errors_smith += res_smith.is_logical_error;
                local_errors_clique_AG += (astreag_weight < res_clique.weight)? res_astreag.is_logical_error :res_clique.is_logical_error;
                local_errors_smith_AG += (astreag_weight < res_smith.weight)? res_astreag.is_logical_error :res_smith.is_logical_error;
                // local_errors_clique += res_clique.islo
                for(uint stg = 0; stg <6; stg++){
                    local_stages[stg] += predecoder->reached_stage[stg];
                }
                
                if(res.is_logical_error){
                    //std::cout << "Error" << k <<"-";
                    flipped_bits = qpd::syndrome_compressed(syndrome);
                }
                if(PAG_err){
                    std::cout << "* ";
                    flipped_bits = qpd::syndrome_compressed(syndrome);
                }

                if(predecoder->total_cycles > MAX_BF6_CYCLE){
                    local_number_above_threshold++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                if(predecoder->total_cycles_parallel_updating > MAX_BF6_CYCLE){
                    local_number_above_threshold_parallel++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                

                //Setting clock related varaibles
                if (local_max_cycle < predecoder->total_cycles){
                    local_max_cycle = predecoder->total_cycles;
                }
                if(local_max_cycle_parallel < predecoder->total_cycles_parallel_updating){
                    local_max_cycle_parallel = predecoder->total_cycles_parallel_updating;
                }
                local_sum_cycle += predecoder->total_cycles;
                local_sum_cycle_parallel += predecoder->total_cycles_parallel_updating;

                if(local_max_round < predecoder->number_of_rounds){
                    local_max_round = predecoder->number_of_rounds;
                }
                if(predecoder->MAX_BF_HW == MAX_BF10_HW && hw > MAX_BF10_HW){
                    local_n_Astrea10_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF8_HW && hw > MAX_BF10_HW){
                    local_n_Astrea8_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF6_HW && hw > MAX_BF10_HW){
                    local_n_Astrea6_used++;
                }
            }
            shots -= shots_this_round;
        }
        uint64_t fault_errors = 0;
        MPI_Allreduce(&local_errors, &fault_errors, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_smith = 0;
        MPI_Allreduce(&local_errors_smith, &fault_errors_smith, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_clique = 0;
        MPI_Allreduce(&local_errors_clique, &fault_errors_clique, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_PAG = 0;
        MPI_Allreduce(&local_errors_PAG, &fault_errors_PAG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_smith_AG = 0;
        MPI_Allreduce(&local_errors_smith_AG, &fault_errors_smith_AG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_clique_AG = 0;
        MPI_Allreduce(&local_errors_clique_AG, &fault_errors_clique_AG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_AG = 0;
        MPI_Allreduce(&local_errors_AG, &fault_errors_AG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);


        uint64_t times_lhw = 0;
        MPI_Allreduce(&local_times_lhw, &times_lhw, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t, 6> stages;
        stages.fill(0);
        MPI_Allreduce(&local_stages, &stages, 6, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle = 0;
        MPI_Allreduce(&local_max_cycle, &max_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle = 0;
        MPI_Allreduce(&local_sum_cycle, &sum_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle_parallel = 0;
        MPI_Allreduce(&local_max_cycle_parallel, &max_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle_parallel = 0;
        MPI_Allreduce(&local_sum_cycle_parallel, &sum_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold = 0;
        MPI_Allreduce(&local_number_above_threshold, &number_above_threshold, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold_parallel = 0;
        MPI_Allreduce(&local_number_above_threshold_parallel, &number_above_threshold_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_round = 0;
        MPI_Allreduce(&local_max_round, &max_round, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        
        uint64_t n_Astrea10_used = 0;
        MPI_Allreduce(&local_n_Astrea10_used, &n_Astrea10_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea8_used = 0;
        MPI_Allreduce(&local_n_Astrea8_used, &n_Astrea8_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea6_used = 0;
        MPI_Allreduce(&local_n_Astrea6_used, &n_Astrea6_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }

        if (fault_errors > 0 || ((number_above_threshold_parallel -  local_number_above_threshold_parallel_prior) > 0)) {
            fp_t failure_rate = ((fp_t)fault_errors/simulator.shots_per_k[k+min_k]);
        
            prev_logical_error_rate = statres.logical_error_rate;
            statres.logical_error_rate += (failure_rate*simulator.prob_k[k+min_k]);
            statres.n_logical_errors += fault_errors;
            fp_t time_failure_rate = (fp_t)(number_above_threshold_parallel -  local_number_above_threshold_parallel_prior)/simulator.shots_per_k[k+min_k];
            time_constrained_ler += (failure_rate*simulator.prob_k[k+min_k]) + (time_failure_rate*simulator.prob_k[k+min_k]);
            local_number_above_threshold_parallel_prior = number_above_threshold_parallel;

            if (world_rank == 0) {
                std::cout << "===== Promatch =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate << std::endl
                        << "\tNumber of errors = " << statres.n_logical_errors << std::endl;
                std::cout << "\tLogical error rate = " << statres.logical_error_rate << std::endl;
                std::cout << "\t\t Logical error rate timed = " << time_constrained_ler << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_PAG > 0) {
            fp_t failure_rate_PAG = ((fp_t)fault_errors_PAG/simulator.shots_per_k[k+min_k]);
        
            ler_PAG += (failure_rate_PAG*simulator.prob_k[k+min_k]);
            number_of_logical_err_PAG += fault_errors_PAG;
        
            if (world_rank == 0) {
                std::cout << "===== Promatch || AG =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_PAG << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_PAG << std::endl;
                std::cout << "\tLogical error rate = " << ler_PAG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_AG > 0) {
            fp_t failure_rate_AG = ((fp_t)fault_errors_AG/simulator.shots_per_k[k+min_k]);
        
            ler_AG += (failure_rate_AG*simulator.prob_k[k+min_k]);
            number_of_logical_err_AG += fault_errors_AG;
        
            if (world_rank == 0) {
                std::cout << "===== AstreaG =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_AG << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_AG << std::endl;
                std::cout << "\tLogical error rate = " << ler_AG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_clique > 0) {
            fp_t failure_rate_clique = ((fp_t)fault_errors_clique/simulator.shots_per_k[k+min_k]);
        
            ler_clique += (failure_rate_clique*simulator.prob_k[k+min_k]);
            number_of_logical_err_clique += fault_errors_clique;
        
            if (world_rank == 0) {
                std::cout << "===== Clique =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_clique << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_clique << std::endl;
                std::cout << "\tLogical error rate = " << ler_clique << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_clique_AG > 0) {
            fp_t failure_rate_clique_AG = ((fp_t)fault_errors_clique_AG/simulator.shots_per_k[k+min_k]);
        
            ler_clique_AG += (failure_rate_clique_AG*simulator.prob_k[k+min_k]);
            number_of_logical_err_clique_AG += fault_errors_clique_AG;
        
            if (world_rank == 0) {
                std::cout << "===== Clique || AG =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_clique_AG << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_clique_AG << std::endl;
                std::cout << "\tLogical error rate = " << ler_clique_AG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_smith > 0) {
            fp_t failure_rate_smith = ((fp_t)fault_errors_smith/simulator.shots_per_k[k+min_k]);
        
            ler_smith += (failure_rate_smith*simulator.prob_k[k+min_k]);
            number_of_logical_err_smith += fault_errors_smith;
        
            if (world_rank == 0) {
                std::cout << "===== Smith =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_smith << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_smith << std::endl;
                std::cout << "\tLogical error rate = " << ler_smith << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_smith_AG > 0) {
            fp_t failure_rate_smith_AG = ((fp_t)fault_errors_smith_AG/simulator.shots_per_k[k+min_k]);
        
            ler_smith_AG += (failure_rate_smith_AG*simulator.prob_k[k+min_k]);
            number_of_logical_err_smith_AG += fault_errors_smith_AG;
        
            if (world_rank == 0) {
                std::cout << "===== Smith || AG =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_smith_AG << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_smith_AG << std::endl;
                std::cout << "\tLogical error rate = " << ler_smith_AG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        
        std::array<fp_t,6> stages_rate;
        stages_rate.fill(0);
        for(uint stg = 0; stg < 6; stg++){
            stages_rate[stg] = ((fp_t)stages[stg]/(simulator.shots_per_k[k+min_k]-times_lhw));
            normalized_stages[stg] += stages_rate[stg]*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
            weighted_stages[stg] += stages_rate[stg]*simulator.prob_k[k+min_k];
 
        }

        // std::cout << world_rank <<" sum_cycle"  <<sum_cycle<< " sum_cycle_parallel: " <<sum_cycle_parallel<<std::endl;

        final_average_cycle += ((fp_t)sum_cycle/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];        
        final_average_cycle_normalized += ((fp_t)sum_cycle/(simulator.shots_per_k[k+min_k]))*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        final_average_cycle_parallel += ((fp_t)sum_cycle_parallel/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];
        final_average_cycle_parallel_normalized += ((fp_t)sum_cycle_parallel/(simulator.shots_per_k[k+min_k]))*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        fp_t percentage_astrea10 = ((fp_t)n_Astrea10_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea10_used +=percentage_astrea10*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
        
        fp_t percentage_astrea8 = ((fp_t)n_Astrea8_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea8_used +=percentage_astrea8*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
        
        fp_t percentage_astrea6 = ((fp_t)n_Astrea6_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea6_used +=percentage_astrea6*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);


        // std::cout << "world_rank: "<<world_rank  << " final_average_cycle: " << final_average_cycle <<
        // " final_average_cycle_normalized: " << final_average_cycle_normalized << " final_average_cycle_parallel: " << final_average_cycle_parallel << 
        // " final_average_cycle_parallel_normalized: " << final_average_cycle_parallel_normalized << std::endl;
        
        if (world_rank == 0) {
            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\t\tS1: " << stages_rate[0] << " S2_1: " << stages_rate[1] << " S2_2: " << stages_rate[2] << " S3: " << stages_rate[3] << " S4_1: " << stages_rate[4] 
                            << " S4_2: " << stages_rate[5] << std::endl;
            std::cout << "\tNormalized Stage Rates:" << std::endl;
            std::cout << "\t\tS1: " << normalized_stages[0] << " S2_1: " << normalized_stages[1] << " S2_2: " << normalized_stages[2] << " S3: " << normalized_stages[3] << " S4_1: " << normalized_stages[4] 
                            << " S4_2: " << normalized_stages[5] << std::endl;
            std::cout << "\tStage Rates:" << std::endl;
            std::cout << "\t\tS1: " << weighted_stages[0] << " S2_1: " << weighted_stages[1] << " S2_2: " << weighted_stages[2] << " S3: " << weighted_stages[3] << " S4_1: " << weighted_stages[4] 
                            << " S4_2: " << weighted_stages[5] << std::endl;

            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\tNormalized:" << std::endl;
            std::cout << "\t\t Avg Cycles (UP) = " << final_average_cycle_parallel_normalized  << "(" << final_average_cycle_normalized <<")"<<std::endl;
            std::cout << "\tNot Normalized:" << std::endl;
            std::cout << "\t\t Avg Cycles(UP) = " << final_average_cycle_parallel << "(" << final_average_cycle <<")"<<std::endl;
            std::cout << "\tMAX Cycles (UP) = " << max_cycle_parallel  << "(" << max_cycle << ") Number above threshold(UP):" << number_above_threshold_parallel << "(" << number_above_threshold << ")"<< std::endl;
            std::cout << "\tMAX Round = " << max_round  << std::endl;
            std::cout << "\t ----------ASTREA PERCENTAGES-------"<< std::endl;
            std::cout << "\t This bucket" << std::endl;
            std::cout << "\t\tAstrea10: " << percentage_astrea10 << " Astrea8: " << percentage_astrea8 << " Astrea6: " << percentage_astrea6 << std::endl;
            std::cout << "\t Total(Normalized)" << std::endl;
            std::cout << "\t\tAstrea10: " << total_Astrea10_used << " Astrea8: " << total_Astrea8_used << " Astrea6: " << total_Astrea6_used << std::endl;

            std::cout << "_______________________________________________"<< std::endl;


        }

    }

}

void PAG_synergies(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, bool& print_time, std::string decoder_name, bool generalized, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string syndrome_clock_challenging_file, uint64_t hshots_replc, uint64_t lshots_rplac){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // std::string syndrome_file = "../NFDecoder/data/challengingSyndromes_d13/challenging_d13.txt";
    if(world_rank == 0){
        std::cout << "round#" <<round_n <<" Promatch-A + " << decoder_name <<std::endl;
    }



    // for saving syndromes
    uint64_t round = 1;
    uint64_t round_t = 1;
    if (world_rank == 0 && save_syndromes){
        std::cout << syndrome_folder_name<< std::endl;
        // To check if we are reading correct files
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
    // if(distance == 13){
    //     max_k = 17;
    // } 
    // else if(distance == 11){
    //     max_k = 12;
    // }
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;
    qrc::benchmark::StatisticalResult statres;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_size*round_n+world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ, generalized);
    predecoder->set_decoder(decoder);

    
    /// For saving time challenging syndromes
    std::vector<vector<vector<uint16_t>>> saved_timechallenging_syndromes;
    saved_timechallenging_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits_t;
    ////



    ///// ================== For parallel running with AstreaG ========================

    qrc::AstreaParams astreaG_param= {};
    qrc::Decoder* AstreaG_decoder;
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    astreaG_param.use_mld = false;

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);


    AstreaG_decoder =  new  qrc::Astrea(circ,
                                        n_detectors_per_round,
                                        weight_filter_cutoff,
                                        astreaG_param);


    ///// ==============================================================================

    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, hshots_replc, lshots_rplac);

    uint64_t max_k_ = max_k;
    fp_t sum_normalized_prob = 0;
    for(uint k = 0; k < max_k-min_k; k++){
        sum_normalized_prob += simulator.prob_k[k+min_k];

        // To eliminate the buckets that their probability is lower than LER and accordinly has zero number of shots
        if(simulator.shots_per_k[k+min_k] == 0){
            max_k_ --;
        }
    }
    max_k = max_k_;

    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << " Physical Err = "<< physcial_error<<std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }
    fp_t prev_logical_error_rate = 0.0;

    
    std::array<fp_t,6> normalized_stages;
    normalized_stages.fill(0);

    std::array<fp_t,6> weighted_stages;
    weighted_stages.fill(0);

    ///////Clock Simulation
    uint64_t local_sum_cycle = 0;

    uint64_t local_max_cycle = 0;

    uint64_t local_max_round = 0;

    uint64_t local_sum_cycle_parallel = 0;

    uint64_t local_max_cycle_parallel = 0;


    fp_t final_average_cycle = 0;
    fp_t final_average_cycle_normalized = 0;
    fp_t final_average_cycle_parallel = 0;
    fp_t final_average_cycle_parallel_normalized = 0;

    uint64_t local_number_above_threshold_parallel = 0;
    uint64_t local_number_above_threshold = 0;

    uint64_t local_n_Astrea10_used = 0;
    uint64_t local_n_Astrea8_used = 0;
    uint64_t local_n_Astrea6_used = 0;

    fp_t total_Astrea10_used = 0;
    fp_t total_Astrea8_used = 0;
    fp_t total_Astrea6_used = 0;

    fp_t time_constrained_ler = 0;

    fp_t ler_PAG = 0;
    fp_t ler_AG = 0;

    uint64_t number_of_logical_err_PAG = 0;
    uint64_t number_of_logical_err_AG = 0;
    uint64_t local_number_above_threshold_parallel_prior = 0;

    std::vector<fp_t> vector_sparcity_weighted;
    std::vector<fp_t> vector_sparcity_unweighted;
    fp_t sparcity_uw;


    //////Clock Simulation


    //////
    // std::vector<std::vector<uint16_t>> vecss = qpd::read_vector_of_vectors(syndrome_file);
    // std::cout << "Size of incorrect pairs = " << vecss.size() << std::endl;
    //////
    // uint k=0;
    // for(auto compressed_synd : vecss){
    //     k++;
    uint count_odd =0;
    for(uint k = 0; k < max_k-min_k; k++){
        std::cout << "***********" << count_odd << std::endl;
        uint64_t local_errors = 0;
        uint64_t local_errors_smith = 0;
        uint64_t local_errors_clique = 0;
        uint64_t local_errors_PAG = 0;
        uint64_t local_errors_smith_AG = 0;
        uint64_t local_errors_clique_AG = 0;
        uint64_t local_errors_AG = 0;
        
        std::array<uint64_t, 6> local_stages;
        local_stages.fill(0);

        uint64_t local_times_lhw = 0;

        local_sum_cycle_parallel = 0;
        local_sum_cycle = 0;

        local_n_Astrea10_used = 0;
        local_n_Astrea8_used = 0;
        local_n_Astrea6_used = 0;

        
        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                uint hw = std::accumulate(syndrome.begin(), syndrome.begin()+predecoder->n_detectors,0);
                if(hw <= MAX_BF10_HW){
                    local_times_lhw++;
                }
            /////////
            // std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);

            /////////
                qrc::DecoderShotResult res;
                qpd::PredecoderDecoderShotResult res_smith;
                qpd::PredecoderDecoderShotResult res_clique;

                qrc::DecoderShotResult res_astreag;


                res = predecoder->adaptively_decode_error(syndrome);
                res_smith = predecoder->smith_decode_error(syndrome);
                res_clique = predecoder->clique_decode_error(syndrome);

                res_astreag = AstreaG_decoder->decode_error(syndrome);
                

                fp_t astreag_weight = predecoder->calc_matching_weight(res_astreag.matching);
                fp_t promatch_weight = predecoder->calc_matching_weight(res.matching);
                local_errors+= res.is_logical_error;
                bool PAG_err = (astreag_weight < promatch_weight)? res_astreag.is_logical_error :res.is_logical_error;
                local_errors_AG += res_astreag.is_logical_error;
                local_errors_PAG += PAG_err;
                // local_errors_clique += res_clique.islo
                if(PAG_err == false){
                    auto flipped_bits = qpd::syndrome_compressed(syndrome, false);
                    fp_t sw = predecoder->graph_sparcity_calculator(flipped_bits, sparcity_uw, distance);
                    if(sw == 1){
                        count_odd++;
                    }
                    vector_sparcity_weighted.push_back(sw);
                    vector_sparcity_unweighted.push_back(sparcity_uw);
                }
                for(uint stg = 0; stg <6; stg++){
                    local_stages[stg] += predecoder->reached_stage[stg];
                }

                if(predecoder->total_cycles > MAX_BF6_CYCLE){
                    local_number_above_threshold++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                if(predecoder->total_cycles_parallel_updating > MAX_BF6_CYCLE){
                    local_number_above_threshold_parallel++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                

                //Setting clock related varaibles
                if (local_max_cycle < predecoder->total_cycles){
                    local_max_cycle = predecoder->total_cycles;
                }
                if(local_max_cycle_parallel < predecoder->total_cycles_parallel_updating){
                    local_max_cycle_parallel = predecoder->total_cycles_parallel_updating;
                }
                local_sum_cycle += predecoder->total_cycles;
                local_sum_cycle_parallel += predecoder->total_cycles_parallel_updating;

                if(local_max_round < predecoder->number_of_rounds){
                    local_max_round = predecoder->number_of_rounds;
                }
                if(predecoder->MAX_BF_HW == MAX_BF10_HW && hw > MAX_BF10_HW){
                    local_n_Astrea10_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF8_HW && hw > MAX_BF10_HW){
                    local_n_Astrea8_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF6_HW && hw > MAX_BF10_HW){
                    local_n_Astrea6_used++;
                }
            }
            shots -= shots_this_round;
        }
        uint64_t fault_errors = 0;
        MPI_Allreduce(&local_errors, &fault_errors, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_PAG = 0;
        MPI_Allreduce(&local_errors_PAG, &fault_errors_PAG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_AG = 0;
        MPI_Allreduce(&local_errors_AG, &fault_errors_AG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);


        uint64_t times_lhw = 0;
        MPI_Allreduce(&local_times_lhw, &times_lhw, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t, 6> stages;
        stages.fill(0);
        MPI_Allreduce(&local_stages, &stages, 6, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle = 0;
        MPI_Allreduce(&local_max_cycle, &max_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle = 0;
        MPI_Allreduce(&local_sum_cycle, &sum_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle_parallel = 0;
        MPI_Allreduce(&local_max_cycle_parallel, &max_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle_parallel = 0;
        MPI_Allreduce(&local_sum_cycle_parallel, &sum_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold = 0;
        MPI_Allreduce(&local_number_above_threshold, &number_above_threshold, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold_parallel = 0;
        MPI_Allreduce(&local_number_above_threshold_parallel, &number_above_threshold_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_round = 0;
        MPI_Allreduce(&local_max_round, &max_round, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        
        uint64_t n_Astrea10_used = 0;
        MPI_Allreduce(&local_n_Astrea10_used, &n_Astrea10_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea8_used = 0;
        MPI_Allreduce(&local_n_Astrea8_used, &n_Astrea8_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea6_used = 0;
        MPI_Allreduce(&local_n_Astrea6_used, &n_Astrea6_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }

        if (fault_errors > 0 || ((number_above_threshold_parallel -  local_number_above_threshold_parallel_prior) > 0)) {
            fp_t failure_rate = ((fp_t)fault_errors/simulator.shots_per_k[k+min_k]);
        
            prev_logical_error_rate = statres.logical_error_rate;
            statres.logical_error_rate += (failure_rate*simulator.prob_k[k+min_k]);
            statres.n_logical_errors += fault_errors;
            fp_t time_failure_rate = (fp_t)(number_above_threshold_parallel -  local_number_above_threshold_parallel_prior)/simulator.shots_per_k[k+min_k];
            time_constrained_ler += (failure_rate*simulator.prob_k[k+min_k]) + (time_failure_rate*simulator.prob_k[k+min_k]);
            local_number_above_threshold_parallel_prior = number_above_threshold_parallel;

            if (world_rank == 0) {
                std::cout << "===== Promatch =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate << std::endl
                        << "\tNumber of errors = " << statres.n_logical_errors << std::endl;
                std::cout << "\tLogical error rate = " << statres.logical_error_rate << std::endl;
                std::cout << "\t\t Logical error rate timed = " << time_constrained_ler << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_PAG > 0) {
            fp_t failure_rate_PAG = ((fp_t)fault_errors_PAG/simulator.shots_per_k[k+min_k]);
        
            ler_PAG += (failure_rate_PAG*simulator.prob_k[k+min_k]);
            number_of_logical_err_PAG += fault_errors_PAG;
        
            if (world_rank == 0) {
                std::cout << "===== Promatch || AG =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_PAG << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_PAG << std::endl;
                std::cout << "\tLogical error rate = " << ler_PAG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_AG > 0) {
            fp_t failure_rate_AG = ((fp_t)fault_errors_AG/simulator.shots_per_k[k+min_k]);
        
            ler_AG += (failure_rate_AG*simulator.prob_k[k+min_k]);
            number_of_logical_err_AG += fault_errors_AG;
        
            if (world_rank == 0) {
                std::cout << "===== AstreaG =====" << std::endl;
                std::cout << "\tFailure rate = " << failure_rate_AG << std::endl
                        << "\tNumber of errors = " << number_of_logical_err_AG << std::endl;
                std::cout << "\tLogical error rate = " << ler_AG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        std::array<fp_t,6> stages_rate;
        stages_rate.fill(0);
        for(uint stg = 0; stg < 6; stg++){
            stages_rate[stg] = ((fp_t)stages[stg]/(simulator.shots_per_k[k+min_k]-times_lhw));
            normalized_stages[stg] += stages_rate[stg]*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
            weighted_stages[stg] += stages_rate[stg]*simulator.prob_k[k+min_k];
 
        }

        // std::cout << world_rank <<" sum_cycle"  <<sum_cycle<< " sum_cycle_parallel: " <<sum_cycle_parallel<<std::endl;

        final_average_cycle += ((fp_t)sum_cycle/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];        
        final_average_cycle_normalized += ((fp_t)sum_cycle/(simulator.shots_per_k[k+min_k]))*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        final_average_cycle_parallel += ((fp_t)sum_cycle_parallel/(simulator.shots_per_k[k+min_k]))*simulator.prob_k[k+min_k];
        final_average_cycle_parallel_normalized += ((fp_t)sum_cycle_parallel/(simulator.shots_per_k[k+min_k]))*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);

        fp_t percentage_astrea10 = ((fp_t)n_Astrea10_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea10_used +=percentage_astrea10*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
        
        fp_t percentage_astrea8 = ((fp_t)n_Astrea8_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea8_used +=percentage_astrea8*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);
        
        fp_t percentage_astrea6 = ((fp_t)n_Astrea6_used/(simulator.shots_per_k[k+min_k]-times_lhw));
        total_Astrea6_used +=percentage_astrea6*((fp_t)simulator.prob_k[k+min_k]/sum_normalized_prob);


        // std::cout << "world_rank: "<<world_rank  << " final_average_cycle: " << final_average_cycle <<
        // " final_average_cycle_normalized: " << final_average_cycle_normalized << " final_average_cycle_parallel: " << final_average_cycle_parallel << 
        // " final_average_cycle_parallel_normalized: " << final_average_cycle_parallel_normalized << std::endl;
        
        if (world_rank == 0) {
            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\t\tS1: " << stages_rate[0] << " S2_1: " << stages_rate[1] << " S2_2: " << stages_rate[2] << " S3: " << stages_rate[3] << " S4_1: " << stages_rate[4] 
                            << " S4_2: " << stages_rate[5] << std::endl;
            std::cout << "\tNormalized Stage Rates:" << std::endl;
            std::cout << "\t\tS1: " << normalized_stages[0] << " S2_1: " << normalized_stages[1] << " S2_2: " << normalized_stages[2] << " S3: " << normalized_stages[3] << " S4_1: " << normalized_stages[4] 
                            << " S4_2: " << normalized_stages[5] << std::endl;
            std::cout << "\tStage Rates:" << std::endl;
            std::cout << "\t\tS1: " << weighted_stages[0] << " S2_1: " << weighted_stages[1] << " S2_2: " << weighted_stages[2] << " S3: " << weighted_stages[3] << " S4_1: " << weighted_stages[4] 
                            << " S4_2: " << weighted_stages[5] << std::endl;

            std::cout << "\t -----------------"<< std::endl;
            std::cout << "\tNormalized:" << std::endl;
            std::cout << "\t\t Avg Cycles (UP) = " << final_average_cycle_parallel_normalized  << "(" << final_average_cycle_normalized <<")"<<std::endl;
            std::cout << "\tNot Normalized:" << std::endl;
            std::cout << "\t\t Avg Cycles(UP) = " << final_average_cycle_parallel << "(" << final_average_cycle <<")"<<std::endl;
            std::cout << "\tMAX Cycles (UP) = " << max_cycle_parallel  << "(" << max_cycle << ") Number above threshold(UP):" << number_above_threshold_parallel << "(" << number_above_threshold << ")"<< std::endl;
            std::cout << "\tMAX Round = " << max_round  << std::endl;
            std::cout << "\t ----------ASTREA PERCENTAGES-------"<< std::endl;
            std::cout << "\t This bucket" << std::endl;
            std::cout << "\t\tAstrea10: " << percentage_astrea10 << " Astrea8: " << percentage_astrea8 << " Astrea6: " << percentage_astrea6 << std::endl;
            std::cout << "\t Total(Normalized)" << std::endl;
            std::cout << "\t\tAstrea10: " << total_Astrea10_used << " Astrea8: " << total_Astrea8_used << " Astrea6: " << total_Astrea6_used << std::endl;

            std::cout << "_______________________________________________"<< std::endl;
            if (!vector_sparcity_unweighted.empty()) {
                std::cout << "========== Unweighted Sparcity ========== "  << std::endl;
                double sum = std::accumulate(vector_sparcity_unweighted.begin(), vector_sparcity_unweighted.end(), 0.0);
                double mean = sum / vector_sparcity_unweighted.size();

                // Calculate variance
                double sq_sum = std::inner_product(vector_sparcity_unweighted.begin(), vector_sparcity_unweighted.end(), vector_sparcity_unweighted.begin(), 0.0);
                double variance = sq_sum / vector_sparcity_unweighted.size() - std::pow(mean, 2);

                // Calculate 2.5th and 97.5th percentiles for the 95% percentile interval
                std::sort(vector_sparcity_unweighted.begin(), vector_sparcity_unweighted.end());
                int lowerIndex = static_cast<int>(std::ceil(0.025 * vector_sparcity_unweighted.size()) - 1);
                int upperIndex = static_cast<int>(std::ceil(0.975 * vector_sparcity_unweighted.size()) - 1);
                double percentile2_5 = vector_sparcity_unweighted[lowerIndex];
                double percentile97_5 = vector_sparcity_unweighted[upperIndex];

                // Output results
                std::cout << "Mean: " << mean << std::endl;
                std::cout << "Variance: " << variance << std::endl;
                std::cout << "95% Percentile Interval: [" << percentile2_5 << ", " << percentile97_5 << "]" << std::endl;
            } 

            if (!vector_sparcity_weighted.empty()) {
                std::cout << "========== Weighted Sparcity ========== "  << std::endl;
                double sum = std::accumulate(vector_sparcity_weighted.begin(), vector_sparcity_weighted.end(), 0.0);
                double mean = sum / vector_sparcity_weighted.size();

                // Calculate variance
                double sq_sum = std::inner_product(vector_sparcity_weighted.begin(), vector_sparcity_weighted.end(), vector_sparcity_weighted.begin(), 0.0);
                double variance = sq_sum / vector_sparcity_weighted.size() - std::pow(mean, 2);

                // Calculate 2.5th and 97.5th percentiles for the 95% percentile interval
                std::sort(vector_sparcity_weighted.begin(), vector_sparcity_weighted.end());
                int lowerIndex = static_cast<int>(std::ceil(0.025 * vector_sparcity_weighted.size()) - 1);
                int upperIndex = static_cast<int>(std::ceil(0.975 * vector_sparcity_weighted.size()) - 1);
                double percentile2_5 = vector_sparcity_weighted[lowerIndex];
                double percentile97_5 = vector_sparcity_weighted[upperIndex];

                // Output results
                std::cout << "Mean: " << mean << std::endl;
                std::cout << "Variance: " << variance << std::endl;
                std::cout << "95% Percentile Interval: [" << percentile2_5 << ", " << percentile97_5 << "]" << std::endl;
            } 

        }

    }

    std::cout << "***********" << count_odd << std::endl;

}

void print_mwpm_mistakes(uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k,  std::ofstream& out,  bool generalized,
 std::string decoder_name, bool write_file){
    std::string syndrome_file = "../mwpm_mistakes_d11_all.txt";//"../NFDecoder/data/ChallengingSyndromes/Rtests/P11/promatch_mistakes_d11_all.txt";//All_PromatchAM_d13.txt";//A13/All_AstreaG_d13.txt" challengingSyndromes_d13/All_AstreaG_d13_final.txt";//Promatch_A_d13_all_Multi_Astrea.txt";////Promatch_A_d11_all_Multi_Astrea.txt";//final_promatchA_mistakes.txt";//ler_miscorrection_d13_ASTREAG_r1_b1416.txt";//All_AstreaG_d13.txt";//all_promatch_u1_d13.txt";//ler_miscorrection_d13_ASTREAG_r1_b1416.txt";//PromatchA_d13_r0_all.txt";//PromatchA_d13_r0.txt";;//P//promatchA_first_trial.txt";//promatchA_first_trial.txt";//;

    // std::cout << "Enter a the path to syndromes: ";
    // std::getline(std::cin, syndrome_file);

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
    if(distance == 13){
        max_k = 17;
    } 
    else if(distance == 11){
        max_k = 12;
    }

    qrc::Decoder* decoder;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ, generalized);
    predecoder->set_decoder(decoder);
 

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(0.01*0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    

    std::vector<std::vector<uint16_t>> vecss = qpd::read_vector_of_vectors(syndrome_file);
    std::cout << "Size of incorrect pairs = " << vecss.size() << " - generalized = "<< predecoder->general<<std::endl;


    qpd::CountResult comp_res;
    for(auto compressed_synd : vecss){
        std::vector<uint8_t> syndrome = qpd::syndrome_decompressed(compressed_synd, distance);        
        auto res_mwpm = decoder->decode_error(syndrome);
        //qrc::DecoderShotResult res_astreag;
        if(res_mwpm.is_logical_error){
            for(auto det: compressed_synd){
                std::cout << det << " ";
                // if(det >= predecoder->n_detectors){
                //     std::cout << " - " << det;
                // }
            }
            std::cout <<std::endl;

        }
    }
}


void testing_all(uint64_t max_shot, uint distance, 
fp_t physcial_error, fp_t meas_er, uint64_t min_k, uint64_t max_k, fp_t threshold_scale, bool& print_time, std::string decoder_name, bool generalized, uint round_n, bool save_syndromes, std::string syndrome_folder_name, std::string syndrome_clock_challenging_file, uint64_t hshots_replc, uint64_t lshots_rplac, bool only_important_bucket){
    uint64_t SHOTS_PER_BATCH = 1'000'000;
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if(only_important_bucket){
        if(distance == 11){
            min_k = 7;
            max_k = 10;
        }
        else if(distance == 13){
            min_k = 9;
            max_k = 12;
        }
    }

    // std::string syndrome_file = "../NFDecoder/data/challengingSyndromes_d13/challenging_d13.txt";
    if(world_rank == 0){
        std::cout << "round#" <<round_n <<" Promatch-A + " << decoder_name <<std::endl;
    }



    // for saving syndromes
    uint64_t round = 1;
    uint64_t round_t = 1;
    if (world_rank == 0 && save_syndromes){
        std::cout << syndrome_folder_name<< std::endl;
        // To check if we are reading correct files
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

    
    bbsim::BucketBasedSim simulator(circ,max_k);
    qrc::Decoder* decoder;
    qrc::benchmark::StatisticalResult statres;

    if (decoder_name.find("MWPM") != std::string::npos){
        decoder = new qrc::MWPMDecoder(circ);
    }
    else{
        std::cout << "NO SPECIFIED DECODER";
        return;
    }
    std::mt19937_64 rng(world_size*round_n+world_rank);
    qpd::Predecoder* predecoder = new qpd::Predecoder(circ, generalized);
    predecoder->set_decoder(decoder);

    //////
    std::vector<vector<vector<uint16_t>>> saved_hhw_syndromes;
    saved_hhw_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits;
    //////
    /// For saving time challenging syndromes
    std::vector<vector<vector<uint16_t>>> saved_timechallenging_syndromes;
    saved_timechallenging_syndromes.resize(world_size);
    std::vector<uint16_t> flipped_bits_t;
    ////



    ///// ================== For parallel running with AstreaG ========================

    qrc::AstreaParams astreaG_param= {};
    qrc::Decoder* AstreaG_decoder;
    astreaG_param.bfu_fetch_width = 2;
    astreaG_param.bfu_priority_queue_size = 8;
    astreaG_param.main_clock_frequency = 250e6;
    astreaG_param.bfu_compute_stages = 2;
    astreaG_param.n_registers = 2000;
    astreaG_param.use_mld = false;

    uint n_detectors_per_round = (distance*distance - 1)/2;
    uint32_t weight_filter_cutoff =  -MWPM_INTEGER_SCALE*log10(threshold_scale*0.03*pow((physcial_error/0.0054), ((distance+1)/2))); // -1*log10(0.01*1e-9);
    if(world_rank == 0){
        std::cout << "threshold_scale:" << threshold_scale << std::endl;
    }


    AstreaG_decoder =  new  qrc::Astrea(circ,
                                        n_detectors_per_round,
                                        weight_filter_cutoff,
                                        astreaG_param);


    ///// ==============================================================================

    fp_t expected_ler = 0.03*pow((physcial_error/0.0054), ((distance+1)/2));
    simulator.set_shots_per_k(expected_ler, max_shot,true, hshots_replc, lshots_rplac);

    uint64_t max_k_ = max_k;
    fp_t sum_normalized_prob = 0;
    for(uint k = 0; k < max_k-min_k; k++){
        sum_normalized_prob += simulator.prob_k[k+min_k];

        // To eliminate the buckets that their probability is lower than LER and accordinly has zero number of shots
        if(simulator.shots_per_k[k+min_k] == 0){
            max_k_ --;
        }
    }
    max_k = max_k_;

    if(world_rank == 0){
        print_time = true;
        std::cout << "Distance = " << distance << " Expected LER = " << expected_ler << " Max Shots per bucket = " << max_shot << " Physical Err = "<< physcial_error<<std::endl;
        for(uint x=0;x<max_k-min_k; x++){
            std::cout << "k = " << x+min_k<< "( " << simulator.prob_k[x+min_k] << ") -itr= " << simulator.shots_per_k[x+min_k]<<std::endl;
        }
    }
    fp_t prev_logical_error_rate = 0.0;

    
    std::array<fp_t,6> normalized_stages;
    normalized_stages.fill(0);

    std::array<fp_t,6> weighted_stages;
    weighted_stages.fill(0);

    ///////Clock Simulation
    uint64_t local_sum_cycle = 0;

    uint64_t local_max_cycle = 0;

    uint64_t local_max_round = 0;

    uint64_t local_sum_cycle_parallel = 0;

    uint64_t local_max_cycle_parallel = 0; 


    fp_t final_average_cycle = 0;
    fp_t final_average_cycle_normalized = 0;
    fp_t final_average_cycle_parallel = 0;
    fp_t final_average_cycle_parallel_normalized = 0;

    uint64_t local_number_above_threshold_parallel = 0;
    uint64_t local_number_above_threshold = 0;

    uint64_t local_n_Astrea10_used = 0;
    uint64_t local_n_Astrea8_used = 0;
    uint64_t local_n_Astrea6_used = 0;

    fp_t total_Astrea10_used = 0;
    fp_t total_Astrea8_used = 0;
    fp_t total_Astrea6_used = 0;

    fp_t time_constrained_ler = 0;

    fp_t ler_PAG = 0;
    fp_t ler_clique = 0;
    fp_t ler_smith = 0;
    fp_t ler_clique_AG = 0;
    fp_t ler_smith_AG = 0;
    fp_t ler_AG = 0;

    uint64_t number_of_logical_err_PAG = 0;
    uint64_t number_of_logical_err_clique = 0;
    uint64_t number_of_logical_err_smith = 0;
    uint64_t number_of_logical_err_clique_AG = 0;
    uint64_t number_of_logical_err_smith_AG = 0;
    uint64_t number_of_logical_err_AG = 0;
    uint64_t local_number_above_threshold_parallel_prior = 0;

    for(uint k = 0; k < max_k-min_k; k++){
        uint64_t local_errors = 0;
        uint64_t local_errors_smith = 0;
        uint64_t local_errors_clique = 0;
        uint64_t local_errors_PAG = 0;
        uint64_t local_errors_smith_AG = 0;
        uint64_t local_errors_clique_AG = 0;
        uint64_t local_errors_AG = 0;
        
        std::array<uint64_t, 6> local_stages;
        local_stages.fill(0);

        uint64_t local_times_lhw = 0;

        local_sum_cycle_parallel = 0;
        local_sum_cycle = 0;

        local_n_Astrea10_used = 0;
        local_n_Astrea8_used = 0;
        local_n_Astrea6_used = 0;


        uint64_t shots = simulator.shots_per_k[k+min_k] / world_size;
        if (world_rank == world_size - 1) {
            shots += simulator.shots_per_k[k+min_k] % world_size;
        }
        while(shots > 0){
            
            uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
            
            auto buffer = simulator.create_syndromes_simd( k+min_k, shots_this_round, rng, false);
            for(uint s = 0; s<shots_this_round; s++){
                std::vector<uint8_t> syndrome = qrc::_to_vector(buffer[s], simulator.n_detectors, simulator.n_observables);
                uint hw = std::accumulate(syndrome.begin(), syndrome.begin()+predecoder->n_detectors,0);
                if(hw <= MAX_BF10_HW){
                    local_times_lhw++;
                }

                qrc::DecoderShotResult res;
                qpd::PredecoderDecoderShotResult res_smith;
                qpd::PredecoderDecoderShotResult res_clique;

                qrc::DecoderShotResult res_astreag;


                res = predecoder->adaptively_decode_error(syndrome);
                res_smith = predecoder->smith_decode_error(syndrome);
                res_clique = predecoder->clique_decode_error(syndrome);

                res_astreag = AstreaG_decoder->decode_error(syndrome);
                

                fp_t astreag_weight = predecoder->calc_matching_weight(res_astreag.matching);
                fp_t promatch_weight = predecoder->calc_matching_weight(res.matching);
                local_errors+= res.is_logical_error;
                bool PAG_err = (astreag_weight < promatch_weight)? res_astreag.is_logical_error :res.is_logical_error;
                local_errors_AG += res_astreag.is_logical_error;
                local_errors_PAG += PAG_err;
                local_errors_clique += res_clique.is_logical_error;
                local_errors_smith += res_smith.is_logical_error;
                local_errors_clique_AG += (astreag_weight < res_clique.weight)? res_astreag.is_logical_error :res_clique.is_logical_error;
                local_errors_smith_AG += (astreag_weight < res_smith.weight)? res_astreag.is_logical_error :res_smith.is_logical_error;
                // local_errors_clique += res_clique.islo
                for(uint stg = 0; stg <6; stg++){
                    local_stages[stg] += predecoder->reached_stage[stg];
                }
                
                if(res.is_logical_error){
                    //std::cout << "Error" << k <<"-";
                    flipped_bits = qpd::syndrome_compressed(syndrome, false);
                }
                if(PAG_err){
                    std::cout << "* ";
                    flipped_bits = qpd::syndrome_compressed(syndrome, false);
                }

                if(predecoder->total_cycles > MAX_BF6_CYCLE){
                    local_number_above_threshold++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                if(predecoder->total_cycles_parallel_updating > MAX_BF6_CYCLE){
                    local_number_above_threshold_parallel++;
                    flipped_bits_t = qpd::syndrome_compressed(syndrome, false);
                    saved_timechallenging_syndromes[world_rank].push_back(flipped_bits_t);
                }
                

                //Setting clock related varaibles
                if (local_max_cycle < predecoder->total_cycles){
                    local_max_cycle = predecoder->total_cycles;
                }
                if(local_max_cycle_parallel < predecoder->total_cycles_parallel_updating){
                    local_max_cycle_parallel = predecoder->total_cycles_parallel_updating;
                }
                local_sum_cycle += predecoder->total_cycles;
                local_sum_cycle_parallel += predecoder->total_cycles_parallel_updating;

                if(local_max_round < predecoder->number_of_rounds){
                    local_max_round = predecoder->number_of_rounds;
                }
                if(predecoder->MAX_BF_HW == MAX_BF10_HW && hw > MAX_BF10_HW){
                    local_n_Astrea10_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF8_HW && hw > MAX_BF10_HW){
                    local_n_Astrea8_used++;
                }
                else if(predecoder->MAX_BF_HW == MAX_BF6_HW && hw > MAX_BF10_HW){
                    local_n_Astrea6_used++;
                }
            }
            shots -= shots_this_round;
        }
        uint64_t fault_errors = 0;
        MPI_Allreduce(&local_errors, &fault_errors, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_smith = 0;
        MPI_Allreduce(&local_errors_smith, &fault_errors_smith, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_clique = 0;
        MPI_Allreduce(&local_errors_clique, &fault_errors_clique, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_PAG = 0;
        MPI_Allreduce(&local_errors_PAG, &fault_errors_PAG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_smith_AG = 0;
        MPI_Allreduce(&local_errors_smith_AG, &fault_errors_smith_AG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_clique_AG = 0;
        MPI_Allreduce(&local_errors_clique_AG, &fault_errors_clique_AG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t fault_errors_AG = 0;
        MPI_Allreduce(&local_errors_AG, &fault_errors_AG, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);


        uint64_t times_lhw = 0;
        MPI_Allreduce(&local_times_lhw, &times_lhw, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        std::array<uint64_t, 6> stages;
        stages.fill(0);
        MPI_Allreduce(&local_stages, &stages, 6, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle = 0;
        MPI_Allreduce(&local_max_cycle, &max_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle = 0;
        MPI_Allreduce(&local_sum_cycle, &sum_cycle, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_cycle_parallel = 0;
        MPI_Allreduce(&local_max_cycle_parallel, &max_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        uint64_t sum_cycle_parallel = 0;
        MPI_Allreduce(&local_sum_cycle_parallel, &sum_cycle_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold = 0;
        MPI_Allreduce(&local_number_above_threshold, &number_above_threshold, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t number_above_threshold_parallel = 0;
        MPI_Allreduce(&local_number_above_threshold_parallel, &number_above_threshold_parallel, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        uint64_t max_round = 0;
        MPI_Allreduce(&local_max_round, &max_round, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);
        
        uint64_t n_Astrea10_used = 0;
        MPI_Allreduce(&local_n_Astrea10_used, &n_Astrea10_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea8_used = 0;
        MPI_Allreduce(&local_n_Astrea8_used, &n_Astrea8_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
        uint64_t n_Astrea6_used = 0;
        MPI_Allreduce(&local_n_Astrea6_used, &n_Astrea6_used, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

        if (world_rank == 0 && simulator.shots_per_k[k+min_k]!=0) {
            std::cout << "Bucket = " << k+min_k << "\n"
                    << "\tProbability of bucket = " << simulator.prob_k[k+min_k] << std::endl;
        }

        if (fault_errors > 0 || ((number_above_threshold_parallel -  local_number_above_threshold_parallel_prior) > 0)) {
            fp_t failure_rate = ((fp_t)fault_errors/simulator.shots_per_k[k+min_k]);
        
            prev_logical_error_rate = statres.logical_error_rate;
            statres.logical_error_rate += (failure_rate*simulator.prob_k[k+min_k]);
            statres.n_logical_errors += fault_errors;
            fp_t time_failure_rate = (fp_t)(number_above_threshold_parallel -  local_number_above_threshold_parallel_prior)/simulator.shots_per_k[k+min_k];
            time_constrained_ler += (failure_rate*simulator.prob_k[k+min_k]) + (time_failure_rate*simulator.prob_k[k+min_k]);
            local_number_above_threshold_parallel_prior = number_above_threshold_parallel;

            if (world_rank == 0) {
                std::cout << "===== Promatch =====" << std::endl;
                std::cout << "\tLogical error rate = " << statres.logical_error_rate << std::endl;
                std::cout << "\t\t Logical error rate timed = " << time_constrained_ler << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_PAG > 0) {
            fp_t failure_rate_PAG = ((fp_t)fault_errors_PAG/simulator.shots_per_k[k+min_k]);
        
            ler_PAG += (failure_rate_PAG*simulator.prob_k[k+min_k]);
            number_of_logical_err_PAG += fault_errors_PAG;
        
            if (world_rank == 0) {
                std::cout << "===== Promatch || AG =====" << std::endl;
                std::cout << "\tLogical error rate = " << ler_PAG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_AG > 0) {
            fp_t failure_rate_AG = ((fp_t)fault_errors_AG/simulator.shots_per_k[k+min_k]);
        
            ler_AG += (failure_rate_AG*simulator.prob_k[k+min_k]);
            number_of_logical_err_AG += fault_errors_AG;
        
            if (world_rank == 0) {
                std::cout << "===== AstreaG =====" << std::endl;
                std::cout << "\tLogical error rate = " << ler_AG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_clique > 0) {
            fp_t failure_rate_clique = ((fp_t)fault_errors_clique/simulator.shots_per_k[k+min_k]);
        
            ler_clique += (failure_rate_clique*simulator.prob_k[k+min_k]);
            number_of_logical_err_clique += fault_errors_clique;
        
            if (world_rank == 0) {
                std::cout << "===== Clique =====" << std::endl;
                std::cout << "\tLogical error rate = " << ler_clique << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        // if (fault_errors_clique_AG > 0) {
        //     fp_t failure_rate_clique_AG = ((fp_t)fault_errors_clique_AG/simulator.shots_per_k[k+min_k]);
        
        //     ler_clique_AG += (failure_rate_clique_AG*simulator.prob_k[k+min_k]);
        //     number_of_logical_err_clique_AG += fault_errors_clique_AG;
        
        //     if (world_rank == 0) {
        //         std::cout << "===== Clique || AG =====" << std::endl;
        //         std::cout << "\tFailure rate = " << failure_rate_clique_AG << std::endl
        //                 << "\tNumber of errors = " << number_of_logical_err_clique_AG << std::endl;
        //         std::cout << "\tLogical error rate = " << ler_clique_AG << std::endl;
        //         std::cout << "====================" << std::endl;
        //     }
        // }
        if (fault_errors_smith > 0) {
            fp_t failure_rate_smith = ((fp_t)fault_errors_smith/simulator.shots_per_k[k+min_k]);
        
            ler_smith += (failure_rate_smith*simulator.prob_k[k+min_k]);
            number_of_logical_err_smith += fault_errors_smith;
        
            if (world_rank == 0) {
                std::cout << "===== Smith =====" << std::endl;
                std::cout << "\tLogical error rate = " << ler_smith << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
        if (fault_errors_smith_AG > 0) {
            fp_t failure_rate_smith_AG = ((fp_t)fault_errors_smith_AG/simulator.shots_per_k[k+min_k]);
        
            ler_smith_AG += (failure_rate_smith_AG*simulator.prob_k[k+min_k]);
            number_of_logical_err_smith_AG += fault_errors_smith_AG;
        
            if (world_rank == 0) {
                std::cout << "===== Smith || AG =====" << std::endl;
                std::cout << "\tLogical error rate = " << ler_smith_AG << std::endl;
                std::cout << "====================" << std::endl;
            }
        }
    }

}

