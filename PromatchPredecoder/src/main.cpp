#include "experiments.h"
#include "predecoder.h"
#include "syndrome_surgery.h"
#include "bb_benchmark.h"

 #include <iostream>
#include <vector>
#include <iostream>
#include <vector>
#include <chrono>

using namespace std;

// vector<double> calculateCoefficients(vector<double> p, int k) {
//     int n = p.size();
//     vector<double> coef(k+1, 0);
//     coef[0] = 1;
    
//     for(int j=0; j<n; j++) {
//         double pj = p[j];
//         for(int i=k; i>=1; i--) {
//             coef[i] = coef[i-1] * pj + coef[i] * (1.0 - pj);
//         }
//         coef[0] *= (1.0 - pj);
//     }
//     return coef;
// }

// int main() {
//     vector<double> p = {((double)1/2), ((double)2/3), ((double)3/4)};
//     int k = 2;
    
//     vector<double> coef = calculateCoefficients(p, k);
    
//     for(int i=0; i<=k; i++) {
//         cout << "Coefficient of x^" << i << " in P(x) is " << coef[i] << endl;
//     }
    
//     return 0;
// }

int main(int argc, char** argv){

    

    bool flag_s = false;
    bool print_logs = false;
    bool print_time = false;
    bool save_syndrome = true;
    bool generalized = false;
    bool important_bucket = false;
    int s = 0;
    int round_n = 1;
    uint64_t max_shots =10'000'000;//1'000'000'000;
    uint64_t high_shots_replacement = 2'000'000;
    uint64_t low_shots_replacement = 1'000'000;
    uint distance = 15;
    fp_t PER_const = 1;

    std::string experiment_type;//";//_notBoundry";
    std::string decoder_name = "MWPM";//;"ASTREAG"
    std::string predecoder_name = "CLIQUE";
    fp_t threshold_scale = 1;

    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if(arg == "-e" && i + 1 < argc){
            experiment_type = "exp" + std::string(argv[i + 1]);
        }
        if(arg == "-d" && i + 1 < argc){
            distance = std::atoi(argv[i + 1]);
        }
        if (arg == "-s" && i + 1 < argc) {
            s = std::stoi(argv[i + 1]);
            if(s == 1){
                flag_s = true; // set the flag_s variable to true
            }
        }
        if (arg == "-r" && i + 1 < argc) {
            round_n = std::atoi(argv[i + 1]);
        }
        if(arg == "-ms" && i + 1 < argc){
            max_shots = parseMagnitude(argv[i + 1]);

        }
        if(arg == "-hs" && i + 1 < argc){
            high_shots_replacement = parseMagnitude(argv[i + 1]);

        }
        if(arg == "-ls" && i + 1 < argc){
            low_shots_replacement = parseMagnitude(argv[i + 1]);

        }
        if(arg == "-p" && i + 1 < argc){
            PER_const = (fp_t)atoi(argv[i + 1]);

        }
        if(arg == "-t" && i + 1 < argc){
            decoder_name = argv[i + 1];

        }
        if(arg == "-pt" && i + 1 < argc){
            predecoder_name = argv[i + 1];

        }
        if(arg == "-g" && i + 1 < argc){
            generalized = std::atoi(argv[i + 1]);

        }
        if(arg == "-ib" && i + 1 < argc){
            important_bucket = std::atoi(argv[i + 1]);

        }
        if(arg == "-th" && i + 1 < argc){
            threshold_scale = std::atof(argv[i + 1]);
        }
    }


    fp_t physical_er = PER_const*1e-4; 
    fp_t m_error = PER_const*1e-4;

    MPI_Init(NULL, NULL);
    
    
    // "exp1" : experiment type 1: HW vs different error chain length in MWPM decoder
    // "exp2" : experiment type 2: HW | # degree one nodes in the decoding graph
    // "exp3" : experiment type 3: Same as 2 but with distance 9 instead of 7
    //         HW | # degree one nodes in the decoding graph
    bool mention_shots = true;
    std::string file_name = " ";
    std::string matching_file = " ";
    std::string data_folder_name = " ";
    uint64_t min_k = std::floor(distance/2);// < 6 ? 6: std::floor(distance/2);
    uint64_t max_k = 35;
    if(flag_s){
        data_folder_name = "../NFDecoder/s_d/";
    }
    else
        data_folder_name = "../NFDecoder/data/";

    
    std::string syndrome_folder_name = data_folder_name  + "syndromes/" + number_of_shots_to_string(max_shots) +"/";
    // std::vector<std::vector<uint16_t>> vecss = qpd::read_vector_of_vectors(syndrome_folder_name + "sgen_p=1.0e-3_m=1.0e-3_d=11_s=0M_c=4_r=1.txt");

    // for( auto vec : vecss){
    //     for( auto x : vec){
    //         std::cout << x << ",";
    //     }
    //     std::cout << std::endl;
    // }
    // std::vector<uint8_t> s = qpd::syndrome_decompressed(vecss[vecss.size()-1], distance);

    // qpd::print_syndrome(s);

    // return 0;
    

    file_name = data_folder_name + generate_experiment_name(experiment_type, distance, physical_er, 
                                                m_error, max_shots, mention_shots);

    bool add_boundry = true;
    std::ofstream outputFile;
    if( experiment_type.find("exp16") == std::string::npos and 
        experiment_type.find("exp19") == std::string::npos and 
        experiment_type.find("exp22") == std::string::npos and 
        experiment_type.find("exp21") == std::string::npos and 
        experiment_type.find("exp23") == std::string::npos and 
        experiment_type.find("exp24") == std::string::npos and 
        experiment_type.find("exp25") == std::string::npos and 
        experiment_type.find("exp26") == std::string::npos and 
        experiment_type.find("exp27") == std::string::npos and 
        experiment_type.find("exp28") == std::string::npos and 
        experiment_type.find("exp29") == std::string::npos and 
        experiment_type.find("exp30") == std::string::npos and 
        experiment_type.find("exp31") == std::string::npos and 
        experiment_type.find("exp33") == std::string::npos and 
        experiment_type.find("exp34") == std::string::npos and 
        experiment_type.find("exp35") == std::string::npos and 
        experiment_type.find("exp36") == std::string::npos and 
        experiment_type.find("exp37") == std::string::npos and 
        experiment_type.find("exp38") == std::string::npos  and 
        experiment_type.find("exp39") == std::string::npos)
        outputFile  = std::ofstream(file_name);

    if (experiment_type.find("exp39") != std::string::npos){
        testing_all(max_shots, distance, physical_er,
            m_error, min_k, max_k, threshold_scale, print_time, "MWPM", generalized, round_n, false , syndrome_folder_name, " ", high_shots_replacement, low_shots_replacement, important_bucket);
    }
    else if (experiment_type.find("exp38") != std::string::npos){
        
        print_mwpm_mistakes( distance, physical_er, m_error, min_k, max_k, outputFile, generalized);
        std::cout << "X";

    }
    if (experiment_type.find("exp37") != std::string::npos){
        PAG_synergies(max_shots, distance, physical_er,
            m_error, min_k, max_k, print_time, "MWPM", generalized, round_n, false , syndrome_folder_name, " ", high_shots_replacement, low_shots_replacement);
    }
    else if (experiment_type.find("exp36") != std::string::npos){
        ler_predecoders(max_shots, distance, physical_er,
            m_error, min_k, max_k, threshold_scale, print_time, "MWPM", generalized, round_n, false , syndrome_folder_name, " ", high_shots_replacement, low_shots_replacement);
    }
    else if (experiment_type.find("exp35") != std::string::npos){
        HW_distributions(outputFile, max_shots, distance, physical_er,m_error, min_k, max_k,
                            print_time, high_shots_replacement, low_shots_replacement);
    }
    else if (experiment_type.find("exp34") != std::string::npos){
        bb_predecoder_ler_calc(max_shots, distance, physical_er, m_error, min_k, max_k, print_time, round_n, save_syndrome, syndrome_folder_name,
                                 decoder_name, predecoder_name, high_shots_replacement, low_shots_replacement);

        //bb_ler_calculation(decoder_name, outputFile, max_shots, distance, physical_er, m_error, min_k, max_k, 5'000'000);
    }
    else if (experiment_type.find("exp33") != std::string::npos){
        edge_distributions(outputFile, max_shots, distance, physical_er,m_error, min_k, max_k,
                            print_time, high_shots_replacement, low_shots_replacement);
    }
    else if(experiment_type.find("exp32") != std::string::npos){
        degree_distribution(outputFile, max_shots, distance, physical_er,m_error, min_k, max_k,
                            print_time, high_shots_replacement, low_shots_replacement);
    }
    else if(experiment_type.find("exp31") != std::string::npos){
        sample_test_degree_priority(max_shots, distance, physical_er, m_error, min_k, max_k, print_time);
    }
    else if(experiment_type.find("exp30") != std::string::npos){
        stats_of_incorrect_syndromes( distance, physical_er, m_error, min_k, max_k, generalized);
    }
    else if(experiment_type.find("exp29") != std::string::npos){
        bb_avg_ec1(max_shots, distance, physical_er, m_error, min_k, max_k, print_time);
    }
    else if(experiment_type.find("exp28") != std::string::npos){
        bool out = true;
        if(out){
            file_name = data_folder_name + "revision_promatch_mistakes_d" + std::to_string(distance) + ".txt";
            outputFile  = std::ofstream(file_name);
            check_the_incorrect_syndrome( distance, physical_er, m_error, min_k, max_k, outputFile, generalized, "MWPM", out);
        }
        else{
            check_the_incorrect_syndrome( distance, physical_er, m_error, min_k, max_k, outputFile, generalized);
        }
        
    }
    else if(experiment_type.find("exp27") != std::string::npos){
        std::string syndrome_folder_name_clock;
        if(save_syndrome){
            syndrome_folder_name = data_folder_name  + "Promatch_A_d" + std::to_string(distance)  
                                                    + "_r" + std::to_string(round_n) + "_Multi_Astrea" +"/";
            syndrome_folder_name_clock = data_folder_name  + "Promatch_A_timing_d" + std::to_string(distance)  
                                                    + "_r" + std::to_string(round_n) +"__Multi_Astrea" +"/";
        }
        bb_ler_calc_promatch(max_shots, distance, physical_er,
            m_error, min_k, max_k, print_time, "MWPM", generalized, round_n, save_syndrome, syndrome_folder_name, syndrome_folder_name_clock, high_shots_replacement, low_shots_replacement);
    }
    else if(experiment_type.find("exp26") != std::string::npos){
        test_size_of_predecoding(max_shots, distance, physical_er, m_error, min_k, max_k);

    }
    else if(experiment_type.find("exp25") != std::string::npos){
        test_priority_groups(max_shots, distance, physical_er, m_error, min_k, max_k, print_time);

    }
    else if(experiment_type.find("exp24") != std::string::npos){
            printing_incorrect_predictions(max_shots, distance, physical_er,
            m_error, min_k, max_k);
        }
    else if(experiment_type.find("exp23") != std::string::npos){
        test_finding_fast_group(max_shots, distance, physical_er,
        m_error, min_k, max_k);
    }
    else if(experiment_type.find("exp22") != std::string::npos){
         printing_samples_of_syndromes_and_matchings(max_shots, distance, 
            physical_er, m_error, min_k, max_k, true);
    }
    else if(experiment_type.find("exp21") != std::string::npos){
        test_b_statistical_ler_AstreaG(max_shots, distance, physical_er, m_error);
    }
    else if(experiment_type.find("exp20") != std::string::npos){
        if(save_syndrome){
            syndrome_folder_name = data_folder_name  + "ler_miscorrection_2_d" + std::to_string(distance)  
                            +"_"+ decoder_name + "_r" + std::to_string(round_n) +"/";
        }
        bb_decoder_ler_calc(max_shots, distance, physical_er, m_error, min_k, max_k, print_time, round_n, save_syndrome, syndrome_folder_name, decoder_name, high_shots_replacement, low_shots_replacement);

        //bb_ler_calculation(decoder_name, outputFile, max_shots, distance, physical_er, m_error, min_k, max_k, 5'000'000);
    }
    else if(experiment_type.find("exp19") != std::string::npos){
        test_b_statistical_ler(max_shots, distance, physical_er, m_error);
    }
    else if(experiment_type.find("exp18") != std::string::npos){
        file_name = data_folder_name + generate_experiment_name(experiment_type, distance, physical_er, 
                                                m_error, max_shots, false);

        if(experiment_type.find("notBoundry") != std::string::npos)
            add_boundry = false;
        bb_ler_calculation_AstreaG(outputFile, syndrome_folder_name, distance, physical_er, m_error, print_logs, add_boundry);
    }
    else if(experiment_type.find("exp17") != std::string::npos){

        if(experiment_type.find("notBoundry") != std::string::npos)
            add_boundry = false;
        bb_ler_calculation_and_stats(outputFile, syndrome_folder_name, max_shots, distance, physical_er, m_error, print_logs, add_boundry);
    }
    else if(experiment_type.find("exp16") != std::string::npos){
        if(experiment_type.find("notBoundry") != std::string::npos)
            add_boundry = false;
        test_sim(outputFile, syndrome_folder_name, max_shots, distance, physical_er, m_error, print_logs, add_boundry);
    }
    else if(experiment_type.find("exp15") != std::string::npos){
        if(experiment_type.find("notBoundry") != std::string::npos)
            add_boundry = false;
        predecoding_hhw_syndromes_astreaG(outputFile, syndrome_folder_name, max_shots, distance, physical_er, m_error, print_logs, add_boundry);
    }
    else if(experiment_type.find("exp14") != std::string::npos){
        if(experiment_type.find("notBoundry") != std::string::npos)
            add_boundry = false;
        predecoding_hhw_syndromes(outputFile, syndrome_folder_name, max_shots, distance, physical_er, m_error, print_logs, add_boundry);
    }
    else if(experiment_type.find("exp12") != std::string::npos){
            saving_hhw_syndromes(outputFile, max_shots,distance,physical_er,m_error, syndrome_folder_name, print_logs);
    }
    else if(experiment_type.find("exp13") != std::string::npos){
        HW_distr(outputFile,  max_shots, distance,  physical_er, m_error);
    }
    else if(experiment_type.find("exp11") != std::string::npos){
        check_subgraphs(distance, physical_er, m_error, 
        1, true);
    }
    else if (experiment_type.find("exp10") != std::string::npos){
        matching_file = "../NFDecoder/data/" + generate_matching_file_name(experiment_type, distance, physical_er, 
                            m_error, max_shots, mention_shots);
        std::ofstream MoutputFile(matching_file);
        high_distance_predecoder_plus_MWPM_decoder(outputFile, max_shots, distance, 
        physical_er, m_error, MoutputFile);
    }
    else if (experiment_type.find("exp1") != std::string::npos)
        error_chain_distr_v2(outputFile,  max_shots, distance,  physical_er, m_error);
    //check_sydromes_hw(max_shots, distance, p);
    else if (experiment_type.find("exp2") != std::string::npos)
        degree_one_distr(outputFile,  max_shots, distance,  physical_er, m_error);
    else if (experiment_type.find("exp3") != std::string::npos){
            distance = 9;
            degree_one_distr(outputFile,  max_shots, distance,  physical_er, m_error);
        }
    else if (experiment_type.find("exp4") != std::string::npos){
        predecoder_plus_MWPM_decoder(outputFile,  max_shots, distance,  physical_er, m_error);
    }
    else if (experiment_type.find("exp5") != std::string::npos)
        HW_remained_distr(outputFile,  max_shots, distance,  physical_er, m_error);
    else if (experiment_type.find("exp6") != std::string::npos)
        HW_remained_distr(outputFile,  max_shots, distance,  physical_er, m_error, 2);
    // Close the file if it was opened successfully
    // if (experiment_type.find("sycamore") == std::string::npos){
    //     if (f != nullptr) {
    //         std::fclose(f);
    //     }
    // }
    else if (experiment_type.find("exp7") != std::string::npos)
        decoding_cost_experiments(outputFile,  max_shots, distance,  physical_er, m_error);
    else if (experiment_type.find("exp8") != std::string::npos)
        HW_remained_distr_v2(outputFile,  max_shots, distance,  physical_er, m_error);
    else if (experiment_type.find("exp9") != std::string::npos){
        high_distance_experiments(outputFile,  max_shots, distance,  physical_er, m_error);
    }

    MPI_Finalize();

    auto end = std::chrono::high_resolution_clock::now();

    // Calculate the duration in hours, minutes, and seconds
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    int hours = duration / 3600;
    int minutes = (duration % 3600) / 60;
    int seconds = duration % 60;

    // Print the runtime in hours, minutes, and seconds from rank 0 process
    if(print_time){
        std::cout << std::endl << "Runtime: " << hours << " hours, " << minutes << " minutes, " << seconds << " seconds" << std::endl;
    }

    return 0;
}


// // #include <mpi.h>
// // #include <stdio.h>
// // #include <iostream>
// // // int main(int argc, char** argv) {
// // //     int version, subversion;
// // //     MPI_Init(&argc, &argv);
// // //     MPI_Get_version(&version, &subversion);
// // //     printf("MPI version: %d.%d\n", version, subversion);
// // //     MPI_Finalize();
// // //     return 0;
// // // }
// // int main(int argc, char** argv) {
// //     MPI_Init(NULL, NULL);

// //     int world_size, world_rank;
// //     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
// //     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
// //     std::cout << world_rank << world_size << std::endl;


// //     MPI_Finalize();
// //     return 0;
// // }

