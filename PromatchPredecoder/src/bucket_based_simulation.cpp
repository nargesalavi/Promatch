#include "bucket_based_simulation.h"
#include "error_event.h"

namespace bbsim{

BucketBasedSim::BucketBasedSim(const stim::Circuit& circ, uint k_max_):
                decoding_graph(qrc::to_decoding_graph(circ)), k_max(k_max_){
                n_detectors = circ.count_detectors();
                n_observables = circ.count_observables();
                create_the_sorted_event_vector_calc_prob_k();
                
}
BucketBasedSim::~BucketBasedSim(){
    for(ErrorEvent* ptr : sorted_events_vector)
        delete ptr;
}


void 
BucketBasedSim::combinationsHelper(const std::vector<fp_t>& nums, int k, int start, std::vector<int>& temp, std::vector<std::vector<int>>& result) {
    // base case: we've chosen k elements
    if (k == 0) {
        result.push_back(temp);
        return;
    }
    
    // recursive case: choose elements from start to end
    for (int i = start; i <= nums.size() - k; ++i) {
        temp.push_back(i);
        combinationsHelper(nums, k-1, i+1, temp, result);
        temp.pop_back();
    }
}

std::vector<std::vector<int>> 
BucketBasedSim::combinations(const std::vector<fp_t>& nums, int k) {
    std::vector<std::vector<int>> result;
    std::vector<int> temp;
    combinationsHelper(nums, k, 0, temp, result);
    return result;
}

void 
BucketBasedSim::create_the_sorted_event_vector_calc_prob_k(){
    probab_sum = 0;
    ErrorEvent* ee_ptr;
    qrc::DecodingGraph::Edge * edge_ptr;
    std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    uint n_vertices = list_of_vertices.size();
    for(uint i = 0 ; i < n_vertices; i++){
        for(uint j = i+1; j < n_vertices; j++){
            edge_ptr  = decoding_graph.get_edge(list_of_vertices[i], list_of_vertices[j]);
    // uint64_t s = 0;
    //std::cout << "n_detectors" << n_detectors << std::endl;
    // for (uint d = 0; d < n_detectors; d++) {
    //     auto v = decoding_graph.get_vertex(d);
    //     for (auto w : decoding_graph.adjacency_list(v)) {
    //         edge_ptr = decoding_graph.get_edge(v, w);
            if(edge_ptr == nullptr){
                continue;
                //std::cout <<"1"<< std::endl;
            }
            // s++;
            // std::cout << s << std::endl;
            probab_sum += edge_ptr->error_probability;
            probability_group[edge_ptr->error_probability]++;
            ee_ptr = new ErrorEvent(edge_ptr->id, edge_ptr->detectors.first, edge_ptr->detectors.second,
                                edge_ptr->edge_weight, edge_ptr->error_probability, edge_ptr->frames);
            //std::cout <<"1"<<std::endl;
            if(ee_ptr == nullptr){
                std::cout << "Error happened in creating an edge event";
                sorted_events_vector.clear();
            }

            sorted_events_vector.push_back(ee_ptr); 
            // std::cout <<"2"<<std::endl;
            list_of_error_probabs.push_back(ee_ptr->error_probability);  
            // std::cout <<"3"<<std::endl; 
        }
        // std::cout <<"6"<< std::endl;

    }
    std::sort(sorted_events_vector.begin(), sorted_events_vector.end(), [](const ErrorEvent* e1, const ErrorEvent* e2) {
    return e1->error_probability < e2->error_probability;});
    fp_t cumulative_prob = 0;
    ///
    // uint eventno = 0;
    // std::array<fp_t, 100'000> log_poly;
    // log_poly.fill(1);
    //
    prob_k.resize(k_max+1, 0);
    prob_k[0] = 1;
    for(ErrorEvent* e : sorted_events_vector){
        cumulative_prob += e->error_probability;
        e->set_cumulative_probability(cumulative_prob);
        
        fp_t pj = e->error_probability;
        for(int kk=k_max; kk>=1; kk--) {
            prob_k[kk] = prob_k[kk-1] * pj + prob_k[kk] * (1.0 - pj);
        }
        prob_k[0] *= (1.0 - pj);

        //For calculating the probability of bucket K and 
        // fp_t ep = e->error_probability;
        //     if (eventno == 0) {
        //         log_poly[0] = log(1-ep);
        //         log_poly[1] = log(ep);
        //     } else {
        //         std::array<fp_t, 100'000> log_pa, log_pb;
        //         log_pa.fill(1);
        //         log_pb.fill(1);
        //         for (uint i = 0; i <= eventno; i++) {
        //             log_pa[i] = log_poly[i] + log(1-ep);
        //             log_pb[i+1] = log_poly[i] + log(ep);
        //         }
        //         for (uint i = 0; i < log_poly.size(); i++) {
        //             if (log_pa[i] == 1 && log_pb[i] == 1) {
        //                 log_poly[i] = 1;
        //             } else if (log_pa[i] == 1) {
        //                 log_poly[i] = log_pb[i];
        //             } else if (log_pb[i] == 1) {
        //                 log_poly[i] = log_pa[i];
        //             } else {
        //                 log_poly[i] = log(pow(M_E, log_pa[i]) + pow(M_E, log_pb[i]));
        //             }
        //         }
        //     }
        //     eventno++;
    }
    
    number_of_events = sorted_events_vector.size();

    // for(int i = 0; i<k_max+1; i++){
    //     prob_k.push_back(pow(M_E, log_poly[i]));

    // }
    // for(uint k = 1; k < k_max+1 ; k++){
    //     std::cout << "k " << k << std::endl;
    //     std::vector<std::vector<int>> possible_errors = combinations(list_of_error_probabs, k);
    //     fp_t temp = 1;
    //     std::cout << "possible_errors size " << possible_errors.size() << std::endl;
    //     for(auto one_of_n_choose_k : possible_errors){
    //         fp_t temp = 1;
    //         for(int l = 0; l < number_of_events; l++){
    //             if(std::find(one_of_n_choose_k.begin(), one_of_n_choose_k.end(), l) != one_of_n_choose_k.end()){
    //                 temp *= list_of_error_probabs[l];
    //             }
    //             else{
    //                 temp *= (1-list_of_error_probabs[l]);
    //             }   
    //         }
    //         buckets_probabs[k] += temp;
    //     }
    // }
}

fp_t BucketBasedSim::return_probab_sum(){return probab_sum;}

void 
BucketBasedSim::print_decoding_graph(){
    std::cout << "__________print_graph___________" << std::endl;
    qrc::DecodingGraph::Edge * edge_ptr;
    qrc::DecodingGraph::Vertex* v1, *v2;
    std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    uint n_vertices = list_of_vertices.size();
    for(uint i = 0 ; i < n_vertices; i++){
        v1 = list_of_vertices[i];
        for(uint j = i+1; j < n_vertices; j++){
            v2 = list_of_vertices[j];
            edge_ptr  = decoding_graph.get_edge(list_of_vertices[i], list_of_vertices[j]);
            if(edge_ptr == nullptr)
                continue;
            std::cout << "(" << v1->detector << ", " << v2->detector << ") - E : [ id=" <<
                        edge_ptr->id << ", w = " << edge_ptr->edge_weight << " , er = " << edge_ptr->error_probability << " ]" << std::endl;    
        }
    }
    std::cout << std::endl << "__________end of print_graph___________" << std::endl;
}

void 
BucketBasedSim::print_vertices(){
    decoding_graph.set_vertices_features();
    //std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    for(uint i = 0; i < n_detectors; i++){
        qrc::DecodingGraph::Vertex* v = decoding_graph.get_vertex(i);
        std::cout <<" N"<<i <<" - " << v->feature << std::endl;
        std::vector<qrc::DecodingGraph::Vertex*> adj_list =  decoding_graph.adjacency_list(v);
        std::cout << "  Neighbors(" << adj_list.size() << "):";
        for(qrc::DecodingGraph::Vertex* v_ : adj_list){
            std::cout << v_->detector << " - ";
            qrc::DecodingGraph::Edge* e = decoding_graph.get_edge(v, v_);
            if(e == nullptr){
                std::cout << "This message should not be printed, existed edge with nulptr!";
            }
            std:: cout << e->error_probability << ", ";
        }
        std::cout << std::endl;
    }
    qrc::DecodingGraph::Vertex* v = decoding_graph.get_vertex(BOUNDARY_INDEX);
        
        std::cout <<" N"<<n_detectors <<" - " << v->feature << "(degree = " << v->n_neighbors << ")" <<std::endl ;
        std::vector<qrc::DecodingGraph::Vertex*> adj_list =  decoding_graph.adjacency_list(v);
        std::cout << "  Neighbors(" << adj_list.size() << "):";
        for(qrc::DecodingGraph::Vertex* v_ : adj_list){
            std::cout << v_->detector << " - ";
            qrc::DecodingGraph::Edge* e = decoding_graph.get_edge(v, v_);
            if(e == nullptr){
                std::cout << "This message should not be printed, existed edge with nulptr!";
            }
            std:: cout << e->error_probability << ", ";
        }
        std::cout << std::endl;
        std::cout << "________________________________________"<< std::endl;
        std::cout <<"Feature : Number of nodes with that Feature" << std::endl;
        for(const auto& f : decoding_graph.feature_group){
            std::cout << f.first << " : " << f.second << std::endl;
        }
        std::cout << "________________________________________"<< std::endl;
        std::cout <<"Degree of nodes : Number of nodes with that degree" << std::endl;
        for(const auto& f : decoding_graph.node_degree_group){
            std::cout << f.first << " : " << f.second << std::endl;
        }
        std::cout << "________________________________________"<< std::endl;
        std::cout <<"(feature, degree of node): Number of nodes with that feature and degree" << std::endl;
        for(const auto& f : decoding_graph.feature_degree_group){
            std::cout << "(" << f.first.first <<", " << f.first.second << "): "
             << f.second;
            // if(f.second == 1){
            //     std::cout << " v " << 
            // }
            std::cout << std::endl;
        }
}

void 
BucketBasedSim::print_weight_matrix(){
    uint n_events = sorted_events_vector.size();
    std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    uint n_vertices = list_of_vertices.size();
    std::vector<std::vector<fp_t>> m(n_vertices, std::vector<fp_t>(n_vertices, 0));

    auto path_table = qrc::compute_path_table(decoding_graph);
    uint i,j;
    ErrorEvent* ee_ptr;
    for(uint e = 0; e < n_events; e++){
        ee_ptr = sorted_events_vector[e];
        i = ee_ptr->detectors.first ;
        j= ee_ptr->detectors.second;
        if(i == BOUNDARY_INDEX){
            i = n_vertices-1;
        }
        if(j == BOUNDARY_INDEX){
            j = n_vertices-1;
        }
        // std::cout << i << " "<<j << std::endl;
        m[i][j] = ee_ptr->error_probability;
        m[j][i] = ee_ptr->error_probability;
    }


    for(uint i = 0; i<n_vertices; i++){
        for(uint j = 0; j< n_vertices; j++){
            std::cout<< m[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
}

void 
BucketBasedSim::print_path_matrix(){
    std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    uint n_vertices = list_of_vertices.size();
    std::vector<std::vector<fp_t>> m(n_vertices, std::vector<fp_t>(n_vertices, 0));

    auto path_table = qrc::compute_path_table(decoding_graph);
    uint i,j;
    for (auto v :decoding_graph.vertices()){
        for(auto w: decoding_graph.vertices()){
            auto v_w = std::make_pair(v,w);
            auto p = path_table[v_w].distance;
            i = v->detector;
            j=w->detector;
            if(v->detector == BOUNDARY_INDEX){
                i = n_vertices-1;
            }
            if(w->detector == BOUNDARY_INDEX){
                j = n_vertices-1;
            }
            // std::cout << i << " "<<j << std::endl;
            fp_t x = pow(10,-1*p);
            m[i][j] = x;;;
        }
    }


    for(uint i = 0; i<n_vertices; i++){
        for(uint j = 0; j< n_vertices; j++){
            std::cout<< m[i][j] << " ";
        }
        std::cout << std::endl;
    }
    
}

void 
BucketBasedSim::print_path_list(){
    std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    uint n_vertices = list_of_vertices.size();
    std::vector<std::vector<fp_t>> m(n_vertices, std::vector<fp_t>(n_vertices, 0));
    std::vector<std::vector<std::vector<qrc::DecodingGraph::Vertex*>>> paths(n_vertices, std::vector<std::vector<qrc::DecodingGraph::Vertex*>>(n_vertices));  
    auto path_table = qrc::compute_path_table(decoding_graph);
    uint i,j;
    for (auto v :decoding_graph.vertices()){
        for(auto w: decoding_graph.vertices()){
            auto v_w = std::make_pair(v,w);
            auto p = path_table[v_w].distance;
            i = v->detector;
            j=w->detector;
            if(v->detector == BOUNDARY_INDEX){
                i = n_vertices-1;
            }
            if(w->detector == BOUNDARY_INDEX){
                j = n_vertices-1;
            }
            // std::cout << i << " "<<j << std::endl;
            fp_t x = pow(10,-1*p);
            m[i][j] = x;
            paths[i][j] = path_table[v_w].path;
            
        }
    }


    for(uint i = 0; i<n_vertices; i++){
        for(uint j = i; j< n_vertices; j++){
            std::cout << "(" << i << ", " << j << "): "<< m[i][j] << " - ";
            
            for(const auto& p :paths[i][j]){
                if(p != nullptr){
                    std::cout << p->detector << "-";
                }
            }
            std::cout <<std::endl;
        }
        std::cout << std::endl;
    }
    
}

void 
BucketBasedSim::print_sorted_event_vector(){
    uint n_events = sorted_events_vector.size();
    ErrorEvent* ee_ptr;
    for(uint i = 0; i < n_events; i++){
        ee_ptr = sorted_events_vector[i];
        std::cout <<" E"<<i <<" - " << ee_ptr->id <<" : [ (" << ee_ptr->detectors.first << ", " << ee_ptr->detectors.second << ") " << " , er = " << ee_ptr->error_probability << " ]" << std::endl;    //" gp = "<< ee_ptr->cumulative_probability<< //"id= " <<ee_ptr->id  //", w = " << ee_ptr->edge_weight <<
    }
    std::cout << "______________________________________"<< std::endl;
    std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    uint n_vertices = list_of_vertices.size();
    std::cout << "Total # of Edges (Error Events)= " << n_events << " - Total # of Vertices (Syndrome Size) = " << n_vertices << std::endl;
    std::cout << "p : probability of edge - n : number of edges with probability p" << std::endl;
    for(auto p : probability_group){
        std::cout << "p: " << p.first << " - n: " << p.second << std::endl;
    }
    // std::cout << "______________________________________"<< std::endl;
    // for( uint i =  0; i< list_of_error_probabs.size(); i++){
    //     std::cout << "i: " << i << " p: " << list_of_error_probabs[i]<< std::endl;
    // }
    // std::cout << "______________________________________"<< std::endl;
    // for( uint i =  0; i< buckets_probabs.size(); i++){
    //     std::cout << "k: " << i+1 << " p: " << buckets_probabs[i]<< std::endl;
    // }
}

std::vector<std::vector<uint8_t>> 
BucketBasedSim::create_syndromes(uint k, uint64_t shots_per_k_, std::mt19937_64& rng, bool print_log){
    std::vector<std::vector<uint8_t>> syndromes_vect;
    std::vector<uint8_t> syndrome;
    uint n_detectors = decoding_graph.vertices().size()-1;
    syndrome.resize(n_detectors+1,((uint8_t)0)); //
    std::vector<bool> selected_event;
    selected_event.resize(sorted_events_vector.size(), false);

    std::uniform_real_distribution<> dis(0.0, probab_sum);
    for(uint i = 0; i< shots_per_k_; i++){
        std::fill(selected_event.begin(), selected_event.end(), false);
        for(uint j = 0; j<k; j++){
            fp_t r;
            uint m;

            // This while is for making sure that we do not select an event multiple times
            while(true){
                r = dis(rng);
                m = 0;
                while (sorted_events_vector[m]->cumulative_probability < r) {
                    m++;
                }
                if(selected_event[m]){
        
                    continue;
                }
                selected_event[m] = true;
                break;
            }
            uint d1 = sorted_events_vector[m]->detectors.first;
            uint d2 = sorted_events_vector[m]->detectors.second;

            if(d1 != BOUNDARY_INDEX){
                syndrome[d1] ^= ((uint8_t)1);//;
            }
            if(d2 != BOUNDARY_INDEX){
                syndrome[d2] ^= ((uint8_t)1);//((uint8_t)1);
            }
            for (uint obs : sorted_events_vector[m]->frames) {
                if (obs >= 0) {
                    syndrome[n_detectors+obs] ^= ((uint8_t)1);//((uint8_t)1);
                }
            }
            if(print_log){
                std::cout << "k:" << j << " d1: " << d1 << " d2: " << d2 << std::endl;
                qpd::print_syndrome(syndrome);
                std::cout << " --" << std::endl;
            }    
        }

        syndromes_vect.push_back(syndrome);

        std::fill(syndrome.begin(), syndrome.end(), ((uint8_t)0)); //((uint8_t)0)
    }

    return syndromes_vect;

}

stim::simd_bit_table
BucketBasedSim::create_syndromes_simd(uint k, uint64_t shots_per_k_, std::mt19937_64& rng, bool print_log){

    uint64_t syndrome_size =  n_detectors+ n_observables;
    stim::simd_bit_table syndromes_table(shots_per_k_, syndrome_size);
    std::vector<bool> selected_event;
    selected_event.resize(sorted_events_vector.size(), false);

    std::uniform_real_distribution<> dis(0.0, probab_sum);
    for(uint i = 0; i< shots_per_k_; i++){
        std::fill(selected_event.begin(), selected_event.end(), false);
        for(uint j = 0; j<k; j++){
            fp_t r;
            uint m;

            // This while is for making sure that we do not select an event multiple times
            while(true){
                r = dis(rng);
                m = 0;
                while (sorted_events_vector[m]->cumulative_probability < r) {
                    m++;
                }
                if(selected_event[m]){
        
                    continue;
                }
                selected_event[m] = true;
                break;
            }
            uint d1 = sorted_events_vector[m]->detectors.first;
            uint d2 = sorted_events_vector[m]->detectors.second;

            if(d1 != BOUNDARY_INDEX){
                syndromes_table[i][d1] ^= 1;//;
            }
            if(d2 != BOUNDARY_INDEX){
                syndromes_table[i][d2] ^= 1;//((uint8_t)1);
            }
            for (uint obs : sorted_events_vector[m]->frames) {
                if (obs >= 0) {
                    syndromes_table[i][n_detectors+obs] ^= 1;//((uint8_t)1);
                }
            }
            
            if(print_log){//(((j+1)*2)) < hw
                std::cout << "********************************&&&&&&&&&&&&&&&&&&"<< std::endl;
                std::vector<uint8_t> sss= qrc::_to_vector(syndromes_table[i], n_detectors, n_observables);
                uint hw = std::accumulate(sss.begin(), sss.begin()+n_detectors, 0);

                std::cout << "k:" << k <<" #err: "<<j << " d1: " << d1 << " d2: " << d2 << " hw:" << hw<< std::endl;
                //qpd::print_syndrome(sss);
                std::cout << " --" << std::endl;
                std::cout << "********************************&&&&&&&&&&&&&&&&&&"<< std::endl;
            }

        }
    }
    return syndromes_table;
}

void 
BucketBasedSim::set_shots_per_k(fp_t expected_ler, uint64_t max_n_iteration, bool set_shots_to_all, uint64_t filling_n_iterations, uint64_t min_n_iteration){
    double result, exponent, base;
    uint64_t shots;
    shots_per_k.clear();
    shots_per_k.resize(prob_k.size());
    for(int k =0; k<prob_k.size(); k++){
        result = prob_k[k] / expected_ler;
        exponent = std::floor(std::log10(result));
        base = std::round(result / std::pow(10, exponent));
        shots = 4*((uint64_t)base*std::pow(10, exponent));
        // std::cout << k << " " << shots << std::endl;
        if(max_n_iteration < shots){
            // if(shots < 10*max_n_iteration){
            //     shots = max_n_iteration;
            // }
            // else if(set_shots_to_all){
            //     shots = filling_n_iterations;
            // }
            // else{
            //     shots = 0;
            // }
            // std::cout << k << " " << shots << "max_n_iteration < shots" << std::endl;
            if(set_shots_to_all){
                shots = filling_n_iterations;
            }
            else{
                shots = 0;
            }
        }
        else if(shots !=0 && shots < min_n_iteration){
            // std::cout << k << " " << shots << "shots !=0 && shots < min_n_iteration" << std::endl;
            shots = min_n_iteration;
        }
        shots_per_k[k] = shots;


        // shots_per_k[k] = 0;
        // if(k == 13){
        //     shots_per_k[k] = 10'000'000'000;
        // }
        // else if(k == 14){
        //     shots_per_k[k] = 1'000'000'000;
        // }
        // else if(k == 15){
        //     shots_per_k[k] = 100'000'000;
        // }
        // else if(k == 16){
        //     shots_per_k[k] = 100'000'000;
        // }
    }

}


} //namespace bbsim