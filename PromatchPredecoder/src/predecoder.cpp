
#include "predecoder.h"
#include "decoding_graph.h"

/*
 *  author: Narges Alavisamani
 *  date:   6 Feburary 2023
 * */
namespace qpd {

Predecoder::Predecoder(const stim::Circuit& circ, bool general_)
            :qrc::Decoder(circ),
            general(general_),
            path_table(),
            prematching_subgraph(),
            decoding_graph(graph),
            n_detectors(circ.count_detectors()),
            n_observables(circ.count_observables())
            {
                path_table = qrc::compute_path_table(graph);
                prematching_subgraph.distance_with_neighbor = 1;
                decoding_graph.set_vertices_features();
                set_simulator_adjacenct_matrix();
                
            }

// void Predecoder::send_sydrome_to_predecoder(std::vector<uint8_t> syndrome_){
//     syndrome = syndrome_;
// }
const PredecoderSubgraph
    Predecoder::get_adjacency_matrix(){
        return prematching_subgraph;
    }

void Predecoder::print_graph(){
    std::cout << "__________print_graph___________" << std::endl;
    std::vector<qrc::DecodingGraph::Vertex*> _vertices = decoding_graph.vertices();
    for(qrc::DecodingGraph::Vertex* v : _vertices){
        std::vector<qrc::DecodingGraph::Vertex*> ad_list = decoding_graph.adjacency_list(v);
        std::cout << "\n___Vertex__ " << v->detector << std::endl;
        for (qrc::DecodingGraph::Vertex* adjecnt_v: ad_list){
            std::cout << adjecnt_v->detector << ", ";
        }
    }

    std::cout << std::endl << "__________end of print_graph___________" << std::endl;
}
std::string
Predecoder::name() {
    return "PREDECODER";
}

bool
Predecoder::is_software() {
    return true;
}

qrc::DecoderShotResult Predecoder::decode_error(const std::vector<uint8_t>& syndrome){
    qrc::DecoderShotResult pre_decode_res= {};

    return pre_decode_res;
}
void 
Predecoder::set_decoder(qrc::Decoder* decoder_){
    decoder = decoder_;
}

uint Predecoder::count_degree_one(){
    uint n_degree_1 = 0;
    auto subgraph = get_adjacency_matrix();
    uint n_vertices = subgraph.detector_list.size();
    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = subgraph.detector_list[vj];
        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);
        if (subgraph.adjacency_matrix[vdj].size() == 1){
            std::vector<qrc::DecodingGraph::Vertex*> neighbors = subgraph.adjacency_matrix[vdj];
            // Checking if the node with degree one is connected to a node with degree one.
            if(subgraph.adjacency_matrix[neighbors[0]].size() == 1)
                n_degree_1++;
        }
    }
    return n_degree_1;
}

PredecoderShotResult 
Predecoder::matching_ec_length_one(const std::vector<uint8_t>& syndrome_, bool add_boundry ){
    PredecoderShotResult result = {};

    // contains the detectors that are in a chain of length one
    std::vector<uint8_t> predecoder_syndrome;
    // contains detectors that are not matched during the predecoder
    std::vector<uint8_t> post_podecoder_syndrome;
    //std::cout << "XXXXXXX" <<syndrome_.size()<<std::endl;
    result.post_syndrome = syndrome_;
    //std::cout << ")----)))" << std::endl;
    result.prematch_syndrome.resize(syndrome_.size());
    //std::cout << ")))))))" << std::endl;
    std::fill(result.prematch_syndrome.begin(), result.prematch_syndrome.end(), 0);
    //std::cout << "*****" << std::endl;
    update_adjacency_matrix_and_detector_list(syndrome_, 1, add_boundry);
    //std::cout << "^^^^^" << std::endl;
    auto subgraph = get_adjacency_matrix();
    uint n_vertices = subgraph.detector_list.size();

    if(subgraph.distance_with_neighbor != 1){
        std::cout << "Cannot use the corrent adjacent matrix becuase it is not built with distance 1" << std::endl;
        return result;
    }

    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = subgraph.detector_list[vj];
        // std::cout << "a" << std::endl;
        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);
        //std::cout << "b" << std::endl;
        if (subgraph.adjacency_matrix[vdj].size() == 1){
            //this neighbors should have only one element as "adjacency_matrix[vdj].size() == 1"
            std::vector<qrc::DecodingGraph::Vertex*> neighbors = subgraph.adjacency_matrix[vdj];
            //std::cout << "c" << std::endl;
            // Checking if the node with degree one is connected to a node with degree one.
            if(subgraph.adjacency_matrix[neighbors[0]] .size() == 1){
                //match two detectors of error chain length one
                result.pre_matches[vdj->detector] = neighbors[0]->detector;
                //std::cout << "d" << std::endl;

                // make the prematch syndrome
                if(vdj->detector!=BOUNDARY_INDEX)
                    result.prematch_syndrome[vdj->detector]=((uint8_t)1);
                //std::cout << "e" << std::endl;
                // for(auto x: neighbors){
                //     std::cout << x->detector << ",";
                // }
                //std::cout<<"end:" <<neighbors[0]->detector<<std::endl;
                if(neighbors[0]->detector!=BOUNDARY_INDEX)
                    result.prematch_syndrome[neighbors[0]->detector]=((uint8_t)1);
                // std::cout << "f" << std::endl;
                // Deleting the detectors that are prematched in the predecoder
                if(vdj->detector!=BOUNDARY_INDEX)
                    result.post_syndrome[vdj->detector]=((uint8_t)0);
                if(neighbors[0]->detector!=BOUNDARY_INDEX)
                    result.post_syndrome[neighbors[0]->detector]=((uint8_t)0);
            }
        }
    }

    return result;

}

void Predecoder::print_subgraph(const std::vector<uint8_t>& syndrome_){    

    update_adjacency_matrix_and_detector_list(syndrome_, 1, true);
    auto subgraph = get_adjacency_matrix();
    std::cout << "__________Subgraph__________" << std::endl;
    for (auto const& entry : subgraph.adjacency_matrix) {
        std::cout << "___Vertex___ " << entry.first->detector << std::endl;
        std::cout << "__Neighbours__ ";
        for (auto const& neighbor : entry.second) {
            std::cout << neighbor->detector << ", ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << "__________Subgraph Info__________" << std::endl;


    auto components = qpd::find_connected_components(subgraph.adjacency_matrix);
    std::cout << "# of connected components is " << components.size() << " :" <<std::endl;
    for(auto const& comp : components){
        std::cout<<comp.size() << ": ";
        for(auto const& v: comp){
            std::cout << v->detector << " - ";
        }
        std::cout << std::endl;
    }
    std::cout << " *Neighbor distance: " << subgraph.distance_with_neighbor <<std::endl;


    std::cout <<std::endl << "__________End of Subgraph__________" << std::endl;
    
}

std::vector<uint8_t>
Predecoder::get_correction_from_matching(const std::map<uint, uint>& matching) {
    // std::set<uint> visited;
    // std::vector<uint8_t> correction(circuit.count_observables(), 0);
    // for (auto di_dj : matching) {
    //     uint di = di_dj.first;
    //     uint dj = di_dj.second;
    //     if (visited.count(di) || visited.count(dj)) {
    //         continue;
    //     }
    //     qrc::DecodingGraph::Vertex * vdi = decoding_graph.get_vertex(di);
    //     qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);
    //     auto edge = graph.get_edge(vdi, vdj);
    //     // The edge should exist.
    //     for (uint obs : edge->frames) {
    //         // Flip the bit.
    //         if (obs >= 0) {
    //             correction[obs] = !correction[obs];
    //         }
    //     }
    //     visited.insert(di);
    //     visited.insert(dj);
    // }
    // return correction;
    std::set<uint> visited;
    std::vector<uint8_t> correction(circuit.count_observables(), 0);
    for (auto di_dj : matching) {
        uint di = di_dj.first;
        uint dj = di_dj.second;
        if (visited.count(di) || visited.count(dj)) {
            continue;
        }
        // Check path between the two detectors.
        // This is examining the error chain.
        qrc::DecodingGraph::Vertex * vdi = graph.get_vertex(di);
        qrc::DecodingGraph::Vertex * vdj = graph.get_vertex(dj);
        auto vdi_vdj = std::make_pair(vdi, vdj);
        std::vector<qrc::DecodingGraph::Vertex*> detector_path(path_table[vdi_vdj].path);
        for (uint i = 1; i < detector_path.size(); i++) {
            // Get edge from decoding graph.
            auto wi = detector_path[i-1];
            auto wj = detector_path[i];
            auto edge = graph.get_edge(wi, wj);
            // The edge should exist.
            for (uint obs : edge->frames) {
                // Flip the bit.
                if (obs >= 0) {
                    correction[obs] = !correction[obs];
                }
            }
        }
        
        visited.insert(di);
        visited.insert(dj);
    }
    return correction;
}

PostPredecoderResult 
Predecoder::predecode_error(PredecoderShotResult prematching_result){
    PostPredecoderResult result = {};
    result.post_syndrome = prematching_result.post_syndrome;
    result.correction = get_correction_from_matching(prematching_result.pre_matches);

    return result;
}

void Predecoder::print_shot_info(const std::vector<uint8_t>& syndrome,
  qpd::PredecoderShotResult& predecoder_results, qrc::DecoderShotResult &predecoder_decoder_results, 
  qrc::DecoderShotResult &decoder_results, int n_1s, bool only_predecoder_matches, bool full_info, bool print_syndromes){
    
    std::array<fp_t,20> weight_arrays_pre;
    weight_arrays_pre.fill(1);
    std::array<fp_t,20> weight_arrays;
    weight_arrays.fill(1);
    uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + n_detectors, 0);
    std::cout<< "_____________________INFO_____________________" << std::endl;
    //std::cout<< "Distance: "<< distance << " - Physical ER: " << physcial_error << " - Meas. ER:" << meas_er << std::endl;
    std::cout<< "HW: " << hw << ",  # (|ec| = 1) : " << n_1s << std::endl;
    std::cout<< "______________________________________________" << std::endl;
    if(print_syndromes){
        std::cout << "Syndrome:\n";
        for (int id =0; id<syndrome.size();id++) {
            std::cout<< id<<":"<< ((int)syndrome[id] )<< ", ";
        }
        std::cout << std::endl;
        std::cout << "Predecoder Syndrome:\n";
        for (int id =0; id<predecoder_results.prematch_syndrome.size();id++){
            std::cout<< id<<":" << ((int)predecoder_results.prematch_syndrome[id] ) << ", ";
        }
        std::cout << std::endl;
        std::cout << "Post Syndrome:\n";
        for (int id =0; id<predecoder_results.post_syndrome.size();id++){
            std::cout<< id<<":" << ((int)predecoder_results.post_syndrome[id] ) << ", ";
        }
    }
    //std::stringstream ss;
    std::cout<< std::endl;
    // std::cout<< std::endl;

    std::cout<< "_____________________Prematching MAP_____________________" << std::endl;
    for (const auto& kv : predecoder_results.pre_matches) {
        std::cout<< "V1: " << kv.first << ", V2: " << kv.second << std::endl;
    }
    std::cout<< "_____________________End of MAP_____________________" << std::endl;

    if(full_info){
        std::cout<< "_____________________Decoder Matching MAP After Predecoding_____________________" << std::endl;
        for (const auto& kv : predecoder_decoder_results.matching) {
            std::cout<< "V1: " << kv.first << ", V2: " << kv.second << std::endl;
        }
        std::cout<< "_____________________End of Pre+Dec MAP_____________________" << std::endl;
    }
    
    std::cout << "_________________Different Matches_________________" << std::endl;
    std::cout << "Detector |   Matched in predecoder: error prob.  | Matched in decoder : error prob. " << std::endl;
    std::array<uint, 2> matching_pair {0, 0};
    for(const auto& kv : predecoder_results.pre_matches){
        // matching_pair[0] = kv.first;
        // matching_pair[1]=kv.second;
        // for (uint idx : matching_pair){

        // This part is to compare the matches for predecoder vs decoders
        uint idx = kv.first;
            if (predecoder_results.pre_matches[idx] != decoder_results.matching[idx]){
                auto v = decoding_graph.get_vertex(idx);

                // this w is the detector that matched with idx in the predecoder
                auto w = decoding_graph.get_vertex(predecoder_results.pre_matches[idx]);

                auto e = decoding_graph.get_edge(v,w);
                if(e == nullptr){
                    std::cout << "Error No Edge in the prematech!!! V1: " <<
                    v->detector << " V2: " <<w->detector;
                    return;
                }
                fp_t predecoder_path_weight = e->error_probability;

                // this w is for the detector that matched with idx in the decoder
                w = decoding_graph.get_vertex(decoder_results.matching[idx]);
                auto v_w = std::make_pair(v, w);
                // Finding error chains from the path table of decoding graph 
                auto ec = path_table[v_w].path;
                const uint path_length = ec.size();
                auto previous_v = v;
                fp_t path_error_probability = 1;
                weight_arrays.fill(1);
                auto mid_edge = decoding_graph.get_edge(previous_v, ec[1]);
                if(mid_edge == nullptr){
                    previous_v = w;
                }

                for(uint k=1; k<path_length; k++ ){
                    auto mid_edge = decoding_graph.get_edge(previous_v, ec[k]);
                    if(mid_edge == nullptr){
                        std::cout << "Error No Edge in the decoder match!!!";
                    return;
                    }
                    weight_arrays[k] = mid_edge->error_probability;
                    previous_v = ec[k];
                    path_error_probability *= weight_arrays[k];
                }
                std::cout<<std::endl;
                std::cout << idx << "  |  " << predecoder_results.pre_matches[idx]  << " : " <<predecoder_path_weight << "  |  " << 
                decoder_results.matching[idx] << " : " << path_error_probability;
                
                if (path_length == 2){
                    std::cout << std::endl;
                    continue;
                }
                
                std::cout << "=" << weight_arrays[1];

                for(uint k=2; k<path_length; k++ ){
                    std::cout << "*"<< weight_arrays[k];
                }
                std::cout << std::endl;
            }
        //} // for matching_pair[0] = kv.first; thingy
    }
    if(only_predecoder_matches)
        return;
    
    std::cout<< "-------------------^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^-------------------" << std::endl;
    std::cout << "Detector |   Matched in decoder after predecoder: error prob.  | Matched in decoder : error prob. " << std::endl;


    //// For finding different matchings for predecoer+"decoder" vs "decoder"
    for(const auto& kv : predecoder_decoder_results.matching){
        // matching_pair[0] = kv.first;matching_pair[1]=kv.second;
        // for (uint idx : matching_pair){
        uint idx = kv.first;
            if (predecoder_decoder_results.matching[idx] != decoder_results.matching[idx]){
                auto v = decoding_graph.get_vertex(idx);

                // w that matched with idx in the decoder after prematching 
                auto w = decoding_graph.get_vertex(predecoder_decoder_results.matching[idx]);
                auto v_w = std::make_pair(v, w);
                // Finding error chains from the path table of decoding graph 
                auto ec = path_table[v_w].path;
                const uint path_length_pre = ec.size();
                auto previous_v = v;
                fp_t path_error_probability_pre = 1;
                weight_arrays_pre.fill(1);
                auto mid_edge = decoding_graph.get_edge(previous_v, ec[1]);
                
                if(mid_edge == nullptr){
                    previous_v = w;
                }

                for(uint k=1; k<path_length_pre; k++ ){
                    mid_edge = decoding_graph.get_edge(previous_v, ec[k]);
                    if(mid_edge == nullptr){
                        std::cout << "Error No Edge in the decoder match in pre+dec!!! V1: "<<
                        previous_v->detector << " V2: "<< ec[k]->detector << std::endl;
                        std::cout << "The chain is: ";
                        for (auto x : ec){
                            std::cout << x->detector << ", ";
                        }
                        // std::cout << std::endl << "The original matching pair is: V1: " << 
                        // v->detector << " V2: " << w->detector; 
                        // std::cout << " first_null_check_taken: "<< first_null_check_taken <<std::endl; 
                    return;
                    }
                    //first_null_check_taken = false;
                    weight_arrays_pre[k] = mid_edge->error_probability;
                    previous_v = ec[k];
                    path_error_probability_pre *= weight_arrays_pre[k];
                }
                //////////for decoder only
                // w that matched with idx in the decoder, when we only use decoder
                w = decoding_graph.get_vertex(decoder_results.matching[idx]);
                v_w = std::make_pair(v, w);
                // Finding error chains from the path table of decoding graph 
                ec = path_table[v_w].path;
                const uint path_length = ec.size();
                previous_v = v;
                fp_t path_error_probability = 1;
                weight_arrays.fill(1);
                mid_edge = decoding_graph.get_edge(previous_v, ec[1]);
                //first_null_check_taken = false;
                if(mid_edge == nullptr){
                    previous_v = w;
                }

                for(uint k=1; k<path_length; k++ ){
                    mid_edge = decoding_graph.get_edge(previous_v, ec[k]);
                    if(mid_edge == nullptr){
                        std::cout << "Error No Edge in the decoder match, in dec only!!!V1: "<<
                        previous_v->detector << " V2: "<< ec[k]->detector << std::endl;
                        for (auto x : ec){
                            std::cout << x->detector << ", ";
                        }
                    return;
                    }
                    weight_arrays[k] = mid_edge->error_probability;
                    previous_v = ec[k];
                    path_error_probability *= weight_arrays[k];
                }


                ///Printing everything out
                std::cout<<std::endl;
                std::cout << idx << "  |  " << predecoder_decoder_results.matching[idx]  << " : " <<path_error_probability_pre;
                
                if(2 < path_length_pre){
                    std::cout << "=" << weight_arrays_pre[1];

                    for(uint k=2; k<path_length_pre; k++ ){
                        std::cout << "*"<< weight_arrays_pre[k];
                    }
                }
                
                std::cout << "  |  " <<  decoder_results.matching[idx] << " : " << path_error_probability;
                
                if (path_length == 2){
                    std::cout << std::endl;
                    continue;
                }
                
                std::cout << "=" << weight_arrays[1];

                for(uint k=2; k<path_length; k++ ){
                    std::cout << "*"<< weight_arrays[k];
                }
                std::cout << std::endl;
            }
        //} //  for thatkv.first
    }
    std::cout <<  "_____________________@END OF THIS ROUND@______________________________" << std::endl;

}

PredecoderShotResult
Predecoder::matching_ec(uint ec_length, const std::vector<uint8_t>& syndrome_, bool add_boundry){
    PredecoderShotResult result = {};

    // contains the detectors that are in a chain of length one
    std::vector<uint8_t> predecoder_syndrome;
    // contains detectors that are not matched during the predecoder
    std::vector<uint8_t> post_podecoder_syndrome;;
    uint max_detector = n_detectors;
    result.post_syndrome = syndrome_;
    result.prematch_syndrome.resize(syndrome_.size());
    std::fill(result.prematch_syndrome.begin(), result.prematch_syndrome.end(), 0);
    update_adjacency_matrix_and_detector_list(syndrome_, ec_length, add_boundry);
    auto subgraph = get_adjacency_matrix();
    uint n_vertices = subgraph.detector_list.size();

    if(ec_length != subgraph.distance_with_neighbor){
        std::cout << "Requested ec_length (" << ec_length 
        << ") is different from the distance of neighbors (" << subgraph.distance_with_neighbor << ")!"<< std::endl;
    }

    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = subgraph.detector_list[vj];
        //std::cout << "a" << std::endl;
        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);
        //std::cout << "b" << std::endl;
        if (subgraph.adjacency_matrix[vdj].size() == 1){
            //this neighbors should have only one element as "adjacency_matrix[vdj].size() == 1"
            std::vector<qrc::DecodingGraph::Vertex*> neighbors = subgraph.adjacency_matrix[vdj];
            // std::cout << "Node degree one found! vdj: " << vdj->detector << 
            //             " # of neighbor's neighbors: " << adjacency_matrix[neighbors[0]] .size() << std::endl;
            //std::cout << "c" << std::endl;
            // Checking if the node with degree one is connected to a node with degree one.
            if(subgraph.adjacency_matrix[neighbors[0]] .size() == 1){
                // std::cout << "One matching found between " << vdj->detector << 
                // " and " << neighbors[0]->detector << std::endl;
                //match two detectors of error chain length one
                result.pre_matches[vdj->detector] = neighbors[0]->detector;
                //std::cout << "d" << std::endl;

                // make the prematch syndrome
                if(vdj->detector!=BOUNDARY_INDEX)
                    result.prematch_syndrome[vdj->detector]=((uint8_t)1);
                //std::cout << "e" << std::endl;
                // for(auto x: neighbors){
                //     std::cout << x->detector << ",";
                // }
                //std::cout<<"end:" <<neighbors[0]->detector<<std::endl;
                if(neighbors[0]->detector!=BOUNDARY_INDEX)
                    result.prematch_syndrome[neighbors[0]->detector]=((uint8_t)1);
                //std::cout << "f" << std::endl;
                // Deleting the detectors that are prematched in the predecoder
                if(vdj->detector!=BOUNDARY_INDEX)
                    result.post_syndrome[vdj->detector]=((uint8_t)0);
                if(neighbors[0]->detector!=BOUNDARY_INDEX)
                    result.post_syndrome[neighbors[0]->detector]=((uint8_t)0);
            }

        }
    }

    return result;
    
}

void Predecoder::update_adjacency_matrix_and_detector_list(const std::vector<uint8_t>& syndrome_
                            , uint ec_length, bool add_boundry){
        
    uint8_t syndrome_is_even = 0x1;
    prematching_subgraph.detector_list.clear();
    prematching_subgraph.distance_with_neighbor = ec_length;
    uint max_detector = n_detectors;
    uint n_observable = circuit.count_observables();
    for (uint di = 0; di < n_detectors; di++) {
        auto syndrome_bit = syndrome_[di];
        if (di > max_detector) {
            syndrome_bit = 0;
        }
        if (syndrome_bit) {
            syndrome_is_even ^= 0x1;
            prematching_subgraph.detector_list.push_back(di);
        }
    }

    if (!syndrome_is_even && add_boundry) {
        // Add boundary to matching graph.
        prematching_subgraph.detector_list.push_back(BOUNDARY_INDEX);
    }
    uint n_vertices = prematching_subgraph.detector_list.size();
    prematching_subgraph.adjacency_matrix.clear();
    for (uint vi = 0; vi < n_vertices; vi++) {
        uint di = prematching_subgraph.detector_list[vi];
        qrc::DecodingGraph::Vertex * vdi = decoding_graph.get_vertex(di);
        //resetting_number of depenedence. This value is calculated in get_fast_matching_group_all_neighbors
        vdi->n_dependents = 0;
        prematching_subgraph.edge_lists[vdi] = decoding_graph.get_vertex_edges_list(vdi);
        for (uint vj = vi + 1; vj < n_vertices; vj++) {
            uint dj = prematching_subgraph.detector_list[vj];
            qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);
            // for the if statement to check if the detector has a neighbor or not
            bool add_neighbor = false;
            // if ec_length is greater than 1, we check the path in path table if the length of the path
            // is =< 1 we add the neighbor
            if(1 < ec_length ){
                auto vdi_vdj = std::make_pair(vdi, vdj);
                auto ec = path_table[vdi_vdj].path;
                add_neighbor = (((ec.size() - 1) < ec_length ) || ((ec.size() - 1) == ec_length));
            }
            // for ec_length == 1, we check if there is a edge between two detectors 
            else{
                qrc::DecodingGraph::Edge* edge = decoding_graph.get_edge(vdi, vdj);
                add_neighbor = (edge != nullptr);
            }
            
            if (add_neighbor){
                // Add vdi to the vector associated with the key of vdj.
                prematching_subgraph.adjacency_matrix[vdj].push_back(vdi);
                // Add vdj to the vector associated with the key of vdi.
                prematching_subgraph.adjacency_matrix[vdi].push_back(vdj);
            }
            add_neighbor = false;
            
        }
    }
    

}
void 
Predecoder::is_matched_parity_isolated(SingleMatchingInfo& m_info){

    //std::cout << "^^^^^" << std::endl;
    auto subgraph = get_adjacency_matrix();
    uint n_vertices = subgraph.detector_list.size();

    if(subgraph.distance_with_neighbor != 1){
        std::cout << "Cannot use the corrent adjacent matrix becuase it is not built with distance 1" << std::endl;
    }

    
    qrc::DecodingGraph::Vertex * v1 = decoding_graph.get_vertex(m_info.first);
    qrc::DecodingGraph::Vertex * v2 = decoding_graph.get_vertex(m_info.second);
    if(m_info.length == 1) {
        //std::cout << "b" << std::endl;
        if (subgraph.adjacency_matrix[v1].size() == 1
            && subgraph.adjacency_matrix[v2].size() == 1){
            //std::cout << "c" << std::endl;
            // Checking if the node with degree one is connected to a node with degree one.
            if(subgraph.adjacency_matrix[v1][0] != v2){
                std::cout << "SOMETHING IS WRONG! Matched nodes with ec length one are not adjacent! " << std::endl;
            }
            m_info.isolated = true;
        }
        else{
            m_info.isolated = false;
        }
    }
    else{
        if (subgraph.adjacency_matrix[v1].size() == 0
            && subgraph.adjacency_matrix[v2].size() == 0){
                m_info.isolated = true;
            }
        else{
            m_info.isolated = false;
            m_info.non_isolated_inf = new NoneIsolatedInf;
            m_info.non_isolated_inf->v1_adjacent_flipped_parity = subgraph.adjacency_matrix[v1];
            m_info.non_isolated_inf->v2_adjacent_flipped_parity = subgraph.adjacency_matrix[v2];
            for(const auto& adj1 : m_info.non_isolated_inf->v1_adjacent_flipped_parity){
                fp_t edge_prob = decoding_graph.get_edge(v1,adj1)->error_probability;
                m_info.non_isolated_inf->v1_adjacent_prob.push_back(edge_prob);
            }
            for(const auto& adj2 : m_info.non_isolated_inf->v2_adjacent_flipped_parity){
                fp_t edge_prob = decoding_graph.get_edge(v2,adj2)->error_probability;
                m_info.non_isolated_inf->v2_adjacent_prob.push_back(edge_prob);
            }

        }

    }
            

}

std::vector<qrc::DecodingGraph::Vertex*> 
Predecoder::get_fast_matching_group(const std::vector<uint8_t>& syndrome_){
    /*
    This function finds the fast matching group based on the 
    comparing the feature of the nodes and the edge probability.
    
    */

    //update the prematching subgraph (the subgraph that only has flipped 
    // parity bits.)
    update_adjacency_matrix_and_detector_list(syndrome_, 1, true);
    sorted_predecoding_edges.clear();

    std::vector<qrc::DecodingGraph::Vertex*> fast_group;
    bool neighbor_added = false;
    qrc::DecodingGraph::Vertex* v1;
    qrc::DecodingGraph::Edge* e;
    for(const auto& element : prematching_subgraph.adjacency_matrix){
        v1 = element.first;
        neighbor_added = false;
        // Boundary vertex is always in group2
        
        for(const auto& v2:  element.second){
            if(prematching_subgraph.adjacency_matrix[v2].size() ==1){
                    v1->n_dependents++;
            }
            // Boundary vertex is always in group2
            if(v2->detector == BOUNDARY_INDEX || v1->detector == BOUNDARY_INDEX){
                continue;
            }
            e = decoding_graph.get_edge(v1,v2);
            if(e == nullptr){
                std::cout <<"THIS MESSAGE SHOULD NOT BE PRINTED. NULLPTR FOR EDGE OF ADJACENT NODES." << std::endl;
            }
            else{
                if((v1->feature <= (e->error_probability*10) ) || (v2->feature <= (e->error_probability*10))){
                    if(!isDuplicate(fast_group, v2->detector)){
                        fast_group.push_back(v2);
                        neighbor_added = true;
                        sorted_predecoding_edges.push_back(e);
                    }
                }
            }

        }
        if(neighbor_added){
            addUniqueVertex(fast_group,v1);
        }

    }
    for(auto e : sorted_predecoding_edges){
        auto v1 = decoding_graph.get_vertex(e->detectors.first);
        auto v2 = decoding_graph.get_vertex(e->detectors.second);
        uint v1_adding_d_zero = 0;
        uint v2_adding_d_zero = 0;
        if(v1->n_dependents == 1){
            v2_adding_d_zero = v2->n_dependents - 1;
        }
        else{
            v2_adding_d_zero = v2->n_dependents;
        }
        if(v2->n_dependents == 1){
            v1_adding_d_zero = v1->n_dependents - 1;
        }
        else{
            v1_adding_d_zero = v1->n_dependents;
        }
        e->probability_impact = e->error_probability*(pow(10,(-4*(v2_adding_d_zero+v1_adding_d_zero))));
    }
    std::sort(sorted_predecoding_edges.begin(), sorted_predecoding_edges.end(), compare_edges_by_future_error_probability);
    return fast_group;  

}

std::vector<qrc::DecodingGraph::Vertex*> 
Predecoder::get_fast_matching_group_all_neighbors(const std::vector<uint8_t>& syndrome_){
    std::vector<qrc::DecodingGraph::Vertex*> fast_group;
    /* 
    This function put all the nodes that has at least one adjacent flipped
    node to the fast group.
    */
    //update the prematching subgraph (the subgraph that only has flipped 
    // parity bits.)
    update_adjacency_matrix_and_detector_list(syndrome_, 1, true);
    sorted_predecoding_edges.clear();
    degree_zero_vertices.clear();

    //auto subgraph = get_adjacency_matrix();
    uint n_vertices = prematching_subgraph.detector_list.size();

    qrc::DecodingGraph::Edge* e;

    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = prematching_subgraph.detector_list[vj];
        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);

        // at the beginning every node just know that it is in the cluster.
        vdj->n_cluster_members = 1;


        //This if exclude the cases that degree is one because of boundary vertex
        if (prematching_subgraph.adjacency_matrix[vdj].size() != 0){
            
            addUniqueVertex(fast_group,vdj);
            std::vector<qrc::DecodingGraph::Vertex*> neighbors = prematching_subgraph.adjacency_matrix[vdj];
            
            // Checking if the node with degree one is connected to a node with degree one.
                //match two detectors of error chain length one
                
            for(const auto& n: neighbors){
                addUniqueVertex(fast_group, n);
                e = decoding_graph.get_edge(vdj,n);
                if(e == nullptr){
                    std::cout <<"THIS MESSAGE SHOULD NOT BE PRINTED. NULLPTR FOR EDGE OF ADJACENT NODES." << std::endl;
                }
                addUniqueEdge(sorted_predecoding_edges, e);
                if(prematching_subgraph.adjacency_matrix[n].size() ==1){
                    vdj->n_dependents++;
                }
                 
            }
        
        }
        else{
            degree_zero_vertices.push_back(vdj);
        }
    }
    
    number_of_edges = sorted_predecoding_edges.size();

    std::sort(sorted_predecoding_edges.begin(), sorted_predecoding_edges.end(), compare_edges_by_error_probability);


    return fast_group;
}

std::vector<qrc::DecodingGraph::Vertex*> 
Predecoder::get_fast_matching_group_v2(const std::vector<uint8_t>& syndrome_){
    std::vector<qrc::DecodingGraph::Vertex*> fast_group;
    /*
    This function create the fast matching group by adding all the nodes that
    has only one adjacent flipped bit (isolated pairs)
    */
    //update the prematching subgraph (the subgraph that only has flipped 
    // parity bits.)
    update_adjacency_matrix_and_detector_list(syndrome_, 1, true);
    auto subgraph = get_adjacency_matrix();
    uint n_vertices = subgraph.detector_list.size();

    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = subgraph.detector_list[vj];
        
        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);

        if (subgraph.adjacency_matrix[vdj].size() == 1){
            
            std::vector<qrc::DecodingGraph::Vertex*> neighbors = subgraph.adjacency_matrix[vdj];
            
            // Checking if the node with degree one is connected to a node with degree one.
            if(subgraph.adjacency_matrix[neighbors[0]] .size() == 1){
                //match two detectors of error chain length one
                addUniqueVertex(fast_group,vdj);
                addUniqueVertex(fast_group, neighbors[0]);
            }
        }
    }

    return fast_group;
}

std::vector<uint8_t>
Predecoder::set_the_priorities(const std::vector<uint8_t>& syndrome_, bool all_in){
    std::vector<uint8_t> md_syndrome(syndrome_.begin(), syndrome_.end());
    //uint hw = std::accumulate(md_syndrome.begin(), md_syndrome.)
    // subgraph is made inside get_fast_matching_group
    std::vector<qrc::DecodingGraph::Vertex*> fast_group = get_fast_matching_group(syndrome_);
    auto subgraph = get_adjacency_matrix();
    pf_pairs.clear();
    pf_pairs.resize(6);
    pf_groups.clear();
    pf_groups.resize(6);
    predecoding_vertices.clear();
    predecoding_edges.clear();
    number_of_edges = 0;
    bool isolated;
    qrc::DecodingGraph::Edge* e;
    qrc::DecodingGraph::Vertex* v1;
    qrc::DecodingGraph::Vertex* v2;

    for(uint i = 0; i<sorted_predecoding_edges.size(); i++){
        isolated = false;
        e = sorted_predecoding_edges[i];
        v1 = decoding_graph.get_vertex(e->detectors.first);
        v2 = decoding_graph.get_vertex(e->detectors.second);
        if(!isDuplicate(fast_group, v1->detector) || !isDuplicate(fast_group,v2->detector)){
            continue;
        }
        // In this point, we are sure that those vertices are in the fast group
        // and have not been used by other pairs.

        // We first pop the vertices from fast_group
        delete_vertex_by_detector(fast_group, v1->detector);
        delete_vertex_by_detector(fast_group, v2->detector);
        predecoding_vertices.push_back(v1);
        predecoding_vertices.push_back(v2);
        predecoding_edges.push_back(e);
        number_of_edges += 2;
        if(v1->detector == BOUNDARY_INDEX || v2->detector == BOUNDARY_INDEX ){
            std::cout << "SOMETHING IS WRONG!! THIS MESSAGE SHOULD BE PRINTED! A BOUNDARY APPEARD IN PRIORITY SETTING!" <<std::endl;
        }
        else{
            md_syndrome[v1->detector] = ((uint8_t)0);
            md_syndrome[v2->detector] = ((uint8_t)0);
        }
        pf_pairs[0].push_back(e);
        pf_groups[0].push_back(v1);
        pf_groups[0].push_back(v2);


    }
    return md_syndrome;


}

std::vector<uint8_t>
Predecoder::prioritize_and_set_potential_matching(const std::vector<uint8_t>& syndrome_){
    std::vector<uint8_t> md_syndrome(syndrome_.begin(), syndrome_.end());
    //uint hw = std::accumulate(md_syndrome.begin(), md_syndrome.)
    // subgraph is made inside get_fast_matching_group
    std::vector<qrc::DecodingGraph::Vertex*> fast_group = get_fast_matching_group(syndrome_);
    auto subgraph = get_adjacency_matrix();
    pf_pairs.clear();
    pf_pairs.resize(6);
    pf_groups.clear();
    pf_groups.resize(6);
    predecoding_vertices.clear();
    predecoding_edges.clear();
    number_of_edges = 0;
    bool isolated;
    qrc::DecodingGraph::Edge* e;
    qrc::DecodingGraph::Vertex* v1;
    qrc::DecodingGraph::Vertex* v2;

    for(uint i = 0; i<sorted_predecoding_edges.size(); i++){
        isolated = false;
        e = sorted_predecoding_edges[i];
        v1 = decoding_graph.get_vertex(e->detectors.first);
        v2 = decoding_graph.get_vertex(e->detectors.second);
        if(!isDuplicate(fast_group, v1->detector) || !isDuplicate(fast_group,v2->detector)){
            continue;
        }
        // In this point, we are sure that those vertices are in the fast group
        // and have not been used by other pairs.

        //Only this added for a check on dependents
        bool no_poor = (((v1->n_dependents+v2->n_dependents) == 0) ||
                        (v1->n_dependents == 1&&prematching_subgraph.adjacency_matrix[v2].size() == 1) ||
                        (v2->n_dependents == 1&&prematching_subgraph.adjacency_matrix[v1].size() == 1) );
        if(!no_poor){
            continue;
        }

        // We first pop the vertices from fast_group
        delete_vertex_by_detector(fast_group, v1->detector);
        delete_vertex_by_detector(fast_group, v2->detector);
        predecoding_vertices.push_back(v1);
        predecoding_vertices.push_back(v2);
        predecoding_edges.push_back(e);
        number_of_edges += 2;
        if(v1->detector == BOUNDARY_INDEX || v2->detector == BOUNDARY_INDEX ){
            std::cout << "SOMETHING IS WRONG!! THIS MESSAGE SHOULD BE PRINTED! A BOUNDARY APPEARD IN PRIORITY SETTING!" <<std::endl;
        }
        else{
            md_syndrome[v1->detector] = ((uint8_t)0);
            md_syndrome[v2->detector] = ((uint8_t)0);
        }
        pf_pairs[0].push_back(e);
        pf_groups[0].push_back(v1);
        pf_groups[0].push_back(v2);


    }
    return md_syndrome;

}

std::vector<EnsembleEntry>
Predecoder::create_syndromes_ensembles(const std::vector<uint8_t>& syndrome_, bool all_in){
    std::vector<EnsembleEntry> ensembles;
    std::vector<uint8_t> md_syndrome = set_the_priorities(syndrome_, all_in);
    uint hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "!!13" << std::endl;
    // std::cout << "md_syndrome before fixing: "<< hw;
    while((std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0) < 3) && (predecoding_edges.size()!=0)){
        qrc::DecodingGraph::Edge* e =predecoding_edges.back();
        md_syndrome[e->detectors.first] = ((uint8_t)1);
        md_syndrome[e->detectors.second] = ((uint8_t)1);
        predecoding_edges.pop_back();
    }
    // hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "md_syndrome after fixing: "<< hw;
    // std::cout << "!!8" << std::endl;
    fp_t predecoding_weights = 0;
    for(auto e: predecoding_edges){
        predecoding_weights += e->edge_weight;
    }
    // std::cout << "!!9" << std::endl;
    std::map<uint, uint> matchings;
    for(auto m: predecoding_edges){
        matchings.insert(std::make_pair(m->detectors.first,m->detectors.second));
        matchings.insert(std::make_pair(m->detectors.second,m->detectors.first));
    }
    // std::cout << "!!10" << std::endl;
    hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    
    uint n_pre_edges = predecoding_edges.size();
    std::vector<std::vector<bool>> chooses;
    
    if(n_pre_edges!=0 && (hw < 5)){
        // This means hw = 4, so we can upair 3 edges and add
        // 6 bits to syndrome
        //std::cout << "!!11!!" << n_pre_edges << ",3) " << std::endl;
        chooses = create_combinations(n_pre_edges,3);
        // std::cout << "!!::" << chooses.size()<<std::endl;
        for(auto c: chooses){
            EnsembleEntry ensemble_entry;
            ensemble_entry.syndrome = md_syndrome;
            ensemble_entry.predecoding_weight = predecoding_weights;
            ensemble_entry.pre_matching = matchings;
            // hw = std::accumulate(ensemble_entry.syndrome.begin(), ensemble_entry.syndrome.begin() + n_detectors, 0);
            // std::cout << "hw before ensemble:" <<hw << "- ";
            for(uint i=0;i<c.size();i++){
                if(c[i]){
                    qrc::DecodingGraph::Edge* e =predecoding_edges[i];
                    ensemble_entry.syndrome[e->detectors.first] = ((uint8_t)1);
                    ensemble_entry.syndrome[e->detectors.second] = ((uint8_t)1);
                    ensemble_entry.predecoding_weight -= e->edge_weight;
                    ensemble_entry.pre_matching.erase(e->detectors.first);
                    ensemble_entry.pre_matching.erase(e->detectors.second);
                }
            }
            // hw = std::accumulate(ensemble_entry.syndrome.begin(), ensemble_entry.syndrome.begin() + n_detectors, 0);
            // std::cout << "hw after ensemble:" <<hw << "- ";
            
            ensembles.push_back(ensemble_entry);
        }

    }
    else if(n_pre_edges!=0 && hw < 7){
        // This means hw = 6, so we can upair 2 edges and add
        // 4 bits to syndrome
        //std::cout << "!!13!!" << n_pre_edges << ",2) " << std::endl;
        chooses = create_combinations(n_pre_edges,2);
        // std::cout << "!!12" << std::endl;
        for(auto c: chooses){
            EnsembleEntry ensemble_entry;
            ensemble_entry.syndrome = md_syndrome;
            ensemble_entry.predecoding_weight = predecoding_weights;
            ensemble_entry.pre_matching = matchings;
            for(uint i=0;i<c.size();i++){
                if(c[i]){
                    qrc::DecodingGraph::Edge* e =predecoding_edges[i];
                    ensemble_entry.syndrome[e->detectors.first] = ((uint8_t)1);
                    ensemble_entry.syndrome[e->detectors.second] = ((uint8_t)1);
                    ensemble_entry.predecoding_weight -= e->edge_weight;
                    ensemble_entry.pre_matching.erase(e->detectors.first);
                    ensemble_entry.pre_matching.erase(e->detectors.second);
                }
            }
            ensembles.push_back(ensemble_entry);
        }

    }
    else if(n_pre_edges!=0 && hw< 9){
        // This means hw = 8, so we can upair 1 edges and add
        // 2 bits to syndrome
        // std::cout << "!!11!!" << n_pre_edges << ",1) " << std::endl;
        chooses = create_combinations(n_pre_edges,1);
        // std::cout << "!!12" << std::endl;
        for(auto c: chooses){
            EnsembleEntry ensemble_entry;
            ensemble_entry.syndrome = md_syndrome;
            ensemble_entry.predecoding_weight = predecoding_weights;
            ensemble_entry.pre_matching = matchings;
            for(uint i=0;i<c.size();i++){
                if(c[i]){
                    qrc::DecodingGraph::Edge* e =predecoding_edges[i];
                    ensemble_entry.syndrome[e->detectors.first] = ((uint8_t)1);
                    ensemble_entry.syndrome[e->detectors.second] = ((uint8_t)1);
                    ensemble_entry.predecoding_weight -= e->edge_weight;
                    ensemble_entry.pre_matching.erase(e->detectors.first);
                    ensemble_entry.pre_matching.erase(e->detectors.second);
                }
            }
            ensembles.push_back(ensemble_entry);
        }

    }
    else if(ensembles.size() == 0){
        // std::cout << "!!12" << std::endl;
        EnsembleEntry ensemble_entry;
        ensemble_entry.syndrome = md_syndrome;
        ensemble_entry.predecoding_weight = predecoding_weights;
        ensemble_entry.pre_matching = matchings;
        ensembles.push_back(ensemble_entry);

    }
    // std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ensembles: "<< ensembles.size()<<", n_pre_edges"<<n_pre_edges<< std::endl;
    return ensembles;

 }

 std::vector<EnsembleEntry>
Predecoder::create_syndromes_ensembles_upair_3(const std::vector<uint8_t>& syndrome_, bool all_in){
    std::vector<EnsembleEntry> ensembles;
    std::vector<uint8_t> md_syndrome = set_the_priorities(syndrome_, all_in);
    uint hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "!!13" << std::endl;
    // std::cout << "md_syndrome before fixing: "<< hw;
    while((std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0) < 3) && (predecoding_edges.size()!=0)){
        qrc::DecodingGraph::Edge* e =predecoding_edges.back();
        md_syndrome[e->detectors.first] = ((uint8_t)1);
        md_syndrome[e->detectors.second] = ((uint8_t)1);
        predecoding_edges.pop_back();
    }
    // hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    
    fp_t predecoding_weights = 0;
    for(auto e: predecoding_edges){
        predecoding_weights += e->edge_weight;
    }
    
    std::map<uint, uint> matchings;
    for(auto m: predecoding_edges){
        matchings.insert(std::make_pair(m->detectors.first,m->detectors.second));
        matchings.insert(std::make_pair(m->detectors.second,m->detectors.first));
    }
    // std::cout << "!!10" << std::endl;
    hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    
    uint n_pre_edges = predecoding_edges.size();
    std::vector<std::vector<bool>> chooses;
    
    if(n_pre_edges!=0){
        // This means hw = 4, so we can upair 3 edges and add
        // 6 bits to syndrome
        //std::cout << "!!11!!" << n_pre_edges << ",3) " << std::endl;
        chooses = create_combinations(n_pre_edges,3);
        // std::cout << "!!::" << chooses.size()<<std::endl;
        for(auto c: chooses){
            EnsembleEntry ensemble_entry;
            ensemble_entry.syndrome = md_syndrome;
            ensemble_entry.predecoding_weight = predecoding_weights;
            ensemble_entry.pre_matching = matchings;
            // hw = std::accumulate(ensemble_entry.syndrome.begin(), ensemble_entry.syndrome.begin() + n_detectors, 0);
            // std::cout << "hw before ensemble:" <<hw << "- ";
            for(uint i=0;i<c.size();i++){
                if(c[i]){
                    qrc::DecodingGraph::Edge* e =predecoding_edges[i];
                    ensemble_entry.syndrome[e->detectors.first] = ((uint8_t)1);
                    ensemble_entry.syndrome[e->detectors.second] = ((uint8_t)1);
                    ensemble_entry.predecoding_weight -= e->edge_weight;
                    ensemble_entry.pre_matching.erase(e->detectors.first);
                    ensemble_entry.pre_matching.erase(e->detectors.second);
                }
            }
            // hw = std::accumulate(ensemble_entry.syndrome.begin(), ensemble_entry.syndrome.begin() + n_detectors, 0);
            // std::cout << "hw after ensemble:" <<hw << "- ";
            
            ensembles.push_back(ensemble_entry);
        }
    }
    return ensembles;

}

std::vector<EnsembleEntry>
Predecoder::create_syndromes_ensembles_unpair_1(const std::vector<uint8_t>& syndrome_, bool all_in){
    std::vector<EnsembleEntry> ensembles;
    std::vector<uint8_t> md_syndrome = set_the_priorities(syndrome_, all_in);
    uint hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "!!13" << std::endl;
    // std::cout << "md_syndrome before fixing: "<< hw;
    while((std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0) < 7) && (predecoding_edges.size()!=0)){
        qrc::DecodingGraph::Edge* e =predecoding_edges.back();
        md_syndrome[e->detectors.first] = ((uint8_t)1);
        md_syndrome[e->detectors.second] = ((uint8_t)1);
        predecoding_edges.pop_back();
    }
    // hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "md_syndrome after fixing: "<< hw;
    // std::cout << "!!8" << std::endl;
    fp_t predecoding_weights = 0;
    for(auto e: predecoding_edges){
        predecoding_weights += e->edge_weight;
    }
    // std::cout << "!!9" << std::endl;
    std::map<uint, uint> matchings;
    for(auto m: predecoding_edges){
        matchings.insert(std::make_pair(m->detectors.first,m->detectors.second));
        matchings.insert(std::make_pair(m->detectors.second,m->detectors.first));
    }
    // std::cout << "!!10" << std::endl;
    hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    
    uint n_pre_edges = predecoding_edges.size();
    std::vector<std::vector<bool>> chooses;
    
    if(n_pre_edges!=0 && hw< 9){
        // This means hw = 8, so we can upair 1 edges and add
        // 2 bits to syndrome
        // std::cout << "!!11!!" << n_pre_edges << ",1) " << std::endl;
        chooses = create_combinations(n_pre_edges,1);
        // std::cout << "!!12" << std::endl;
        for(auto c: chooses){
            EnsembleEntry ensemble_entry;
            ensemble_entry.syndrome = md_syndrome;
            ensemble_entry.predecoding_weight = predecoding_weights;
            ensemble_entry.pre_matching = matchings;
            for(uint i=0;i<c.size();i++){
                if(c[i]){
                    qrc::DecodingGraph::Edge* e =predecoding_edges[i];
                    ensemble_entry.syndrome[e->detectors.first] = ((uint8_t)1);
                    ensemble_entry.syndrome[e->detectors.second] = ((uint8_t)1);
                    ensemble_entry.predecoding_weight -= e->edge_weight;
                    ensemble_entry.pre_matching.erase(e->detectors.first);
                    ensemble_entry.pre_matching.erase(e->detectors.second);
                }
            }
            ensembles.push_back(ensemble_entry);
        }

    }
    else if(ensembles.size() == 0){
        // std::cout << "!!12" << std::endl;
        EnsembleEntry ensemble_entry;
        ensemble_entry.syndrome = md_syndrome;
        ensemble_entry.predecoding_weight = predecoding_weights;
        ensemble_entry.pre_matching = matchings;
        ensembles.push_back(ensemble_entry);

    }
    // std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ensembles: "<< ensembles.size()<<", n_pre_edges"<<n_pre_edges<< std::endl;
    return ensembles;

 }

std::vector<EnsembleEntry>
Predecoder::create_syndromes_ensembles_unpair_2(const std::vector<uint8_t>& syndrome_, bool all_in){
    std::vector<EnsembleEntry> ensembles;
    std::vector<uint8_t> md_syndrome = set_the_priorities(syndrome_, all_in);
    uint hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "!!13" << std::endl;
    // std::cout << "md_syndrome before fixing: "<< hw;
    while((std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0) < 5) && (predecoding_edges.size()!=0)){
        qrc::DecodingGraph::Edge* e =predecoding_edges.back();
        md_syndrome[e->detectors.first] = ((uint8_t)1);
        md_syndrome[e->detectors.second] = ((uint8_t)1);
        predecoding_edges.pop_back();
    }
    // hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "md_syndrome after fixing: "<< hw;
    // std::cout << "!!8" << std::endl;
    fp_t predecoding_weights = 0;
    for(auto e: predecoding_edges){
        predecoding_weights += e->edge_weight;
    }
    // std::cout << "!!9" << std::endl;
    std::map<uint, uint> matchings;
    for(auto m: predecoding_edges){
        matchings.insert(std::make_pair(m->detectors.first,m->detectors.second));
        matchings.insert(std::make_pair(m->detectors.second,m->detectors.first));
    }
    // std::cout << "!!10" << std::endl;
    hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    
    uint n_pre_edges = predecoding_edges.size();
    std::vector<std::vector<bool>> chooses;
    
    
    if(n_pre_edges!=0 && hw < 7){
        // This means hw = 6, so we can upair 2 edges and add
        // 4 bits to syndrome
        //std::cout << "!!13!!" << n_pre_edges << ",2) " << std::endl;
        chooses = create_combinations(n_pre_edges,2);
        // std::cout << "!!12" << std::endl;
        for(auto c: chooses){
            EnsembleEntry ensemble_entry;
            ensemble_entry.syndrome = md_syndrome;
            ensemble_entry.predecoding_weight = predecoding_weights;
            ensemble_entry.pre_matching = matchings;
            for(uint i=0;i<c.size();i++){
                if(c[i]){
                    qrc::DecodingGraph::Edge* e =predecoding_edges[i];
                    ensemble_entry.syndrome[e->detectors.first] = ((uint8_t)1);
                    ensemble_entry.syndrome[e->detectors.second] = ((uint8_t)1);
                    ensemble_entry.predecoding_weight -= e->edge_weight;
                    ensemble_entry.pre_matching.erase(e->detectors.first);
                    ensemble_entry.pre_matching.erase(e->detectors.second);
                }
            }
            ensembles.push_back(ensemble_entry);
        }

    }
    else if(n_pre_edges!=0 && hw< 9){
        // This means hw = 8, so we can upair 1 edges and add
        // 2 bits to syndrome
        // std::cout << "!!11!!" << n_pre_edges << ",1) " << std::endl;
        chooses = create_combinations(n_pre_edges,1);
        // std::cout << "!!12" << std::endl;
        for(auto c: chooses){
            EnsembleEntry ensemble_entry;
            ensemble_entry.syndrome = md_syndrome;
            ensemble_entry.predecoding_weight = predecoding_weights;
            ensemble_entry.pre_matching = matchings;
            for(uint i=0;i<c.size();i++){
                if(c[i]){
                    qrc::DecodingGraph::Edge* e =predecoding_edges[i];
                    ensemble_entry.syndrome[e->detectors.first] = ((uint8_t)1);
                    ensemble_entry.syndrome[e->detectors.second] = ((uint8_t)1);
                    ensemble_entry.predecoding_weight -= e->edge_weight;
                    ensemble_entry.pre_matching.erase(e->detectors.first);
                    ensemble_entry.pre_matching.erase(e->detectors.second);
                }
            }
            ensembles.push_back(ensemble_entry);
        }

    }
    else if(ensembles.size() == 0){
        // std::cout << "!!12" << std::endl;
        EnsembleEntry ensemble_entry;
        ensemble_entry.syndrome = md_syndrome;
        ensemble_entry.predecoding_weight = predecoding_weights;
        ensemble_entry.pre_matching = matchings;
        ensembles.push_back(ensemble_entry);

    }
    // std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ensembles: "<< ensembles.size()<<", n_pre_edges"<<n_pre_edges<< std::endl;
    return ensembles;

}
void 
Predecoder::stage1_decoding(std::vector<uint8_t>& syndrome_, uint round){
    uint hw;
    qrc::DecodingGraph::Edge* e;
    qrc::DecodingGraph::Vertex* v1;
    qrc::DecodingGraph::Vertex* v2;

    
    //Calculate params adaptively based on the number of rounds
    double choosing_param = ((double)0.1/round);

    for (auto it = pf_pairs[0].begin(); it != pf_pairs[0].end();) {

        e = *it;
        v1 = decoding_graph.get_vertex(e->detectors.first);
        v2 = decoding_graph.get_vertex(e->detectors.second);
        if(adaptive_predecoding_map.count(v1->detector) !=0 ||
            adaptive_predecoding_map.count(v2->detector) !=0){
            ++it;
            continue;
        }
        
        if (((v1->feature*choosing_param) <= e->error_probability ) || (v2->feature*choosing_param <= e->error_probability)) {
            predecoding_edges.push_back(e);
            adaptive_predecoding_map[v1->detector] = v2->detector;
            adaptive_predecoding_map[v2->detector] = v1->detector;

            it = pf_pairs[0].erase(it);  // Erase the element from the original vector
            if(v1->detector != BOUNDARY_INDEX){
                syndrome_[v1->detector] = ((uint8_t)0);
            }
            if(v2->detector != BOUNDARY_INDEX){
                syndrome_[v2->detector] = ((uint8_t)0);
            }
            hw = std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0);
            if(hw <= MAX_BF_HW){
                ready = true;
                break;
            }
            
        } else {
            ++it;
        }
    }

}
void 
Predecoder::stage2_decoding(std::vector<uint8_t>& syndrome_, uint& round, bool no_poor){
    uint hw;
    qrc::DecodingGraph::Edge* e;
    qrc::DecodingGraph::Vertex* v1;
    qrc::DecodingGraph::Vertex* v2;
    
    //Calculate params adaptively based on the number of rounds
    double choosing_param = ((double)0.5/round);

    for (auto it = pf_pairs[1].begin(); it != pf_pairs[1].end();) {

        e = *it;
        v1 = decoding_graph.get_vertex(e->detectors.first);
        v2 = decoding_graph.get_vertex(e->detectors.second);
        if(adaptive_predecoding_map.count(v1->detector) !=0 ||
            adaptive_predecoding_map.count(v2->detector) !=0){
            ++it;
            continue;
        }
        if(no_poor){
            if(general){
                if(v1->n_cluster_members % 2 != 0){
                    ++it;
                    continue;
                }
            }
            else if(v1->n_dependents+v2->n_dependents != 1){
                ++it;
                continue;
            }
        }

        if (((v1->feature*choosing_param) <= e->error_probability ) || (v2->feature*choosing_param <= e->error_probability)) {
            predecoding_edges.push_back(e);

            adaptive_predecoding_map[v1->detector] = v2->detector;
            adaptive_predecoding_map[v2->detector] = v1->detector;

            it = pf_pairs[1].erase(it);  // Erase the element from the original vector
            if(v1->detector != BOUNDARY_INDEX){
                syndrome_[v1->detector] = ((uint8_t)0);
            }
            if(v2->detector != BOUNDARY_INDEX){
                syndrome_[v2->detector] = ((uint8_t)0);
            }
            hw = std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0);
            if(hw <= MAX_BF_HW){
                ready = true;
            }
            break;
        } else {
            ++it;
        }
    }

}
void
Predecoder::stage3_decoding(std::vector<uint8_t>& syndrome_, uint& round, bool no_poor){
    uint hw;
    qrc::DecodingGraph::Edge* e;
    qrc::DecodingGraph::Vertex* v1;
    qrc::DecodingGraph::Vertex* v2;
    
    //Calculate params adaptively based on the number of rounds
    double choosing_param = ((double)0.5/round);

    for (auto it = pf_pairs[2].begin(); it != pf_pairs[2].end();) {

        e = *it;
        v1 = decoding_graph.get_vertex(e->detectors.first);
        v2 = decoding_graph.get_vertex(e->detectors.second);
        if(adaptive_predecoding_map.count(v1->detector) !=0 ||
            adaptive_predecoding_map.count(v2->detector) !=0){
            ++it;
            continue;
        }

        if(no_poor){
            if(general){
                if(v1->n_cluster_members % 2 != 0){
                    ++it;
                    continue;
                }
            }
            if(v1->n_dependents+v2->n_dependents != 0){
                ++it;
                continue;
            }
        }

        if (((v1->feature*choosing_param) <= e->error_probability ) || (v2->feature*choosing_param <= e->error_probability)) {
            predecoding_edges.push_back(std::move(*it));

            adaptive_predecoding_map[v1->detector] = v2->detector;
            adaptive_predecoding_map[v2->detector] = v1->detector;

            it = pf_pairs[2].erase(it);  // Erase the element from the original vector
            if(v1->detector != BOUNDARY_INDEX){
                syndrome_[v1->detector] = ((uint8_t)0);
            }
            if(v2->detector != BOUNDARY_INDEX){
                syndrome_[v2->detector] = ((uint8_t)0);
            }
            hw = std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0);
            if(hw <= MAX_BF_HW){
                ready = true;
            }
            break;
        } else {
            ++it;
        }
    }
}

std::vector<std::pair<uint,uint>>
Predecoder::degree_zero_sorted_queue_paths(){

    std::vector<std::pair<uint,uint>> stage4_matching;
    uint n_subgraph_vertices = prematching_subgraph.detector_list.size();
    
    for( uint i = 0; i< degree_zero_vertices.size(); i++){
        auto v = degree_zero_vertices[i];
        if(prematching_subgraph.adjacency_matrix[v].size() != 0){
            std::cout << "SOMTHING IS WRONG! AT STAGE 4 A NONE ZERO DEGREE REMAINED!"<< std::endl;
            continue;
        }
        for( uint j = 0; j< n_subgraph_vertices; j++){
            auto w = decoding_graph.get_vertex(prematching_subgraph.detector_list[j]);
            if(adaptive_predecoding_map.count(w->detector) != 0){
                continue;
            }
            if((v->detector == w->detector)){
                continue;
            }

            stage4_matching.push_back(std::make_pair(v->detector,w->detector));
            

        }
    }

    //std::cout <<"sorting " <<stage4_matching.size() << " items in stage4_matching" << std::endl;

    std::sort(stage4_matching.begin(), stage4_matching.end(),
        [&](const auto& lhs, const auto& rhs) {
            qrc::DecodingGraph::Vertex * lhs_v1 = decoding_graph.get_vertex(lhs.first);
            qrc::DecodingGraph::Vertex * lhs_v2 = decoding_graph.get_vertex(lhs.second);

            qrc::DecodingGraph::Vertex * rhs_v1 = decoding_graph.get_vertex(rhs.first);
            qrc::DecodingGraph::Vertex * rhs_v2 = decoding_graph.get_vertex(rhs.second);

            auto lhs_vpair = std::make_pair(lhs_v1, lhs_v2);
            auto rhs_vpair = std::make_pair(rhs_v1, rhs_v2);

            if(path_table[lhs_vpair].path.size() < path_table[rhs_vpair].path.size()){
                return true;
            }
            else if(path_table[lhs_vpair].path.size() == path_table[rhs_vpair].path.size()){
                if(general){
                    if(path_table[lhs_vpair].distance < path_table[rhs_vpair].distance){
                        return true;
                    }
                }
                else{
                    if(lhs_v1->n_dependents+lhs_v2->n_dependents < rhs_v1->n_dependents+rhs_v2->n_dependents ){
                        return true;
                    }
                    else if(lhs_v1->n_dependents+lhs_v2->n_dependents == rhs_v1->n_dependents+rhs_v2->n_dependents ){
                        uint d_lhs_v1 = prematching_subgraph.adjacency_matrix.count(lhs_v1) == 0 ? 0: 
                        prematching_subgraph.adjacency_matrix[lhs_v1].size();
                        uint d_lhs_v2 = prematching_subgraph.adjacency_matrix.count(lhs_v2) == 0 ? 0: 
                        prematching_subgraph.adjacency_matrix[lhs_v2].size();

                        uint d_rhs_v1 = prematching_subgraph.adjacency_matrix.count(rhs_v1) == 0 ? 0: 
                        prematching_subgraph.adjacency_matrix[rhs_v1].size();
                        uint d_rhs_v2 = prematching_subgraph.adjacency_matrix.count(rhs_v2) == 0 ? 0: 
                        prematching_subgraph.adjacency_matrix[rhs_v2].size();

                        if(d_lhs_v1+d_lhs_v2 < d_rhs_v1+d_rhs_v2 ){
                            return true;
                        }
                    }
                }
            }
            
            return false;
            // Sort in descending order
        }
    );

    // std:: cout << "All possible sorted matching of d0: " << stage4_matching.size() << std::endl;
    // for(const auto& d0: stage4_matching){
    //     qrc::DecodingGraph::Vertex * v1 = decoding_graph.get_vertex(d0.first);
    //     qrc::DecodingGraph::Vertex * v2 = decoding_graph.get_vertex(d0.second);

    //     qrc::DecodingGraph::Vertex * rhs_v1 = decoding_graph.get_vertex(d0.first);
        
    //     auto v1v2_pair = std::make_pair(v1, v2);

    //     std::cout <<"(" <<d0.first << "(" << v1->n_dependents<< ")" <<","<<d0.second << "(" << v2->n_dependents<< ")" <<")"<< ": "<< path_table[v1v2_pair].distance << std::endl;
    // }

    return stage4_matching;
}

void
Predecoder::stage0_decoding(std::vector<uint8_t>& syndrome_, uint& round, bool no_poor){
    uint hw = std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0);
    std::map<uint,uint> matching;

    // std::cout << "Finding degree 0 nodes ... "<< std::endl;
    std::vector<std::pair<uint,uint>> degree_zero_queue = degree_zero_sorted_queue_paths();
    // std::cout << "Done with degree_zero_queue sorting " << degree_zero_queue.size() << std::endl;

    // if(degree_zero_queue.size() == 0){
    //     //std::cout << "SOMETHING WRONG HAPPENED! WENT TO STAGE 4 WITHOUT ANY DEGREE 0 NODE";
    // }
    number_paths_singleton_matching_candidate = degree_zero_queue.size();
    uint loop_count = 0;
    // std::cout << "Starting the loop - hw = " << hw << std::endl;
    for(auto it = degree_zero_queue.begin(); it != degree_zero_queue.end();) {
        hw = std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0); 
        if(hw <= MAX_BF_HW){
            ready = true;
            break;
        }

        loop_count++;
        // std::cout << "_____Stage 4: round " << loop_count << std::endl;
        std::pair<uint,uint> m = *it;

        // std::cout <<"Round "<<loop_count<< ": Fetched "<< m.first << 
        // " and " << m.second << std::endl;
        if(no_poor){
            auto v1 = decoding_graph.get_vertex(m.first);
            auto v2 = decoding_graph.get_vertex(m.second);
            if(general){
                if((v1->n_cluster_members + v2->n_cluster_members)%2 != 0){
                    ++it;
                    continue;
                }
            }
            else if(v1->n_dependents+v2->n_dependents !=0){
                ++it;
                continue;
            }
        }

        if(adaptive_predecoding_map.count(m.first) != 0 || 
           adaptive_predecoding_map.count(m.second) != 0){
            ++it;
            continue;
        }     
        // std::cout << "Found best matching between " << m.first << 
        // " and " << m.second << std::endl;

        if(m.first != BOUNDARY_INDEX){
            syndrome_[m.first] = ((uint8_t)0);
        }
        if(m.second != BOUNDARY_INDEX){
            syndrome_[m.second] = ((uint8_t)0);
        }

        adaptive_predecoding_map[m.first] = m.second;
        adaptive_predecoding_map[m.second] = m.first;
       
    }

}

qrc::DecoderShotResult 
Predecoder::adaptively_decode_error(const std::vector<uint8_t>& syndrome_){
    // std::cout << "START adaptively_decode_error";

    std::vector<uint8_t> syndrome_copy(syndrome_.begin(), syndrome_.end());

    qrc::DecoderShotResult decoding_result;

    uint decoding_round = 0;

    ////Simulating the clock
    reset_clocks();

    ////Simulating the clock


    predecoding_edges.clear();
    adaptive_predecoding_map.clear();
    number_of_edges = 0;
    reached_stage.fill(false);
    uint hw = std::accumulate(syndrome_copy.begin(), syndrome_copy.begin()+n_detectors,0); 
    uint previous_hw = hw;
    ready =  hw <= MAX_BF10_HW ? true : false;


    while(!ready){
        hw = std::accumulate(syndrome_copy.begin(), syndrome_copy.begin()+n_detectors,0 );
        uint n_edges_matched = ((uint)ceil(((fp_t)(previous_hw-hw)/2)));
        previous_hw = hw;
        // total_cycles += n_edges_matched;

        
        // total_cycles_parallel_updating += n_edges_matched;
        // std::cout << "___________Adaptive Predecoding MATCHING__________" << std::endl;
        // for(const auto& m: adaptive_predecoding_map){
        //         std::cout << "(" << m.first << ", " << m.second << ")-";   
        // }   
        // std::cout << std::endl<< "_____________________________" << std::endl;
        decoding_round ++;
        number_of_rounds++;
        // std::cout << "===================== Start of Round " << decoding_round << " hw = "<< hw <<" ===============" << std::endl;
        // std::cout << "total_cycles_parallel_updating: " << total_cycles_parallel_updating << " check_condition: " << check_condition << " n_edges_matched: " << n_edges_matched<< std::endl;
        // std::cout << "Set the Prioritization..." << std::endl;

        prioritize_edges_by_degree(syndrome_copy);
        // if(decoding_round == 1){
        //     starting_number_of_edges = number_of_edges;
        // }
        // std::cout << "decoding_round: "<< decoding_round << "number_of_edges: " <<number_of_edges<< std::endl;
        if(decoding_round != 1){
            // starting_number_of_edges = number_of_edges;
            total_cycles += 3;
            total_cycles_parallel_updating += 3;
        }
        total_cycles += 3 + number_of_edges;
        total_cycles_parallel_updating += 4 + ((uint)ceil(((fp_t)(number_of_edges)/4)));;
        if(total_cycles_parallel_updating > MAX_BF10_CYCLE){
            MAX_BF_HW = MAX_BF8_HW;
        }
        if(total_cycles_parallel_updating > MAX_BF8_CYCLE){
            MAX_BF_HW = MAX_BF6_HW;
        }


        //stage 1
        // std::cout << "Stage 1 ..." << std::endl;
        // This stage decodes stand alone pairs
        stage1_decoding(syndrome_copy, decoding_round);
        if(reached_stage[1]  == false &&
            reached_stage[2]  == false &&
            reached_stage[3]  == false &&
            reached_stage[4]  == false &&
            reached_stage[5]  == false){
                    
            reached_stage[0]  = true;
        }

        if(ready){
            break;
        }

        

        if(exist_degree_one_no_poor){
            // stage 2.1 
            // std::cout << "Stage 2.1 No Poor..." << std::endl;
            // This stage with true no poor flag to just 
            // match min degree 1 if there will be no stand alon 
            // node
            stage2_decoding(syndrome_copy, decoding_round, true);
            if(reached_stage[2]  == false &&
                reached_stage[3]  == false &&
                reached_stage[4]  == false  &&
                reached_stage[5]  == false){
                    reached_stage[0]  = false;
                    reached_stage[1]  = true;
                }
            
            // std::cout << "Done with Stage 2.1 No Poor..." << std::endl;
            continue;
        }
        // if(stage1_queue_size!= 0){
        //     continue;
        // }

        // stage 2.2
        // std::cout << "Stage 3 No Poor..." << std::endl;
        if(exist_multi_degree_no_poor){
            stage3_decoding(syndrome_copy, decoding_round, true);
            if(reached_stage[3]  == false &&
                reached_stage[4]  == false  &&
                reached_stage[5]  == false){
                reached_stage[0] = false;
                reached_stage[1] = false;
                reached_stage[2] = true;
            }
            continue;
        }
        uint prior_matching_size = adaptive_predecoding_map.size();

        if(degree_zero_vertices.size()!=0){
            // stage 0 
            // std::cout << "Stage 0-Degree No Poor..." << std::endl;
            stage0_decoding(syndrome_copy, decoding_round, true);
            total_cycles += std::max(number_of_edges, number_paths_singleton_matching_candidate);
            total_cycles -= number_of_edges;

            total_cycles_parallel_updating += std::max(number_of_edges, number_paths_singleton_matching_candidate);
            total_cycles_parallel_updating -= number_of_edges;
            // std::cout << "Done with Stage 0-Degree No Poor..." << std::endl;
            if(reached_stage[4]  == false   &&
                reached_stage[5]  == false ){
                reached_stage[0] = false;
                reached_stage[1] = false;
                reached_stage[2] = false;
                reached_stage[3] = true;
            }
            uint posterior_matching_size = adaptive_predecoding_map.size();
            if(prior_matching_size!= posterior_matching_size){
                continue;
            }
            
            
        }

        // stage 2 
        // std::cout << "Stage 2 Regular..." << std::endl;

        if(pf_pairs[1].size()!= 0){
            stage2_decoding(syndrome_copy, decoding_round, false);
            if(reached_stage[5]  == false){
                reached_stage[0] = false;
                reached_stage[1] = false;
                reached_stage[2] = false;
                reached_stage[3]  = false;
                reached_stage[4]  = true;
            }
                continue;
        }

        // stage 3 
        // std::cout << "Stage 3 Regular..." << std::endl;
        stage3_decoding(syndrome_copy, decoding_round, false);
        reached_stage[0] = false;
        reached_stage[1] = false;
        reached_stage[2] = false;
        reached_stage[3] = false;
        reached_stage[4] = false;
        reached_stage[5] = true;

    }
    
    // total_cycles += 8;
    // total_cycles_parallel_updating += 8;
    // std:: cout << "All stages completed." << std::endl;

    if(std::accumulate(syndrome_copy.begin(), syndrome_copy.begin() + n_detectors, 0) > MAX_BF_HW ){
        std::cout << "SOMETHING WRONG HAPPENED! HW MORE THAN " << MAX_BF_HW << std::endl;
    }

    decoding_result = decoder->decode_error(syndrome_copy);

    // std:: cout << "Decoded hw = " << MAX_BF_HW << std::endl;

    // for(auto edge : predecoding_edges){
    //     if(decoding_result.matching.count(edge->detectors.first) || decoding_result.matching.count(edge->detectors.second)
    //      || adaptive_predecoding_map.count(edge->detectors.first) || adaptive_predecoding_map.count(edge->detectors.second) ){
    //         std::cout << "SOMETHING IS WRONG! DUPLICATED MATCHING";
    //     }

    //     adaptive_predecoding_map[edge->detectors.first] = edge->detectors.second;
    //     adaptive_predecoding_map[edge->detectors.second] = edge->detectors.first;

    // }

    decoding_result.matching.insert(adaptive_predecoding_map.begin(),adaptive_predecoding_map.end());

    // std::cout << "Attached all the matchings" << std::endl;

    decoding_result.correction = get_correction_from_matching(decoding_result.matching);
    decoding_result.is_logical_error = 
        qrc::is_logical_error(decoding_result.correction, syndrome_, n_detectors, n_observables);
    decoding_result.weight = calc_matching_weight(decoding_result.matching);

    return decoding_result;


}

void 
Predecoder::reset_clocks(){
    MAX_BF_HW = MAX_BF10_HW;
    total_cycles = 0;
    total_cycles_parallel_updating = 0;
    updating_clk_parallel = 0;
    updating_clk = 0;
    number_of_rounds = 0;
    check_condition = 0;
    number_paths_singleton_matching_candidate = 0;
}

    
std::vector<EnsembleEntry> 
Predecoder::create_syndromes_no_unpairing(const std::vector<uint8_t>& syndrome_, bool all_in){
    std::vector<EnsembleEntry> ensembles;
    std::vector<uint8_t> md_syndrome = set_the_priorities(syndrome_, all_in);
    // std::cout << "!!13" << std::endl;
    // std::cout << "md_syndrome before fixing: "<< hw;
    while((std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0)< MAX_BF_HW-1) && (predecoding_edges.size()!=0)){
        qrc::DecodingGraph::Edge* e =predecoding_edges.back();
        md_syndrome[e->detectors.first] = ((uint8_t)1);
        md_syndrome[e->detectors.second] = ((uint8_t)1);
        predecoding_edges.pop_back();
    }
    // hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "md_syndrome after fixing: "<< hw;
    // std::cout << "!!8" << std::endl;
    fp_t predecoding_weights = 0;
    for(auto e: predecoding_edges){
        predecoding_weights += e->edge_weight;
    }
    // std::cout << "!!9" << std::endl;
    std::map<uint, uint> matchings;
    for(auto m: predecoding_edges){
        matchings.insert(std::make_pair(m->detectors.first,m->detectors.second));
        matchings.insert(std::make_pair(m->detectors.second,m->detectors.first));
    }
    uint n_pre_edges = predecoding_edges.size();
    
    // std::cout << "!!12" << std::endl;
    EnsembleEntry ensemble_entry;
    ensemble_entry.syndrome = md_syndrome;
    ensemble_entry.predecoding_weight = predecoding_weights;
    ensemble_entry.pre_matching = matchings;
    ensembles.push_back(ensemble_entry);

    // std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ensembles: "<< ensembles.size()<<", n_pre_edges"<<n_pre_edges<< std::endl;
    return ensembles;
    
}

std::vector<EnsembleEntry> 
Predecoder::create_syndromes_selected_unpairing(const std::vector<uint8_t>& syndrome_){
    std::vector<EnsembleEntry> ensembles;
    std::vector<uint8_t> md_syndrome = set_the_priorities(syndrome_);
    uint hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "!!13" << std::endl;
    // std::cout << "md_syndrome before fixing: "<< hw;
    while((hw < MAX_BF_HW-1) && (predecoding_edges.size()!=0)){
        qrc::DecodingGraph::Edge* e =predecoding_edges.back();
        md_syndrome[e->detectors.first] = ((uint8_t)1);
        md_syndrome[e->detectors.second] = ((uint8_t)1);
        predecoding_edges.pop_back();
    }
    // hw = std::accumulate(md_syndrome.begin(), md_syndrome.begin() + n_detectors, 0);
    // std::cout << "md_syndrome after fixing: "<< hw;
    // std::cout << "!!8" << std::endl;
    fp_t predecoding_weights = 0;
    for(auto e: predecoding_edges){
        predecoding_weights += e->edge_weight;
    }
    // std::cout << "!!9" << std::endl;
    std::map<uint, uint> matchings;
    for(auto m: predecoding_edges){
        matchings.insert(std::make_pair(m->detectors.first,m->detectors.second));
        matchings.insert(std::make_pair(m->detectors.second,m->detectors.first));
    }
    uint n_pre_edges = predecoding_edges.size();
    
    // std::cout << "!!12" << std::endl;
    EnsembleEntry ensemble_entry;
    ensemble_entry.syndrome = md_syndrome;
    ensemble_entry.predecoding_weight = predecoding_weights;
    ensemble_entry.pre_matching = matchings;
    ensembles.push_back(ensemble_entry);

    // std::cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ ensembles: "<< ensembles.size()<<", n_pre_edges"<<n_pre_edges<< std::endl;
    return ensembles;
    
}

qrc::DecoderShotResult 
Predecoder::ensemble_decode_error(const std::vector<uint8_t>& syndrome_, uint max_unpair_number, bool all_in){


    std::vector<EnsembleEntry> ensembles;
    if(max_unpair_number == 3){
        ensembles = create_syndromes_ensembles(syndrome_, all_in);//create_syndromes_ensembles(syndrome_, all_in);
    }
    else if(max_unpair_number == 2){
        ensembles = create_syndromes_ensembles_unpair_2(syndrome_, all_in);
    }
    else if(max_unpair_number == 1){
        ensembles = create_syndromes_ensembles_unpair_1(syndrome_, all_in);
    }
    else if(max_unpair_number == 0){
        ensembles = create_syndromes_no_unpairing(syndrome_, all_in);
    }

    fp_t min_weight = std::numeric_limits<double>::max();
    int solution_index = -1;
    qrc::DecoderShotResult decoding_result;
    qrc::DecoderShotResult solution_shot_res;
    for(int i = 0; i< ensembles.size(); i++){
        if(10<std::accumulate(ensembles[i].syndrome.begin(), ensembles[i].syndrome.begin()+n_detectors, 0) ){
            std::cout << "Too high" << std::endl;
            continue;
        }
        decoding_result = decoder->decode_error(ensembles[i].syndrome);
        decoding_result.weight = calc_matching_weight(decoding_result.matching);

        

        if((ensembles[i].predecoding_weight+decoding_result.weight )< min_weight){
            min_weight = ensembles[i].predecoding_weight+decoding_result.weight;
            solution_index = i;
            solution_shot_res = decoding_result;
        }
    }
   if(solution_index== -1){
        std::cout << "HW more than 10 - ensemble size: " << ensembles.size()<<std::endl;
        // uint hw = std::accumulate(syndrome_.begin(), syndrome_.begin()+n_detectors, 0);
        // std::cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<< pf_groups[0].size() << "," <<hw <<std::endl;
        // for(auto x : pf_groups[0]){
        //     std::cout << x->detector << ", ";
        // }
        // std::cout << std::endl<< "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<< std::endl;
        solution_shot_res.is_logical_error = true;
        return solution_shot_res;
    }


    std::map<uint, uint> temp = ensembles[solution_index].pre_matching;
    // Attaching two matchings
    solution_shot_res.matching.insert(ensembles[solution_index].pre_matching.begin(),ensembles[solution_index].pre_matching.end());

    
    
    // Fixing the correction for the soltion 
    solution_shot_res.correction = get_correction_from_matching(solution_shot_res.matching);
    // std::cout << "!!4" << std::endl;
    // Determining if it is a logical error or not
    solution_shot_res.is_logical_error = 
        qrc::is_logical_error(solution_shot_res.correction, syndrome_, n_detectors, n_observables);
    
    // Assigning the weight
    solution_shot_res.weight = min_weight;
    return solution_shot_res;

}

fp_t 
Predecoder::calc_matching_weight(std::map<uint,uint> matching_){
    std::map<uint, uint> visited_pairs;
    fp_t match_weight = 0;
    for(auto m: matching_){
        uint first = m.first;
        uint second = m.second;
        if(visited_pairs.count(first) == 0 && visited_pairs.count(second) == 0){
            visited_pairs[first] = second;

            auto v = decoding_graph.get_vertex(m.first);
            auto w = decoding_graph.get_vertex(m.second);
            auto pair = std::make_pair(v,w);
            match_weight += path_table[pair].distance;

        }
    }
    return match_weight;
}



void
Predecoder::prioritize_edges_by_degree(const std::vector<uint8_t>& syndrome_, bool all_in) {
    // std::cout << "start prioritize_edges_by_degree";
    std::vector<uint8_t> md_syndrome(syndrome_.begin(), syndrome_.end());

    // subgraph is made inside get_fast_matching_group
    std::vector<qrc::DecodingGraph::Vertex*> fast_group = get_fast_matching_group_all_neighbors(syndrome_);
    auto subgraph = get_adjacency_matrix();

    if(general){
        /// find the connected components and set the 
        auto components = qpd::find_connected_components(subgraph.adjacency_matrix);
        for(const auto& comp : components){
            for(const auto& v_c : comp){
                v_c->n_cluster_members = comp.size();
            }
        }
    }
    exist_multi_degree_no_poor = false;
    exist_degree_one_no_poor = false;

    uint edge_count = 0;

    std::map<uint, uint> degree_map;
    for (const auto& entry : subgraph.adjacency_matrix) {
        degree_map[entry.first->detector] =((uint)(entry.second.size()));
    }

    // Sort the edges based on the provided criteria
    uint sort_count =0;
    custom_sort(sorted_predecoding_edges, degree_map, decoding_graph, subgraph.adjacency_matrix);   


    //creates pf_pairs
    pf_pairs.clear();
    pf_groups.clear();
    pf_pairs.resize(3);
    pf_groups.clear();
    pf_groups.resize(3);
    for(const auto& sorted_edge : sorted_predecoding_edges){
        // if(sorted_edge->detectors.first == BOUNDARY_INDEX ||sorted_edge->detectors.second == BOUNDARY_INDEX){
        //     continue;
        // }
        bool e1_degree_1 = (degree_map.at(sorted_edge->detectors.first) == 1 && degree_map.at(sorted_edge->detectors.second) == 1);
        if(e1_degree_1){
            pf_pairs[0].push_back(sorted_edge);
        }
        else if(std::min(degree_map.at(sorted_edge->detectors.first), degree_map.at(sorted_edge->detectors.second)) == 1){
            pf_pairs[1].push_back(sorted_edge);
            auto v1 = decoding_graph.get_vertex(sorted_edge->detectors.first);
            auto v2 = decoding_graph.get_vertex(sorted_edge->detectors.second);
            if(general){
                if(v1->n_cluster_members % 2 == 0){
                    exist_degree_one_no_poor =  true;
                }
            }
            else if(v1->n_dependents+v2->n_dependents == 1){
                exist_degree_one_no_poor =  true;
            }
        }
        else{
            pf_pairs[2].push_back(sorted_edge);
            auto v1 = decoding_graph.get_vertex(sorted_edge->detectors.first);
            auto v2 = decoding_graph.get_vertex(sorted_edge->detectors.second);
            if(general){
                if(v1->n_cluster_members % 2 == 0){
                    exist_multi_degree_no_poor =  true;
                }
            }
            else if(v1->n_dependents+v2->n_dependents == 0){
                exist_multi_degree_no_poor =  true;
            }
        }
    } 

    // std::cout << "____PF edges____" << std::endl;
    // edge_count = 0;
    // for(int r = 0; r<degree_zero_vertices.size(); r++){
    //     std::cout << "Node " << degree_zero_vertices[r]->detector << std::endl;
    // }
    // for(int i =0 ; i < pf_pairs.size(); i++){
    //     std::cout <<" Group = " <<i+1 << std::endl;
    //     for(auto se : pf_pairs[i]){
    //         edge_count++;
    //         std::cout << edge_count<< ": ";
    //         qrc::DecodingGraph::Vertex* vertex1 = decoding_graph.get_vertex(se->detectors.first);
    //         qrc::DecodingGraph::Vertex* vertex2 = decoding_graph.get_vertex(se->detectors.second);
            
    //         size_t degree1 = subgraph.adjacency_matrix.at(vertex1).size();
    //         size_t degree2 = subgraph.adjacency_matrix.at(vertex2).size();
    //         std::cout << se->detectors.first <<"(" << degree1 << " - " <<  vertex1->feature<< " - " << vertex1->n_dependents << ")"<< ", " <<
    //         se->detectors.second <<"(" << degree2 <<"- " <<  vertex2->feature<< " - " << vertex2->n_dependents <<  ") - " << se->error_probability << std::endl;
            
    //     }
    // }
    // std::cout << "____PF edges____" << std::endl; 


    
}


/*
This part is for simulation on hardware of Promatch
*/


void 
Predecoder::update_degrees(){
     /* 
    This function put all the nodes that has at least one adjacent flipped
    node to the fast group.
    */
    //update the prematching subgraph (the subgraph that only has flipped 
    // parity bits.)
    update_adjacency_matrix_and_detector_list(syndrome_register, 1, true);
    sorted_predecoding_edges.clear();
    degree_zero_vertices.clear();

    //auto subgraph = get_adjacency_matrix();
    uint n_vertices = prematching_subgraph.detector_list.size();

    qrc::DecodingGraph::Edge* e;

    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = prematching_subgraph.detector_list[vj];
        
        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);


        //This if exclude the cases that degree is one because of boundary vertex
        if (prematching_subgraph.adjacency_matrix[vdj].size() != 0){
            
            std::vector<qrc::DecodingGraph::Vertex*> neighbors = prematching_subgraph.adjacency_matrix[vdj];
            
            // Checking if the node with degree one is connected to a node with degree one.
                //match two detectors of error chain length one
                
            for(const auto& n: neighbors){
                e = decoding_graph.get_edge(vdj,n);
                if(e == nullptr){
                    std::cout <<"THIS MESSAGE SHOULD NOT BE PRINTED. NULLPTR FOR EDGE OF ADJACENT NODES." << std::endl;
                }
                addUniqueEdge(sorted_predecoding_edges, e);
                if(prematching_subgraph.adjacency_matrix[n].size() ==1){
                    vdj->n_dependents++;
                }
                 
            }
        
        }
        else{
            degree_zero_vertices.push_back(vdj);
        }
    }
}

void 
Predecoder::set_simulator_adjacenct_matrix(){
    std::vector<qrc::DecodingGraph::Vertex*> list_of_vertices = decoding_graph.vertices();
    uint n_vertices = list_of_vertices.size();
    simulator_adjacenct_matrix.assign(n_vertices, std::vector<uint8_t>(n_vertices, 0));

    for(uint i = 0; i<n_vertices; i++){
        for(uint j = i+1; j<n_vertices; j++){
            auto e = decoding_graph.get_edge(list_of_vertices[i],
                        list_of_vertices[j]);
            if(e != nullptr){
                uint first_indx = list_of_vertices[i]->detector ==  BOUNDARY_INDEX?n_vertices-1:list_of_vertices[i]->detector;
                uint second_indx = list_of_vertices[j]->detector ==  BOUNDARY_INDEX?n_vertices-1:list_of_vertices[j]->detector;

                simulator_adjacenct_matrix[i][j] = 1;
                simulator_adjacenct_matrix[j][i] = 1;
            }
            
        }
    }
    
} 



PredecoderShotResult_v2 
Predecoder::smith_predecoder(const std::vector<uint8_t>& syndrome_, bool add_boundry ){
    PredecoderShotResult_v2 result = {};

    result.post_syndrome = syndrome_;

    update_adjacency_matrix_and_detector_list(syndrome_, 1, add_boundry);

    auto subgraph = get_adjacency_matrix();
    uint n_vertices = subgraph.detector_list.size();

    if(subgraph.distance_with_neighbor != 1){
        std::cout << "Cannot use the current adjacent matrix becuase it is not built with distance 1" << std::endl;
        return result;
    }

    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = subgraph.detector_list[vj];
        
        if(dj == BOUNDARY_INDEX){
            continue;
        }

        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);
        

        std::vector<qrc::DecodingGraph::Vertex*> neighbors = subgraph.adjacency_matrix[vdj];

        // Checking if the node with degree one is connected to a node with degree one.
        for(auto const& vdi : neighbors){

            if(vdi->detector == BOUNDARY_INDEX){
                continue;
            }

            // TO check if the pair alreay exit in the pre_matches vector.
            bool found = false;
            for (const auto& p : result.pre_matches) {
                if ((p.first == vdj->detector && p.second == vdi->detector) || 
                    (p.first == vdi->detector && p.second == vdj->detector)){
                    found = true;
                    break;
                }
            }
            if(found){
                continue;
            }

            //counting the edges and making sure that the min detector will be the first element of the pair.
            result.edge_counting[std::make_pair(std::min(vdj->detector, vdi->detector), std::max(vdj->detector, vdi->detector))] ++;
            
            result.pre_matches.push_back(std::make_pair(std::min(vdj->detector, vdi->detector), std::max(vdj->detector, vdi->detector)));

            if(vdj->detector!=BOUNDARY_INDEX)
                result.post_syndrome[vdj->detector]=!(result.post_syndrome[vdj->detector]);
            if(vdi->detector!=BOUNDARY_INDEX)
                result.post_syndrome[vdi->detector]=!(result.post_syndrome[vdi->detector]);
        }
    }

    return result;

}

PredecoderDecoderShotResult 
Predecoder::smith_decode_error(const std::vector<uint8_t>& syndrome_){
    PredecoderDecoderShotResult decoding_result;
    if(std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0) < 11){
        qrc::DecoderShotResult r = decoder->decode_error(syndrome_);
        
        decoding_result.correction = r.correction; 
        for(const auto& x : r.matching){
            decoding_result.dec_matching.push_back(x);
        }
        
        decoding_result.weight = calc_matching_weight(r.matching);

        decoding_result.is_logical_error = 
            qrc::is_logical_error(decoding_result.correction, syndrome_, n_detectors, n_observables);
    
        return decoding_result;
    }
    PredecoderShotResult_v2 predecode_results = smith_predecoder(syndrome_, true);

    // if(std::accumulate(predecode_results.post_syndrome.begin(), predecode_results.post_syndrome.begin() + n_detectors, 0) > 10 ){
    //     std::cout << "HW MORE THAN 10" << std::endl;
    //     auto flipped_bits = qpd::syndrome_compressed(syndrome_);
    // }

    if(decoder->name() == "MWPMDecoder" && std::accumulate(predecode_results.post_syndrome.begin(), predecode_results.post_syndrome.begin() + n_detectors, 0) > 10){
        decoding_result.weight = DBL_MAX;
        decoding_result.is_logical_error = true;
        return decoding_result;

    }

    qrc::DecoderShotResult r = decoder->decode_error(predecode_results.post_syndrome);


    
    decoding_result.pre_matching = predecode_results.pre_matches;
    decoding_result.correcting_edge = predecode_results.edge_counting;
    decoding_result.weight = 0;


    //adding the matchings from the decoder to the 
    std::set<uint> visited;
    visited.clear();
    for (const auto& di_dj : r.matching) {
        uint di = di_dj.first;
        uint dj = di_dj.second;

        if (visited.count(di) || visited.count(dj)) {
            continue;
        }
        decoding_result.dec_matching.push_back(std::make_pair(std::min(dj, di), std::max(dj, di)));
        qrc::DecodingGraph::Vertex * vdi = graph.get_vertex(di);
        qrc::DecodingGraph::Vertex * vdj = graph.get_vertex(dj);
        auto vdi_vdj = std::make_pair(vdi, vdj);
        std::vector<qrc::DecodingGraph::Vertex*> detector_path(path_table[vdi_vdj].path);
        for (uint i = 1; i < detector_path.size(); i++) {
            // Get edge from decoding graph.
            auto wi = detector_path[i-1];
            auto wj = detector_path[i];
            decoding_result.correcting_edge[std::make_pair(std::min(wi->detector, wj->detector), std::max(wi->detector, wj->detector))]++;
            // The edge should exist.
            
        }
        visited.insert(di);
        visited.insert(dj);
        
    }

    // Finding the correction based on the counts of edges
    decoding_result.correction.clear();
    decoding_result.correction.resize(n_observables, 0);


    for(auto const& pair : decoding_result.correcting_edge){
        uint count = pair.second;
        std::pair<uint,uint> v1_v2 = pair.first;
        auto edge = graph.get_edge(v1_v2.first, v1_v2.second);
        if(edge == nullptr){
            std::cout << "SOMETHING IS WRONG! EDGE IS NOT FOUND!" << std::endl;
            return decoding_result;
        }

        decoding_result.weight += (count*edge->edge_weight);

        for(uint c =0; c<count; c++){
            for (uint obs : edge->frames) {
                // Flip the bit.
                if (obs >= 0) {
                    decoding_result.correction[obs] = !decoding_result.correction[obs];
                }
            }

        }
    }


    decoding_result.is_logical_error = 
        qrc::is_logical_error(decoding_result.correction, syndrome_, n_detectors, n_observables);

    return decoding_result;
    

}



PredecoderShotResult_v2 
Predecoder::clique_predecoder(const std::vector<uint8_t>& syndrome_, bool add_boundry ){
    PredecoderShotResult_v2 result = {};

    // contains detectors that are not matched during the predecoder
    std::vector<uint8_t> post_podecoder_syndrome;

    result.post_syndrome = syndrome_;
    post_podecoder_syndrome = syndrome_;

    update_adjacency_matrix_and_detector_list(syndrome_, 1, add_boundry);

    auto subgraph = get_adjacency_matrix();
    uint n_vertices = subgraph.detector_list.size();

    if(subgraph.distance_with_neighbor != 1){
        std::cout << "Cannot use the current adjacent matrix becuase it is not built with distance 1" << std::endl;
        return result;
    }

    for (uint vj = 0; vj < n_vertices; vj++) {
        uint dj = subgraph.detector_list[vj];
        
        if(dj == BOUNDARY_INDEX){
            continue;
        }

        qrc::DecodingGraph::Vertex * vdj = decoding_graph.get_vertex(dj);
        

        std::vector<qrc::DecodingGraph::Vertex*> neighbors = subgraph.adjacency_matrix[vdj];

        if(neighbors.size() % 2 != 1){
            result.edge_counting.clear();
            result.pre_matches.clear();
            return result;
        }

        // Checking if the node with degree one is connected to a node with degree one.
        for(auto const& vdi : neighbors){

            if(vdi->detector == BOUNDARY_INDEX){
                continue;
            }

            // TO check if the pair alreay exit in the pre_matches vector.
            bool found = false;
            for (const auto& p : result.pre_matches) {
                if ((p.first == vdj->detector && p.second == vdi->detector) || 
                    (p.first == vdi->detector && p.second == vdj->detector)){
                    found = true;
                    break;
                }
            }
            if(found){
                continue;
            }

            //counting the edges and making sure that the min detector will be the first element of the pair.
            result.edge_counting[std::make_pair(std::min(vdj->detector, vdi->detector), std::max(vdj->detector, vdi->detector))] ++;
            
            result.pre_matches.push_back(std::make_pair(std::min(vdj->detector, vdi->detector), std::max(vdj->detector, vdi->detector)));

            if(vdj->detector!=BOUNDARY_INDEX)
                post_podecoder_syndrome[vdj->detector]=((uint8_t)0);
            if(vdi->detector!=BOUNDARY_INDEX)
                post_podecoder_syndrome[vdi->detector]=((uint8_t)0);
        }
    }
    result.post_syndrome = post_podecoder_syndrome;
    
    return result;

}

PredecoderDecoderShotResult 
Predecoder::clique_decode_error(const std::vector<uint8_t>& syndrome_){
    PredecoderDecoderShotResult decoding_result;
    if(std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0) < 11){
        qrc::DecoderShotResult r = decoder->decode_error(syndrome_);
        
        decoding_result.correction = r.correction; 
        for(const auto& x : r.matching){
            decoding_result.dec_matching.push_back(x);
        }
        
        decoding_result.weight = calc_matching_weight(r.matching);

        decoding_result.is_logical_error = 
            qrc::is_logical_error(decoding_result.correction, syndrome_, n_detectors, n_observables);
    
        return decoding_result;
    }
    PredecoderShotResult_v2 predecode_results = clique_predecoder(syndrome_, false);
    
    // if(std::accumulate(predecode_results.post_syndrome.begin(), predecode_results.post_syndrome.begin() + n_detectors, 0) > 10 ){
    //     std::cout << "HW MORE THAN 10" << std::endl;
    //     auto flipped_bits = qpd::syndrome_compressed(syndrome_);
    // }

    if (predecode_results.pre_matches.size() != 0){
        
        decoding_result.pre_matching = predecode_results.pre_matches;
        decoding_result.correcting_edge = predecode_results.edge_counting;
        decoding_result.weight = 0;


        

        // Finding the correction based on the counts of edges
        decoding_result.correction.clear();
        decoding_result.correction.resize(n_observables, 0);


        for(auto const& pair : decoding_result.correcting_edge){
            uint count = pair.second;
            std::pair<uint,uint> v1_v2 = pair.first;
            auto edge = graph.get_edge(v1_v2.first, v1_v2.second);
            if(edge == nullptr){
                std::cout << "SOMETHING IS WRONG! EDGE IS NOT FOUND!" << std::endl;
                return decoding_result;
            }

            decoding_result.weight += (count*edge->edge_weight);

            for(uint c =0; c<count; c++){
                for (uint obs : edge->frames) {
                    // Flip the bit.
                    if (obs >= 0) {
                        decoding_result.correction[obs] = !decoding_result.correction[obs];
                    }
                }

            }
        }

    }
    else if(decoder->name() == "MWPMDecoder" && std::accumulate(syndrome_.begin(), syndrome_.begin() + n_detectors, 0) > 10){
        decoding_result.weight = DBL_MAX;
        decoding_result.is_logical_error = true;
        return decoding_result;

    }
    else{
        qrc::DecoderShotResult r = decoder->decode_error(syndrome_);
        
        decoding_result.correction = r.correction; 
        for(const auto& x : r.matching){
            decoding_result.dec_matching.push_back(x);
        }
        
        decoding_result.weight = calc_matching_weight(r.matching);
        
    }
    
    decoding_result.is_logical_error = 
            qrc::is_logical_error(decoding_result.correction, syndrome_, n_detectors, n_observables);
    
    return decoding_result;
}

uint
Predecoder::print_paths_matching(std::map<uint, uint> matching){
    std::set<uint> visited;
    uint total_path_d = 0;
    for (const auto& di_dj : matching) {
        uint di = di_dj.first;
        uint dj = di_dj.second;
        if (visited.count(di) || visited.count(dj)) {
            continue;
        }
        qrc::DecodingGraph::Vertex * vdi = graph.get_vertex(di);
        qrc::DecodingGraph::Vertex * vdj = graph.get_vertex(dj);
        auto vdi_vdj = std::make_pair(vdi, vdj);
        std::vector<qrc::DecodingGraph::Vertex*> detector_path(path_table[vdi_vdj].path);
        std::cout << "(" << di << ", " << dj << ") -- " << 
                        path_table[vdi_vdj].distance << " -- " <<
                        detector_path.size() - 1 << " -" ;
        total_path_d += detector_path.size() - 1; 
        for (uint i = 0; i < detector_path.size(); i++) {
            // Get edge from decoding graph.
            auto wi = detector_path[i];
            std::cout <<  "- " << wi->detector;
        }
        std::cout << std::endl;

        visited.insert(di);
        visited.insert(dj);
    }
    std::cout << " Total distance: " << total_path_d << std::endl;
    return total_path_d;


}

fp_t 
Predecoder::graph_sparcity_calculator(std::vector<uint16_t> compressed_synd, fp_t& sparcity_unweighted, uint distance){
    fp_t sparcity_weighted = 0;
    sparcity_unweighted = 0;
    fp_t min_ = 0;
    fp_t min_p = 0;

    auto syndrome_ = syndrome_decompressed(compressed_synd, distance);

    update_adjacency_matrix_and_detector_list(syndrome_, 1, true);
    auto subgraph = get_adjacency_matrix();
    auto components = qpd::find_connected_components(subgraph.adjacency_matrix);

    for(uint i =0; i< components.size(); i++){
        if(components[i].size() % 2 ==1)
            return 1; 
    }
    return 0;

    for(uint i =0; i< components.size(); i++){        
        for(uint j = i+1; j< components.size(); j++){
            min_ = DBL_MAX;
            min_p = DBL_MAX;
            for(const auto& vdi: components[i]){
                for(const auto& vdj: components[j]){

                    auto vdi_vdj = std::make_pair(vdi, vdj);
                    std::vector<qrc::DecodingGraph::Vertex*> detector_path(path_table[vdi_vdj].path);
                    // sparcity_unweighted += ((fp_t)1/detector_path.size());
                    // sparcity_weighted += ((fp_t)1/path_table[vdi_vdj].distance);
                    if(min_ > ((fp_t )detector_path.size())){
                        min_ = ((fp_t )detector_path.size());
                    }
                    if(min_p > path_table[vdi_vdj].distance){
                        min_p = path_table[vdi_vdj].distance;
                    }

                }
            }
            sparcity_unweighted += ((fp_t)min_);
            sparcity_weighted += (min_p);

        }
    }
    uint n = components.size();
    if(n == 1){
        sparcity_unweighted = 0;
        sparcity_weighted = 0 ;
        return sparcity_weighted;

    }
    sparcity_unweighted = ((fp_t)sparcity_unweighted/(n*(n-1)/2));
    sparcity_weighted = ((fp_t)sparcity_weighted/(n*(n-1)/2));
    
    // fp_t sparcity_weighted = 0;
    // sparcity_unweighted = 0;
    // bool extra_node = false;
    // uint max_size = ((distance+1)*(distance*distance - 1)/2);
    // for(uint i =0; i<compressed_synd.size(); i++){
    //     if(compressed_synd[i] == max_size){
    //             extra_node = true;
    //             continue;
    //     }
    //     for (uint j = i+1; j < compressed_synd.size(); j++){
    //         if(compressed_synd[j] == max_size){
    //             extra_node = true;
    //             continue;
    //         }
    //         qrc::DecodingGraph::Vertex * vdi = graph.get_vertex(compressed_synd[i]);
    //         qrc::DecodingGraph::Vertex * vdj = graph.get_vertex(compressed_synd[j]);
    //         auto vdi_vdj = std::make_pair(vdi, vdj);
    //         std::vector<qrc::DecodingGraph::Vertex*> detector_path(path_table[vdi_vdj].path);
    //         // sparcity_unweighted += ((fp_t)1/detector_path.size());
    //         // sparcity_weighted += ((fp_t)1/path_table[vdi_vdj].distance);
    //         sparcity_unweighted += ((fp_t)detector_path.size());
    //         sparcity_weighted += (path_table[vdi_vdj].distance);

    //         // std::cout <<compressed_synd[i]<< "-"<< compressed_synd[j] << "sparcity_unweighted " << sparcity_unweighted <<
    //         //             " detector_path.size(): " <<detector_path.size() <<
    //         //             " sparcity_weighted: " << sparcity_weighted <<  
    //         //             " path_table[vdi_vdj].distance: " << path_table[vdi_vdj].distance << std::endl;
    //     }  
    // }
    // uint n = compressed_synd.size();
    // if(extra_node){
    //     n -= 1;
    // }
    // sparcity_unweighted = ((fp_t)sparcity_unweighted/(n*(n-1)/2));
    // sparcity_weighted = ((fp_t)sparcity_weighted/(n*(n-1)/2));
    // std::cout << "==========end======" << std::endl;
    return sparcity_weighted;
}



} //namesapce qpd