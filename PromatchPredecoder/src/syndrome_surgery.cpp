/*
 *  author: Narges Alavisamani
 *  date:   28 March 2023
 * */

#include "syndrome_surgery.h"

namespace qpd{

void print_graph(qrc::DecodingGraph decoding_graph){
    std::cout << "__________print_graph___________" << std::endl;
    std::vector<qrc::DecodingGraph::Vertex*> _vertices = decoding_graph.vertices();
    for(qrc::DecodingGraph::Vertex* v : _vertices){
        std::vector<qrc::DecodingGraph::Vertex*> ad_list = decoding_graph.adjacency_list(v);
        std::cout << "\n___Vertex__ " << v->detector;
        if(v->detector!=BOUNDARY_INDEX)
            std::cout<< "(" << v->coord[0]<<","<<v->coord[1]<<","<<v->coord[2]<<")";
        for (int i=3;i<N_COORD;i++){
            if(v->coord[i] && v->coord[i]!=BOUNDARY_INDEX)
                std::cout << i << ":" << v->coord[i] << "^" << std::endl;
        }
        std::cout << std::endl;
        for (qrc::DecodingGraph::Vertex* adjecnt_v: ad_list){
            std::cout << adjecnt_v->detector;
            if(adjecnt_v->detector!=BOUNDARY_INDEX)
                std::cout <<"(" <<adjecnt_v->coord[0] << "," << adjecnt_v->coord[1]<<","<<adjecnt_v->coord[2]<< ")";
            std::cout << ",";
        }
    }

    std::cout << std::endl << "__________end of print_graph___________" << std::endl;
    auto path_table = qrc::compute_path_table(decoding_graph);
    for(auto v1 : _vertices){
        for(auto v2 : _vertices){
            auto vdi_vdj = std::make_pair(v1, v2);
            if(path_table[vdi_vdj].path.size() == 2)
                std::cout << "*" << v1->detector << "," <<
                v2->detector << ",";

        }
    }
    std::cout << std::endl;
}

void print_graph_and_syndrome(qrc::DecodingGraph graph, const std::vector<uint8_t>& syndrome){

    print_graph(graph);
    std::cout << "\nSyndrome:\n";
            for (int id =0; id<syndrome.size();id++) {
                std::cout<< id<<":"<< ((int)syndrome[id] )<< ", ";
            }
    std::cout << std::endl;
}

void testing_syndroms_and_graphs(uint64_t total_shots, uint distance, 
        fp_t physcial_error, fp_t meas_er){
       
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
        1,
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
    qrc::DecodingGraph decoding_graph = qrc::to_decoding_graph(circ);

    const uint n_detectors = circ.count_detectors();

    while (shots > 0) {
        uint32_t shots_this_round = shots > SHOTS_PER_BATCH ? SHOTS_PER_BATCH : shots;
        // Make syndromes
        auto sample_buffer = stim::detector_samples(circ, shots_this_round, false, true, rng);
        sample_buffer = sample_buffer.transposed();

        for (uint64_t i = 0; i < shots_this_round; i++) {
            auto syndrome = qrc::_to_vector(sample_buffer[i], 
                                    circ.count_detectors(),
                                    circ.count_observables());
            uint hw = std::accumulate(syndrome.begin(), syndrome.begin() + circ.count_detectors(), 0);
            std::cout << "HW: " << hw << std::endl;
            print_graph_and_syndrome(decoding_graph,syndrome);

        }
        shots -= shots_this_round;
    }
}

std::vector<uint16_t> syndrome_compressed(const std::vector<uint8_t>& syndrome, bool print_syndrome){
    std::vector<uint16_t> s;

    for(uint16_t i=0; i<syndrome.size();i++ ){
        if((int)syndrome[i]){
            s.push_back(i);
            if(print_syndrome){
                std::cout << i << " ";
            }
        }
    }
    if(print_syndrome){
        std::cout << std::endl;
    }
    return s;
}

std::vector<uint8_t> syndrome_decompressed(std::vector<uint16_t> synd_vec, uint d){
    std::vector<uint8_t> syndrome;
    // std::cout << "Resizing started";
    syndrome.resize(((d+1)*(d*d - 1)/2)+1, ((uint8_t)0));
    // std::cout << " decomp start";
    for( auto indx : synd_vec){
        syndrome[indx] = ((uint8_t)1);
    }
    // std::cout << " decomp end" << std::endl;
    return syndrome;
}

void write_vector(const std::vector<uint16_t>& vec, const std::string& filename) {
    // Open the file for binary output
    // std::ofstream file(filename, std::ios::binary);

    // // Write the size of the vector as a 32-bit unsigned integer
    // uint32_t size = static_cast<uint32_t>(vec.size());
    // file.write(reinterpret_cast<const char*>(&size), sizeof(size));

    // // Write the elements of the vector
    // file.write(reinterpret_cast<const char*>(vec.data()), vec.size() * sizeof(uint16_t));

    // file.close();

    std::ofstream outfile(filename);
    for (const auto& val : vec) {
        outfile << val << " ";
    }
    outfile.close();
}

std::vector<uint16_t> read_vector(const std::string& filename) {
    // Open the file for binary input
    // std::ifstream file(filename, std::ios::binary);

    // // Read the size of the vector as a 32-bit unsigned integer
    // uint32_t size;
    // file.read(reinterpret_cast<char*>(&size), sizeof(size));

    // // Resize the vector and read the elements
    // std::vector<uint16_t> vec(size);
    // file.read(reinterpret_cast<char*>(vec.data()), size * sizeof(uint16_t));
    // file.close();

    // return vec;

    std::vector<uint16_t> vec;
    std::ifstream infile(filename);
    uint16_t val;
    while (infile >> val) {
        vec.push_back(val);
    }
    infile.close();
    return vec;
}


// Function to write a std::vector<std::vector<uint16_t>> to a binary file
void write_vector_of_vectors(const std::vector<std::vector<uint16_t>>& data, const std::string& filename) {
    // Open the file in binary mode
    // std::ofstream file(filename, std::ios::binary);
    
    // // Write the number of rows and columns to the file
    // std::uint16_t num_rows = static_cast<std::uint16_t>(data.size());
    // std::uint16_t num_cols = static_cast<std::uint16_t>(data[0].size());
    // file.write(reinterpret_cast<const char*>(&num_rows), sizeof(num_rows));
    // file.write(reinterpret_cast<const char*>(&num_cols), sizeof(num_cols));

    // // Write the data to the file
    // for (const auto& row : data) {
    //     for (const auto& val : row) {
    //         std::uint16_t value = static_cast<std::uint16_t>(val);
    //         file.write(reinterpret_cast<const char*>(&value), sizeof(value));
    //     }
    // }

    // // Close the file
    // file.close();

    std::ofstream ofs(filename);
    if (!ofs.is_open()) {
        std::cout << "Error opening file " << filename << " for writing" << std::endl;
        return;
    }

    for (const auto& vec : data) {
        for (const auto& elem : vec) {
            ofs << elem << " ";
        }
        ofs << std::endl;
    }
}

// Function to read a std::vector<std::vector<uint16_t>> from a binary file
std::vector<std::vector<uint16_t>> read_vector_of_vectors(const std::string& filename) {
    // Open the file in binary mode
    // std::ifstream file(filename, std::ios::binary);

    // // Read the number of rows and columns from the file
    // std::uint16_t num_rows, num_cols;
    // file.read(reinterpret_cast<char*>(&num_rows), sizeof(num_rows));
    // file.read(reinterpret_cast<char*>(&num_cols), sizeof(num_cols));

    // // Read the data from the file
    // std::vector<std::vector<uint16_t>> data(num_rows, std::vector<uint16_t>(num_cols));
    // for (auto& row : data) {
    //     for (auto& val : row) {
    //         std::uint16_t value;
    //         file.read(reinterpret_cast<char*>(&value), sizeof(value));
    //         val = static_cast<uint16_t>(value);
    //     }
    // }

    std::ifstream ifs(filename);
    if (!ifs.is_open()) {
        std::cout << "Error opening file " << filename << " for reading" << std::endl;
        return {};
    }

    std::vector<std::vector<uint16_t>> vec_of_vecs;
    std::string line;
    while (std::getline(ifs, line)) {
        std::vector<uint16_t> vec;
        std::istringstream iss(line);
        uint16_t elem;
        while (iss >> elem) {
            vec.push_back(elem);
        }
        vec_of_vecs.push_back(vec);
    }

    return vec_of_vecs;

    // // Close the file
    // file.close();

    // return data;
}

void print_syndrome(const std::vector<uint8_t>& syndrome){
    std::cout << "Syndrome:" << std::endl;
    for (int id =0; id<syndrome.size();id++) {
        std::cout<< id<<":"<< ((int)syndrome[id] )<< ", ";
    }

    std::cout << std::endl;
}


}//namespace qpd