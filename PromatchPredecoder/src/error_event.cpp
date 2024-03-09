#include "error_event.h"


namespace bbsim{

ErrorEvent::ErrorEvent(): qrc::DecodingGraph::Edge(){}
ErrorEvent::ErrorEvent(int32_t id, uint di, uint dj, fp_t w, fp_t p, std::set<uint> frames): 
                    qrc::DecodingGraph::Edge(id, di, dj, w, p, frames){}
ErrorEvent::ErrorEvent(int32_t id, uint di, uint dj, fp_t w, fp_t p, std::set<uint> frames, fp_t cumulative_probability): 
                    qrc::DecodingGraph::Edge(id, di, dj, w, p, frames), cumulative_probability(cumulative_probability){}

ErrorEvent::ErrorEvent(const ErrorEvent& other){
    id = other.id;
    detectors = other.detectors;
    edge_weight = other.edge_weight; 
    error_probability = other.error_probability;
    frames = other.frames;
    cumulative_probability = other.cumulative_probability;
}
ErrorEvent::~ErrorEvent(){}

void ErrorEvent::set_cumulative_probability(fp_t g_probab_sum){
    //cumulative_probability = ((fp_t)error_probability/g_probab_sum);
    cumulative_probability = g_probab_sum;
}

bool compare_events_prob(const ErrorEvent& e1, const ErrorEvent& e2){
    return e1.error_probability < e2.error_probability;
}


}