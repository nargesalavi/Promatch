/*
    author: Narges Alavisamani
*/

#ifndef ERROR_EVENT_h
#define ERROR_EVENT_h

#include "decoding_graph.h"

namespace bbsim{

class ErrorEvent : public qrc::DecodingGraph::Edge{
    public:
        fp_t cumulative_probability;
        
        ErrorEvent();
        ErrorEvent(int32_t id, uint di, uint dj, fp_t w, fp_t p, std::set<uint> frames);
        ErrorEvent(const ErrorEvent& other);
        ErrorEvent(int32_t id, uint di, uint dj, fp_t w, fp_t p, std::set<uint> frames, fp_t cumulative_probability);
        ~ErrorEvent();

        void set_cumulative_probability(fp_t g_probab_sum);

        bool compare_events_prob(const ErrorEvent& e1, const ErrorEvent& e2);

};


}

#endif