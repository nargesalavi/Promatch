/*
 *  author: Suhas Vittal
 *  date:   8 September 2022
 * */

#ifndef __RHTRACKER_H
#define __RHTRACKER_H

#include "common.h"

#include <map>

namespace dramsim3_ext {

/* Everything here is not intended for a hardware
 * implementation. This is just tracking data structures
 * for other uses. */

#define T_RH 2048
    
struct RHTrackerEntry {
    RHTrackerEntry();

    uint32_t n_rh_flips;
    uint64_t n_row_acts;    // Resets on refresh.
    
    bool increment(void);  // Returns true if the RH threshold is reached.
};

class RHTracker {
public:
    RHTracker();

    void register_activation(const dramsim3::Address&);
    void refresh_rank(int);
    void refresh_bank(int bank, int rank);

    uint64_t n_rh_flips;
    uint64_t n_activations;
private:
    std::map<dramsim3::Address, RHTrackerEntry> rh_table;
};


}   // dramsim3_ext

#endif
