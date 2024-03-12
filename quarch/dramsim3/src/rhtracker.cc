/*
 *  author: Suhas Vittal
 *  date:   8 September 2022
 * */

#include "rhtracker.h"

namespace dramsim3_ext {

RHTrackerEntry::RHTrackerEntry()
    :n_rh_flips(0),
    n_row_acts(0)
{}

bool
RHTrackerEntry::increment() {
    n_row_acts++;
    if (n_row_acts % T_RH == 0) {
        n_rh_flips++;
        return true;
    } else {
        return false;
    }
}

RHTracker::RHTracker()
    :n_rh_flips(0),
    n_activations(0),
    rh_table()
{}

void
RHTracker::register_activation(const dramsim3::Address& address) {
    dramsim3::Address row_addr(address);
    row_addr.column = 0;
    if (rh_table.count(row_addr) == 0) {
        rh_table[row_addr] = RHTrackerEntry();
    }
    
    if (rh_table[row_addr].increment()) {
        n_rh_flips++;
    }
    n_activations++;
}

void
RHTracker::refresh_rank(int rank) {
    for (auto kv_entry : rh_table) {
        dramsim3::Address address = kv_entry.first;
        RHTrackerEntry& e = kv_entry.second;

        if (address.rank == rank) {
            e.n_row_acts = 0;
        }
    }
}

void
RHTracker::refresh_bank(int bank, int rank) {
    for (auto kv_entry : rh_table) {
        dramsim3::Address address = kv_entry.first;
        RHTrackerEntry& e = kv_entry.second;

        if (address.rank == rank && address.bank == bank) {
            e.n_row_acts = 0;
        }
    }
}

} // dramsim3_ext
