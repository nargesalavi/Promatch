/*
 *  author: Suhas Vittal
 *  date:   3 February 2023
 * */

#ifndef BENCHMARK_STATBENCH_h
#define BENCHMARK_STATBENCH_h

#include "defs.h"

#include <functional>
#include <map>

namespace qrc {

class Decoder;
/*
 *  Decoder pointers from a dgf_t should be
 *  allocated on the heap.
 * */
typedef std::function<Decoder*(fp_t)> dgf_t;

namespace benchmark {

struct StatisticalResult {
    fp_t logical_error_rate = 0;
    uint64_t n_logical_errors = 0;
    fp_t mean_execution_time = 0;
    fp_t max_execution_time = 0;
};

}   // benchmark
}   // qrc

#endif // BENCHMARK_STATBENCH
