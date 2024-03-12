/*
 *  author: Suhas Vittal
 *  date:   6 August 2022
 * */

#include "defs.h"

namespace qrc {

qfp_t 
quantize(fp_t orig, fp_t fp_max, qfp_t qfp_max) {
    fp_t mid = (orig/fp_max) * qfp_max;
    return (qfp_t)mid;
}

void
safe_create_directory(const std::filesystem::path& path) {
    if (!std::filesystem::exists(path)) {
        std::filesystem::create_directory(path);
    }
}

}  // qrc
