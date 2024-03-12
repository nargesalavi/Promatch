/*
 *  author: Suhas Vittal
 *  date:   2 August 2022
 * */

#ifndef DEFS_h
#define DEFS_h

#include <limits>
#include <filesystem>

#include <stdint.h>

typedef double fp_t;   // floating point type
typedef int32_t qfp_t; // quantized floating point type

#ifdef __APPLE__
typedef uint32_t uint;
#endif
typedef int32_t sint;

typedef uint64_t addr_t;

#define KB 1024.0
#define MB (KB*1024.0)
#define GB (MB*1024.0)

namespace qrc {

qfp_t
quantize(fp_t, fp_t fp_max, qfp_t qfp_max);
void 
safe_create_directory(const std::filesystem::path&);

}  // qrc

#endif
