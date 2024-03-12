# Declare CMAKE version
cmake_minimum_required(VERSION 3.20.2)
#2a4e7ca86be5b659c7227ed05c2faaf0c287e807 Declare Project and Version number
project(Quarch VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Declare executable.
set(QUARCH_PROGRAM_FILES
    # Main code.
    src/benchmark.cpp
    src/decoding_graph.cpp
    src/decoder.cpp
    src/mwpm_decoder.cpp
    src/defs.cpp
    # ASTREA
    src/astrea.cpp
    src/astrea/simulator.cpp
    src/astrea/mld_decoder.cpp
    # TMR
    src/tmr_decoder.cpp
    # FLEECE
    src/fleece.cpp
    src/fleece/lattice_graph.cpp
    src/fleece/rtanalysis.cpp
    src/fleece/hldecoder.cpp)
add_library(quarch ${QUARCH_PROGRAM_FILES})
# Add compile options
if(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(quarch PUBLIC -O3)
else()
    target_compile_options(quarch PUBLIC -ggdb3)
endif()

find_package(MPI REQUIRED)

# Add include directory to include path.
target_include_directories(quarch PUBLIC "include" ${MPI_INCLUDE_PATH})
# Link OpenMP to Quarch
# Add Stim to Quarch.
add_subdirectory(stim)
target_link_libraries(quarch libstim)
# Add Blossom V to Quarch
add_subdirectory(blossom5)
target_link_libraries(quarch libblossom5)
# Add DRAMSim to Quarch
#add_subdirectory(dramsim3)
#target_link_libraries(quarch dramsim3)

target_link_libraries(quarch ${MPI_CXX_LIBRARIES})
