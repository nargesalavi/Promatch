cmake_minimum_required(VERSION 3.20.2)
project(Promatch VERSION 0.1)
# Declare CPP Version
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_subdirectory(quarch)

###################
#### PromatchPredecoder ####
###################

set(PROMATCHPREDECODER_FILES
    PromatchPredecoder/src/main.cpp
    PromatchPredecoder/src/predecoder.cpp
    PromatchPredecoder/src/experiments.cpp
    PromatchPredecoder/src/predecoding_graph_helper.cpp
    PromatchPredecoder/src/syndrome_surgery.cpp
    PromatchPredecoder/src/bucket_based_simulation.cpp
    PromatchPredecoder/src/error_event.cpp
    PromatchPredecoder/src/bb_benchmark.cpp)
add_executable(promatch ${PROMATCHPREDECODER_FILES})
if (CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(promatch PUBLIC -O3)
else()
    target_compile_options(promatch PUBLIC -ggdb3)
endif()

target_link_libraries(promatch PRIVATE quarch)
target_include_directories(promatch PRIVATE "PromatchPredecoder/include")
find_package(MPI REQUIRED)
target_compile_definitions(promatch PRIVATE HOME_DIRECTORY="${CMAKE_SOURCE_DIR}")
target_include_directories(promatch PRIVATE ${MPI_C_INCLUDE_PATH})
target_include_directories(promatch PRIVATE ${MPI_CXX_INCLUDE_PATH})
#target_link_libraries(promatch PRIVATE MPI::MPI_CXX)
target_link_libraries(promatch PRIVATE ${MPI_C_LIBRARIES})
target_link_libraries(promatch PRIVATE ${MPI_CXX_LIBRARIES})

# Set number of processes as an argument for mpirun
# if(MPI_COMPILE_FLAGS)
#   set_target_properties(promatch PROPERTIES
#     COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
# endif()

# if(MPI_LINK_FLAGS)
#   set_target_properties(promatch PROPERTIES
#     LINK_FLAGS "${MPI_LINK_FLAGS}")
# endif()
# ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} PROCS
#   ${MPIEXEC_PREFLAGS} EXECUTABLE ${MPIEXEC_POSTFLAGS} ARGS

# add_custom_target(run
#   COMMAND mpirun -np ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/promatch ${MPIEXEC_POSTFLAGS} ${MPIEXEC_NUMPROC_FLAG}
#   DEPENDS nfd)