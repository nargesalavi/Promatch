cmake_minimum_required(VERSION 3.0.0)
project(dramsim3)

set(default_build_type "Release")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

add_library(inih INTERFACE)
target_include_directories(inih INTERFACE ext/headers)

add_library(format INTERFACE)
target_include_directories(format INTERFACE ext/fmt/include)
target_compile_definitions(format INTERFACE FMT_HEADER_ONLY=1)

# argparsing library, only used in main program not the library
add_library(args INTERFACE)
target_include_directories(args INTERFACE ext/headers)

add_library(json INTERFACE)
target_include_directories(json INTERFACE ext/headers)

# Main DRAMSim Lib
add_library(dramsim3 SHARED
    src/bankstate.cc
    src/channel_state.cc
    src/command_queue.cc
    src/common.cc
    src/configuration.cc
    src/controller.cc
    src/dram_system.cc
    src/hmc.cc
    src/refresh.cc
    src/simple_stats.cc
    src/timing.cc
    src/memory_system.cc
    src/rhtracker.cc
)

if (THERMAL)
    # dependency check
    if (APPLE) 
        find_package(BLAS REQUIRED
            HINTS /usr/local/atlas/lib)
    else()
        # sudo apt-get install libatlas-base-dev on ubuntu
        find_package(BLAS REQUIRED)
    endif()
    find_package(OpenMP REQUIRED)
    # YOU need to build superlu on your own. Do the following:
    # git submodule update --init
    # cd ext/SuperLU_MT_3.1 && make lib
    find_library(SUPERLU
        NAME superlu_mt_OPENMP libsuperlu_mt_OPENMP
        HINTS ${PROJECT_SOURCE_DIR}/ext/SuperLU_MT_3.1/lib/
    )

    target_link_libraries(dramsim3
        PRIVATE ${SUPERLU} f77blas atlas m ${OpenMP_C_FLAGS}
    )
    target_sources(dramsim3
        PRIVATE src/thermal.cc src/sp_ienv.c src/thermal_solver.c
    )
    target_compile_options(dramsim3 PRIVATE -DTHERMAL -D_LONGINT -DAdd_ ${OpenMP_C_FLAGS})

    add_executable(thermalreplay src/thermal_replay.cc)
    target_link_libraries(thermalreplay dramsim3 inih)
    target_compile_options(thermalreplay PRIVATE -DTHERMAL -D_LONGINT -DAdd_ ${OpenMP_C_FLAGS})
endif (THERMAL)

if (ROWHAMMER)
    target_compile_definitions(dramsim3 PRIVATE ROWHAMMER)
endif()

if (CMD_TRACE)
    target_compile_options(dramsim3 PRIVATE -DCMD_TRACE)
endif (CMD_TRACE)

if (ADDR_TRACE)
    target_compile_options(dramsim3 PRIVATE -DADDR_TRACE)
endif (ADDR_TRACE)


target_include_directories(dramsim3 INTERFACE src)
target_compile_options(dramsim3 PRIVATE -Wall)
target_link_libraries(dramsim3 PUBLIC inih format)
set_target_properties(dramsim3 PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)
