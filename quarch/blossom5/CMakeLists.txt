cmake_minimum_required(VERSION 3.20.2)
# Define project
project(BlossomV VERSION 2.05)
# Define C++ version.
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
# Define program files
set(BLOSSOMV_PROGRAM_FILES
    PMduals.cpp
    PMexpand.cpp
    PMinit.cpp
    PMinterface.cpp
    PMmain.cpp
    PMrepair.cpp
    PMshrink.cpp
    MinCost/MinCost.cpp
    GEOM/GPMinit.cpp
    GEOM/GPMinterface.cpp
    GEOM/GPMkdtree.cpp)
# Declared library.
add_library(libblossom5 ${BLOSSOMV_PROGRAM_FILES})
set_target_properties(libblossom5 PROPERTIES PREFIX "")
target_include_directories(libblossom5 PUBLIC .)
target_compile_options(libblossom5 PRIVATE -O3 -D_NDEBUG)
if(NOT APPLE)
    target_link_libraries(libblossom5 PUBLIC rt)
endif()

