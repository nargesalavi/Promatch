cmake_minimum_required(VERSION 3.23.2)
project(DramSim2 VERSION 1.0)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(DRAMSIM2_PROGRAM_FILES
    AddressMapping.cpp
    Bank.cpp
    BankState.cpp
    BusPacket.cpp
    ClockDomain.cpp
    CommandQueue.cpp
    IniReader.cpp
    MemoryController.cpp
    MemorySystem.cpp
    MultiChannelMemorySystem.cpp
    PrintMacros.cpp
    Rank.cpp
    SimulatorObject.cpp
    Transaction.cpp)
add_library(libdramsim2 ${DRAMSIM2_PROGRAM_FILES})
set_target_properties(libdramsim2 PROPERTIES PREFIX "")
target_include_directories(libdramsim2 PUBLIC .)
target_compile_options(libdramsim2 PRIVATE -O3 -DNO_STORAGE -Wall -DDEBUG_BUILD)
