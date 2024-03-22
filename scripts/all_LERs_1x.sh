#!/bin/bash

# Navigate to the PromatchDecoder directory (adjust the path as necessary)

# Function to create directory if it doesn't exist
create_dir_if_not_exists() {
    if [ ! -d "$1" ]; then
        echo "Creating $1 directory..."
        mkdir -p "$1"
    else
        echo "$1 directory already exists."
    fi
}

# Check and create necessary directories
create_dir_if_not_exists "../AEResults/PER_d11/Predecoders/1x"
create_dir_if_not_exists "../AEResults/PER_d13/Predecoders/1x"
create_dir_if_not_exists "../AEResults/PER_d11/MWPM/1x"
create_dir_if_not_exists "../AEResults/PER_d13/MWPM/1x"

# Check if the number of cores is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <number_of_cores>"
    exit 1
fi

NUM_CORES=$1

# Distance 11 for MWPM
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 1 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/1x/mwpm_d11_r0_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 1 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/1x/mwpm_d11_r1_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 1 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/1x/mwpm_d11_r2_1xp.txt
# Distance 11 for all predecoders
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 1 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/1x/predecoders_d11_r0_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 1 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/1x/predecoders_d11_r1_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 1 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/1x/predecoders_d11_r2_1xp.txt

# Distance 13 for MWPM
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 0 -ms 200M -hs 20M -ls 1M -p 1 -t MWPM | tee ../AEResults/PER_d13/MWPM/1x/mwpm_d13_r0_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 1 -ms 200M -hs 20M -ls 1M -p 1 -t MWPM | tee ../AEResults/PER_d13/MWPM/1x/mwpm_d13_r1_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 2 -ms 200M -hs 20M -ls 1M -p 1 -t MWPM | tee ../AEResults/PER_d13/MWPM/1x/mwpm_d13_r2_1xp.txt
# Distance 13 for all predecoders
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 0 -ms 200M -hs 20M -ls 1M -p 1 -g 0 | tee ../AEResults/PER_d13/Predecoders/1x/predecoders_d13_r0_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 1 -ms 200M -hs 20M -ls 1M -p 1 -g 0 | tee ../AEResults/PER_d13/Predecoders/1x/predecoders_d13_r1_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 2 -ms 200M -hs 20M -ls 1M -p 1 -g 0 | tee ../AEResults/PER_d13/Predecoders/1x/predecoders_d13_r2_1xp.txt
