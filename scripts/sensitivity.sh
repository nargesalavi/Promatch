#!/bin/bash

# Navigate to the PromatchDecoder directory (adjust the path as necessary)
cd ..

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
create_dir_if_not_exists "AEResults/PER_d11/Predecoders/2x"
create_dir_if_not_exists "AEResults/PER_d13/Predecoders/2x"
create_dir_if_not_exists "AEResults/PER_d11/MWPM/2x"
create_dir_if_not_exists "AEResults/PER_d13/MWPM/2x"

create_dir_if_not_exists "AEResults/PER_d11/Predecoders/3x"
create_dir_if_not_exists "AEResults/PER_d13/Predecoders/3x"
create_dir_if_not_exists "AEResults/PER_d11/MWPM/3x"
create_dir_if_not_exists "AEResults/PER_d13/MWPM/3x"

create_dir_if_not_exists "AEResults/PER_d11/Predecoders/4x"
create_dir_if_not_exists "AEResults/PER_d13/Predecoders/4x"
create_dir_if_not_exists "AEResults/PER_d11/MWPM/4x"
create_dir_if_not_exists "AEResults/PER_d13/MWPM/4x"

create_dir_if_not_exists "AEResults/PER_d11/Predecoders/5x"
create_dir_if_not_exists "AEResults/PER_d13/Predecoders/5x"
create_dir_if_not_exists "AEResults/PER_d11/MWPM/5x"
create_dir_if_not_exists "AEResults/PER_d13/MWPM/5x"


# Check if the number of cores is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <number_of_cores>"
    exit 1
fi

NUM_CORES=$1

# Distance 11 for MWPM p = 2x
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 2 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/2x/mwpm_d11_r0_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 2 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/2x/mwpm_d11_r1_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 2 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/2x/mwpm_d11_r2_2xp.txt
# Distance 11 for all predecoders p = 2x
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 2 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/2x/predecoders_d11_r0_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 2 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/2x/predecoders_d11_r1_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 2 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/2x/predecoders_d11_r2_2xp.txt
# Distance 11 for MWPM p = 3x
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 3 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/3x/mwpm_d11_r0_3xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 3 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/3x/mwpm_d11_r1_3xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 3 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/3x/mwpm_d11_r2_3xp.txt
# Distance 11 for all predecoders p = 3x
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 3 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/3x/predecoders_d11_r0_3xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 3 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/3x/predecoders_d11_r1_3xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 3 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/3x/predecoders_d11_r2_3xp.txt
# Distance 11 for MWPM p = 4x
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 4 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/4x/mwpm_d11_r0_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 4 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/4x/mwpm_d11_r1_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 4 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/4x/mwpm_d11_r2_4xp.txt
# Distance 11 for all predecoders p = 4x
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 4 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/4x/predecoders_d11_r0_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 4 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/4x/predecoders_d11_r1_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 4 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/4x/predecoders_d11_r2_4xp.txt
# Distance 11 for MWPM p = 5x
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 5 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/5x/mwpm_d11_r0_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 5 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/5x/mwpm_d11_r1_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 5 -t MWPM -ib 1 | tee ../AEResults/PER_d11/MWPM/5x/mwpm_d11_r2_5xp.txt
# Distance 11 for all predecoders p = 5x
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 0 -ms 100M -hs 100M -ls 100M -p 5 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/5x/predecoders_d11_r0_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 1 -ms 100M -hs 100M -ls 100M -p 5 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/5x/predecoders_d11_r1_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 11 -r 2 -ms 100M -hs 100M -ls 100M -p 5 -g 0 -ib 1 | tee ../AEResults/PER_d11/Predecoders/5x/predecoders_d11_r2_5xp.txt


# Distance 13 for MWPM p = 2x
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 0 -ms 100M -hs 100M -ls 100M -p 2 -t MWPM | tee ../AEResults/PER_d13/MWPM/2x/mwpm_d13_r0_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 1 -ms 100M -hs 100M -ls 100M -p 2 -t MWPM | tee ../AEResults/PER_d13/MWPM/2x/mwpm_d13_r1_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 2 -ms 100M -hs 100M -ls 100M -p 2 -t MWPM | tee ../AEResults/PER_d13/MWPM/2x/mwpm_d13_r2_2p.txt
# Distance 13 for all predecoders p = 2x
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 0 -ms 200M -hs 20M -ls 1M -p 2 -g 0 | tee ../AEResults/PER_d13/Predecoders/2x/predecoders_d13_r0_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 1 -ms 200M -hs 20M -ls 1M -p 2 -g 0 | tee ../AEResults/PER_d13/Predecoders/2x/predecoders_d13_r1_2xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 2 -ms 200M -hs 20M -ls 1M -p 2 -g 0 | tee ../AEResults/PER_d13/Predecoders/2x/predecoders_d13_r2_2xp.txt
# Distance 13 for MWPM p = 3x
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 0 -ms 100M -hs 100M -ls 100M -p 3 -t MWPM | tee ../AEResults/PER_d13/MWPM/1x/mwpm_d13_r0_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 1 -ms 100M -hs 100M -ls 100M -p 3 -t MWPM | tee ../AEResults/PER_d13/MWPM/1x/mwpm_d13_r1_1xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 2 -ms 100M -hs 100M -ls 100M -p 3 -t MWPM | tee ../AEResults/PER_d13/MWPM/1x/mwpm_d13_r2_1xp.txt
# Distance 13 for all predecoders p = 3x
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 0 -ms 200M -hs 20M -ls 1M -p 3 -g 0 | tee ../AEResults/PER_d13/Predecoders/3x/predecoders_d13_r0_3xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 1 -ms 200M -hs 20M -ls 1M -p 3 -g 0 | tee ../AEResults/PER_d13/Predecoders/3x/predecoders_d13_r1_3xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 2 -ms 200M -hs 20M -ls 1M -p 3 -g 0 | tee ../AEResults/PER_d13/Predecoders/3x/predecoders_d13_r2_3xp.txt
# Distance 13 for MWPM p = 4x
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 0 -ms 100M -hs 100M -ls 100M -p 4 -t MWPM | tee ../AEResults/PER_d13/MWPM/4x/mwpm_d13_r0_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 1 -ms 100M -hs 100M -ls 100M -p 4 -t MWPM | tee ../AEResults/PER_d13/MWPM/4x/mwpm_d13_r1_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 2 -ms 100M -hs 100M -ls 100M -p 4 -t MWPM | tee ../AEResults/PER_d13/MWPM/4x/mwpm_d13_r2_4xp.txt
# Distance 13 for all predecoders p = 4x
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 0 -ms 200M -hs 20M -ls 1M -p 4 -g 0 | tee ../AEResults/PER_d13/Predecoders/4x/predecoders_d13_r0_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 1 -ms 200M -hs 20M -ls 1M -p 4 -g 0 | tee ../AEResults/PER_d13/Predecoders/4x/predecoders_d13_r1_4xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 2 -ms 200M -hs 20M -ls 1M -p 4 -g 0 | tee ../AEResults/PER_d13/Predecoders/4x/predecoders_d13_r2_4xp.txt
# Distance 13 for MWPM p = 3x
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 0 -ms 100M -hs 100M -ls 100M -p 5 -t MWPM | tee ../AEResults/PER_d13/MWPM/5x/mwpm_d13_r0_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 1 -ms 100M -hs 100M -ls 100M -p 5 -t MWPM | tee ../AEResults/PER_d13/MWPM/5x/mwpm_d13_r1_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 40 -d 13 -r 2 -ms 100M -hs 100M -ls 100M -p 5 -t MWPM | tee ../AEResults/PER_d13/MWPM/5x/mwpm_d13_r2_5xp.txt
# Distance 13 for all predecoders p = 3x
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 0 -ms 200M -hs 20M -ls 1M -p 5 -g 0 | tee ../AEResults/PER_d13/Predecoders/5x/predecoders_d13_r0_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 1 -ms 200M -hs 20M -ls 1M -p 5 -g 0 | tee ../AEResults/PER_d13/Predecoders/5x/predecoders_d13_r1_5xp.txt
mpirun -np $NUM_CORES ./promatch -e 39 -d 13 -r 2 -ms 200M -hs 20M -ls 1M -p 5 -g 0 | tee ../AEResults/PER_d13/Predecoders/5x/predecoders_d13_r2_5xp.txt

