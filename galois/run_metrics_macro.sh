#!/bin/bash

# Those lines make the script stop in case of failure.
set -e
set -o pipefail # needed in case we use a pipe in the script.

#export LD_LIBRARY_PATH=/usr/lib/libblas:/usr/lib/atlas-base:${LD_LIBRARY_PATH}
#export LIBRARY_PATH=/usr/lib/libblas:/usr/lib/atlas-base:${LIBRARY_PATH}

# Compile the files
chmod 777 build_macro.sh
./build_macro.sh

# Constants
GRAPH_LOCATION="/specific/disk1/home/mad/Galois-PQ/inputs"
runs_per_benchmark=1
out_dir="tempTest3"
echo $out_dir
standard_parameters="-m jemalloc -D 4 -d $out_dir sssp"
echo $standard_parameters


for threads_number in 176 88 44 22 11 5 3 1
do
    for g in scalefree/rmat16p-2e24.gr road/USA-road-t.USA.gr clique4000.gr scalefree/rmat16p-2e24.gr scalefree/rmat16p-2e27.gr road/USA-road-d.USA.gr road/USA-road-t.USA.gr ljournal-2008.gr twitter40.gr random/r4-2e24.gr # planar10M.bin scalefree/rmat-large.gr road/USA-road-d.USA.gr road/USA-road-t.USA.gr # clique4000.bin
    do
        for worklist_type in capq klsm256_capq capq klsm256 klsm256_capq klsm2 klsm4 klsm8 klsm16 klsm32 klsm64 klsm128 klsm256 klsm512 klsm1024 klsm2048 klsm4096
        do
            echo -e "\e[33mCommand: python run.py -t ${threads_number} -r $runs_per_benchmark -g "$GRAPH_LOCATION/${g}" -v ${worklist_type} $standard_parameters\e[0m"
            python run.py -t ${threads_number} -r $runs_per_benchmark -g "$GRAPH_LOCATION/${g}" -v ${worklist_type} $standard_parameters
        done
    done
done

