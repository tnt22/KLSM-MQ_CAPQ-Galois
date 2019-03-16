#!/bin/bash

# Those lines make the script stop in case of failure.
set -e
set -o pipefail # needed in case we use a pipe in the script.

# parameter for build.
export CACHE_LINE_SIZE="128" # This is used for the CAPQ code compilation.

# Build it (with debug flag).
cd build
mkdir -p debug
cd debug
cmake -DCMAKE_C_COMPILER=gcc-5.3.0 -DCMAKE_BUILD_TYPE=Debug ../..
make sssp -j # the -j makes it parallel

# if completed succesfully, can go on now.
echo -e "Build finished \e[32mSuccessfully\e[0m"

