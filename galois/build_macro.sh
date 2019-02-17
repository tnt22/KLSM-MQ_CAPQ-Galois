#!/bin/bash

# Those lines make the script stop in case of failure.
set -e
set -o pipefail # needed in case we use a pipe in the script.

# Build it (with debug flag).
cd build
mkdir debug
cd debug
cmake -DCMAKE_BUILD_TYPE=Debug ../..
make install

# if completed succesfully, can go on now.