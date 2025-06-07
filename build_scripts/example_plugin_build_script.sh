#!/bin/bash

#This script assumes that it is in the same directory as the ML-MIX .cpp/.h files

mkdir build
cd build
cmake ../cmake -D LAMMPS_SOURCE_DIR=/path/to/lammps/src \

cmake --build . -j 1

