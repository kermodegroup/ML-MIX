#!/bin/bash

# Example build script to build LAMMPS on 20 processors
# Assumes that the LAMMPS source code is in a directory called 'lammps'
# Assumes that LAMMPS is being built into a virtual environment

cd lammps

mkdir build

cd build

wget -O libpace.tar.gz https://github.com/wcwitt/lammps-user-pace/archive/main.tar.gz

cmake ../cmake  -D CMAKE_BUILD_TYPE=Release \
                -D BUILD_SHARED_LIBS=yes \
                -D CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV \
                -D PKG_PYTHON=yes \
                -D PKG_ML-PACE=yes \
                -D PKG_OPENMP=yes \
                -D BUILD_MPI=yes \ 
                -D PACELIB_MD5=$(md5sum libpace.tar.gz | awk '{print $1}') \
                -D PKG_RIGID=yes \
                -D PKG_ML-UF3=yes \
                -D PKG_PLUGIN=yes \
                -D PKG_EXTRA-PAIR=yes \
                -D PKG_MOLECULE=yes \
                -D PKG_MANYBODY=yes \
                -D PKG_REPLICA=yes \
               ../cmake

cmake --build . -j 20

cmake --install .

make install-python


#get current directory
current_dir=$(pwd)

# add to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$current_dir