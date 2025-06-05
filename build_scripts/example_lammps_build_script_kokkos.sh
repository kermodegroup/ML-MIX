#!/bin/bash

# Example build script to build LAMMPS on 20 processors
# Assumes that the LAMMPS source code is in a directory called 'lammps'
# Assumes that LAMMPS is being built into a virtual environment

# !!!!!!! IMPORTANT !!!!!!! #
# Make sure that you have run the ./install.sh script in the LAMMPS plugin directory
# before running this script, otherwise the ML-MIX package will not install.

# The following script is customised for building LAMMPS with A100 GPUs
# If you are using a different GPU architecture, you need to modify -D Kokkos_ARCH_AMPERE86
# accordingly. Follow the instructions in the LAMMPS documentation: https://docs.lammps.org/Build_extras.html#kokkos

cd lammps

mkdir build

cd build

wget -O libpace.tar.gz https://github.com/wcwitt/lammps-user-pace/archive/main.tar.gz

cmake ../cmake \
  -D CMAKE_CXX_COMPILER=/home/eng/phrffk/lammps-kokkos/lib/kokkos/bin/nvcc_wrapper \
  -D CMAKE_BUILD_TYPE=Release \
  -D BUILD_MPI=on \
  -D PKG_KOKKOS=yes \
  -D Kokkos_ENABLE_SERIAL=ON \
  -D Kokkos_ARCH_AMPERE86=yes \
  -D BUILD_SHARED_LIBS=yes \
  -D BUILD_OMP=ON \
  -D Kokkos_ENABLE_CUDA=yes \
  -D CMAKE_INSTALL_PREFIX=$VIRTUAL_ENV \
  -D PACELIB_MD5=$(md5sum libpace.tar.gz | awk '{print $1}') \
  -D PKG_ML-UF3=yes \
  -D PKG_ML-PACE=yes \
  -D PKG_ML-SNAP=yes \
  -D PKG_RIGID=yes \
  -D PKG_MANYBODY=yes \
  -D PKG_MOLECULE=yes \
  -D PKG_EXTRA-PAIR=yes \


cmake --build . -j 20

cmake --install .

make install-python


#get current directory
current_dir=$(pwd)

# add to LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$current_dir
