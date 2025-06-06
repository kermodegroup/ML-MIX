#!/bin/bash

# Adapted from https://github.com/mir-group/flare/blob/master/lammps_plugins/install.sh

set -e

if [ "$#" -ne 1 ]; then
    echo "Usage:    ./install.sh path/to/lammps"
    exit 1
fi

lammps=$1

# add new lammps source files
ln -sf $(pwd)/src/fix_mlml.h ${lammps}/src/fix_mlml.h
ln -sf $(pwd)/src/fix_mlml.cpp ${lammps}/src/fix_mlml.cpp
ln -sf $(pwd)/src/KOKKOS/fix_mlml_kokkos.h ${lammps}/src/KOKKOS/fix_mlml_kokkos.h
ln -sf $(pwd)/src/KOKKOS/fix_mlml_kokkos.cpp ${lammps}/src/KOKKOS/fix_mlml_kokkos.cpp


ln -sf $(pwd)/src/pair_hybrid_overlay_mlml.h ${lammps}/src/pair_hybrid_overlay_mlml.h
ln -sf $(pwd)/src/pair_hybrid_overlay_mlml.cpp ${lammps}/src/pair_hybrid_overlay_mlml.cpp
ln -sf $(pwd)/src/KOKKOS/pair_hybrid_overlay_mlml_kokkos.h ${lammps}/src/KOKKOS/pair_hybrid_overlay_mlml_kokkos.h
ln -sf $(pwd)/src/KOKKOS/pair_hybrid_overlay_mlml_kokkos.cpp ${lammps}/src/KOKKOS/pair_hybrid_overlay_mlml_kokkos.cpp
