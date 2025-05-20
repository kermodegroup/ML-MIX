/* ----------------------------------------------------------------------
   This file is part of a plugin for LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2025) Fraser Birks
   This file was written by Fraser Birks and is distributed under the 
   GNU General Public License.

   This file is not part of the original LAMMPS distribution but is designed 
   to be used with LAMMPS. It is provided under the same GPLv2 license as LAMMPS 
   to ensure compatibility.

   See the LICENSE file for details.
------------------------------------------------------------------------- */

#include "fix_mlml_kokkos.h"
#include "atom_kokkos.h"

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
template<class DeviceType>
FixMLMLKokkos::FixMLMLKokkos(LAMMPS *lmp, int narg, char **arg) : FixMLML(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK //only reads atom positions from device, does not modify anything
}

FixMLMLKokkos::~FixMLMLKokkos(){}

void FixMLMLKokkos::allocate_regions(double**){
  // kokkos version of allocate regions which handles device memory
  atomKK->sync(LMPHostType,datamask_read);
  x = atomKK->k_x.view<LMPHostType>();
  FixMLML::allocate_regions(x);
}

namespace LAMMPS_NS {
template class FixLearnKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixLearnKokkos<LMPHostType>;
#endif
}