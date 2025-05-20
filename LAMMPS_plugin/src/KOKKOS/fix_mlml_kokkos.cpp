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
#include "atom_masks.h"
#include "kokkos_type.h"
#include "neigh_list_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
template<class DeviceType>
FixMLMLKokkos<DeviceType>::FixMLMLKokkos(LAMMPS *lmp, int narg, char **arg) : FixMLML(lmp, narg, arg)
{
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK; //only reads atom positions from device, does not modify anything
}

template<class DeviceType>
FixMLMLKokkos<DeviceType>::~FixMLMLKokkos(){}

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::init()
{
  double max_cutoff = fmax(rqm, fmax(bw, rblend));
  auto request = neighbor->add_request(this, NeighConst::REQ_FULL);
  request->set_cutoff(max_cutoff);
  request->set_kokkos_host(true);
  request->set_kokkos_device(false);
}

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::allocate_regions(){
  // kokkos version of allocate regions which handles device memory
  atomKK->sync(Host,datamask_read);
  x = atomKK->k_x.view<LMPHostType>();
  FixMLML::allocate_regions();
}

template<class DeviceType>
double FixMLMLKokkos<DeviceType>::get_x(int i, int j){
  return x(i, j);
}

namespace LAMMPS_NS {
template class FixMLMLKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMLMLKokkos<LMPHostType>;
#endif
}
