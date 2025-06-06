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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(hybrid/overlay/mlml/kk,PairHybridOverlayMLMLKokkos<LMPDeviceType>);
PairStyle(hybrid/overlay/mlml/kk/device,PairHybridOverlayMLMLKokkos<LMPDeviceType>);
PairStyle(hybrid/overlay/mlml/kk/host,PairHybridOverlayMLMLKokkos<LMPHostType>);
// clang-format on
#else

#ifndef LMP_PAIR_HYBRID_OVERLAY_MLML_KOKKOS_H
#define LMP_PAIR_HYBRID_OVERLAY_MLML_KOKKOS_H

#include "pair_hybrid_overlay_mlml.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {
template<class DeviceType>
class PairHybridOverlayMLMLKokkos : public PairHybridOverlayMLML {
  friend class FixGPU;
  friend class FixIntel;
  friend class FixOMP;
  friend class Force;
  friend class Respa;
  friend class Info;
 public:
  typedef ArrayTypes<DeviceType> AT;
  PairHybridOverlayMLMLKokkos(class LAMMPS *);
  ~PairHybridOverlayMLMLKokkos();

  void settings(int, char**) override;
  void coeff(int, char **) override;
  void compute(int, int) override;
  void init_style() override;
  void allocate_mem() override;
  void resize_arrays() override;
  void modify_neighbor_list(int) override;
  void restore_neighbor_list(int) override;

 private:
   int potential_1_idx, potential_2_idx;
   int eval_1_idx, eval_2_idx;

   typename AT::t_int_1d d_keep_flags;
   typename AT::t_int_1d d_prefix_sum;
   DAT::tdual_int_1d k_ilist_copy;
   DAT::tdual_int_1d k_ilist_temp;
   DAT::tdual_double_2d k_f_copy;
   DAT::tdual_int_1d k_pot_eval_arr;
   DAT::tdual_int_1d k_iskip;

//   void init_svector() override;
//   void copy_svector(int, int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif