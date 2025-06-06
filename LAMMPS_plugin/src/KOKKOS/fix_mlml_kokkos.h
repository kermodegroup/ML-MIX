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
#ifdef FIX_CLASS
// clang-format off
FixStyle(mlml/kk,FixMLMLKokkos<LMPDeviceType>);
FixStyle(mlml/kk/device,FixMLMLKokkos<LMPDeviceType>);
FixStyle(mlml/kk/host,FixMLMLKokkos<LMPHostType>);
// clang-format on
#else
// clang-format off

#ifndef LMP_FIX_MLML_KOKKOS_H
#define LMP_FIX_MLML_KOKKOS_H

#include "fix_mlml.h"
#include "neigh_list_kokkos.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {
template <class DeviceType>
class FixMLMLKokkos : public FixMLML {
 public:
  typedef ArrayTypes<DeviceType> AT;
  NeighListKokkos<DeviceType> *k_list;
  FixMLMLKokkos(class LAMMPS *, int, char **);
  ~FixMLMLKokkos();
  void init() override;
  void allocate_regions() override;
  void setup_pre_force(int) override;
  double get_x(int i, int j) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  int potential_1_idx, potential_2_idx;
  int eval_1_idx, eval_2_idx;
  double *d_eval_prev_1, *d_eval_prev_2;
  
};

}    // namespace LAMMPS_NS

#endif

    /* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
#endif
