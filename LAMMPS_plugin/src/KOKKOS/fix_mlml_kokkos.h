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
#include "kokkos_type.h"

namespace LAMMPS_NS {
template <class DeviceType>
class FixMLMLKokkos : public FixMLML {
 public:
  FixMLMLKokkos(class LAMMPS *, int, char **);
  ~FixMLMLKokkos();
  void init() override;
  void init_list(int, class NeighList *) override;
  void allocate_regions() override;
  double get_x(int, int) override;
 private:
  typename ArrayTypes<LMPHostType>::t_x_array x;
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
