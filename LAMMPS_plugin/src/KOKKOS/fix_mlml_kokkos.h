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
FixStyle(mlml/kk,FixLearnKokkos<LMPDeviceType>);
FixStyle(mlml/kk/device,FixLearnKokkos<LMPDeviceType>);
FixStyle(mlml/kk/host,FixLearnKokkos<LMPHostType>);
// clang-format on
#else
// clang-format off

#ifndef LMP_FIX_MLML_KOKKOS_H
#define LMP_FIX_MLML_KOKKOS_H

#include "fix_mlml.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class FixMLML : public Fix {
 public:
  FixMLMLKokkos(class LAMMPS *, int, char **);
  ~FixMLMLKokkos();

  void allocate_regions(double **) override;
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