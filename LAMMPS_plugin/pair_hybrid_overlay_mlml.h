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

#ifndef LMP_PAIR_HYBRID_OVERLAY_MLML_H
#define LMP_PAIR_HYBRID_OVERLAY_MLML_H

#include "pair_hybrid_overlay.h"
#include "args_manager.h"

namespace LAMMPS_NS {

class PairHybridOverlayMLML : public PairHybridOverlay {
 public:
  PairHybridOverlayMLML(class LAMMPS *);
  ~PairHybridOverlayMLML();


  void settings(int, char**) override;
  void coeff(int, char **) override;
  void compute(int, int) override;
  void allocate_mem();
  void resize_arrays();

  ArgsManager manager;

 protected:
  int* pot_eval_arr;
  double** f_copy;
  double* f_summed;
  int* ilist_temp;
  int* ilist_copy;
  bool zero_flag;

  int last_ntot;
  int inum_copy;
  int inum_new;
  int last_nlocal;
  int style_counter;
  bool store_args;
  bool on_fly_flag;
  int fit_pot_tstep;
  char *path_to_mlmix;
  char *path_to_config;



//   void init_svector() override;
//   void copy_svector(int, int) override;
};

}    // namespace LAMMPS_NS

#endif