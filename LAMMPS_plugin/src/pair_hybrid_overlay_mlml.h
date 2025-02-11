/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
#ifndef LMP_PAIR_HYBRID_OVERLAY_MLML_H
#define LMP_PAIR_HYBRID_OVERLAY_MLML_H

#include "pair_hybrid_overlay.h"

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


//   void init_svector() override;
//   void copy_svector(int, int) override;
};

}    // namespace LAMMPS_NS

#endif