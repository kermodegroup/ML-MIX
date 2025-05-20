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

#ifndef LMP_FIX_MLML_H
#define LMP_FIX_MLML_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMLML : public Fix {
 public:
  FixMLML(class LAMMPS *, int, char **);
  ~FixMLML();

  int setmask() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void end_of_step() override;
  void setup_pre_force(int) override;
  void min_setup(int) override;
  void min_post_force(int) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  void grow_arrays(int) override;

  bool check_cutoff(double *, double *, double);
  double blend(double *, double *);
  virtual void allocate_regions(double **);
  void update_global_QM_list();

 protected:
  class NeighList *list;
  // class NeighList *rqm_list;
  // class NeighList *rbl_list;
  // class NeighList *bw_list;

  // class NeighRequest *requestQM;
  // class NeighRequest *requestBL;
  // class NeighRequest *requestB;

  char *group2;
  char *fix_id;
  int igroup2, group2bit;
  int nfreq;
  int blend_type;
  
  double lb, ub;

  bool gflag, fflag, setup_only, init_flag, all_pot_one_flag;
  bool first_set;

  double dtv, dtf;
  double *step_respa;
  int mass_require, is_type, type_val, prev_ntot, prev_qm_tot;
  int tot_qm;
  double rqm, bw, rblend;

  Fix *classify_fix;
  double *classify_vec;
  int *local_qm_atom_list, *all_qm, *core_qm_atom_idx;

  double** d2_eval_prev;
  bool time_decay_hysteresis;
  double time_decay_constant_in, time_decay_constant_out;
  bool initial_allocation;

};

}    // namespace LAMMPS_NS

#endif

    /* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
