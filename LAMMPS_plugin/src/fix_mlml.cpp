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

#include "fix_mlml.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "respa.h"
#include "update.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "comm.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "memory.h"
#include <cstring>

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixMLML::FixMLML(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  local_qm_atom_list = nullptr;
  all_qm = nullptr;
  gflag = false;
  fflag = false;
  init_flag = false;
  comm_forward = 4;
  comm_reverse = 4;
  setup_only = false;
  all_pot_one_flag = false;
  prev_ntot = -1;
  prev_qm_tot = -1;
  first_set = true;
  time_decay_hysteresis = false;
  time_decay_constant_in = 0.0;
  initial_allocation = false;
  blend_type = 0; // linear blending by default
  bool check_kwargs = false;
  set_pair_only = false;

  memory->create(core_qm_atom_idx, atom->nmax, "FixMLML: core_qm_atom_idx");
  memory->create(local_qm_atom_list, atom->nmax, "FixMLML: local_qm_atom_list");

  atom->add_callback(Atom::GROW);

  bool no_init_group = true;
  bool blend_set = false;
  // fix 1 all mlml nevery rqm bw rblend type
  if (narg < 9) utils::missing_cmd_args(FLERR, "fix mlml", error);

  nevery = utils::inumeric(FLERR,arg[3],false,lmp);
  if (nevery <= 0){
    // std::cout<<nevery<<std::endl;
    if (nevery == 0) setup_only = true;
    else error->all(FLERR,"Illegal fix mlml nevery value: {}", nevery);
  }
  
  rqm = utils::numeric(FLERR,arg[4],false,lmp);
  if (rqm < 0.0)
    error->all(FLERR,"Illegal fix mlml rqm value: {}", rqm);

  bw = utils::numeric(FLERR,arg[5],false,lmp);
  if (bw < 0.0)
    error->all(FLERR,"Illegal fix mlml bw value: {}", bw);
  
  rblend = utils::numeric(FLERR,arg[6],false,lmp);
  if (rblend < 0.0)
    error->all(FLERR,"Illegal fix mlml rblend value: {}", rblend);
  // now check for keyword arguments
  // group
  
  int iarg = 7;
  if (strcmp(arg[iarg],"group") == 0) {
    group2 = utils::strdup(arg[iarg+1]);
    igroup2 = group->find(arg[iarg+1]);
    if (igroup2 == -1)
      error->all(FLERR,"Group ID does not exist");
    group2bit = group->bitmask[igroup2];    
    gflag = true;
    iarg = 9;
    if (narg > 9) check_kwargs = true;
  // classify using the output of a different fix
  } else if (strcmp(arg[iarg], "fix_classify") == 0){
    // if (iarg + 5 > narg) error->all(FLERR,"Illegal fix mlml fix_classify command");
    fix_id = utils::strdup(arg[iarg+1]);
    // error checking for fix is done when it is needed

    // next argument is nfreq, integer (how often to check connected fix)
    nfreq = utils::inumeric(FLERR,arg[iarg+2],false,lmp);

    // pass remaining arguments as doubles
    
    // check if lb is -inf
    if (strcmp(arg[iarg+3], "-inf") == 0){
      lb = -1.0e100; // switch to flag system later if this doesn't work
    } else lb = utils::numeric(FLERR,arg[iarg+3],false,lmp);

    // check if ub is inf
    if (strcmp(arg[iarg+4], "inf") == 0){
      ub = 1.0e100; // switch to flag system later if this doesn't work
    } else ub = utils::numeric(FLERR,arg[iarg+4],false,lmp);

    // check coord_ub > coord_lb
    if (ub < lb){
      error->all(FLERR, "FixMLML: fix upper bound must be greater than lower bound");
    }
    fflag=true;
    if (narg>12){
    // now check if we are using an initialisation group
      if (strcmp(arg[iarg+5], "init_group")==0){
        no_init_group = false;
        init_flag=true;
        first_set=false;
        group2 = utils::strdup(arg[iarg+6]);
        igroup2 = group->find(group2);
        if (igroup2 == -1)
          error->all(FLERR,"Group ID does not exist");
        group2bit = group->bitmask[igroup2];
        iarg = 14;
        if (narg > 14) check_kwargs = true;
      } else {
        check_kwargs = true;
        iarg = 12;
      }
    }
    if (no_init_group){
      all_pot_one_flag = true;
      first_set = false;
    }
  } else error->all(FLERR,"FixMLML: Illegal fix mlml command");
  if (check_kwargs){
    while (true){
      if (iarg >= narg) break;
      if (strcmp(arg[iarg], "hysteresis-time") == 0){
        if (time_decay_hysteresis) error->all(FLERR, "FixMLML: hysteresis-time specified multiple times");
        // check there are two more arguments
        if (iarg+2 > narg) error->all(FLERR,"FixMLML: Illegal fix mlml on-fly command");
        time_decay_hysteresis=true;
        time_decay_constant_in = utils::numeric(FLERR,arg[iarg+1],false,lmp);
        if (time_decay_constant_in <= 0.0){
          error->all(FLERR,"FixMLML: Illegal fix mlml hysteresis-time value: {}", time_decay_constant_in);
        }
        time_decay_constant_out = utils::numeric(FLERR,arg[iarg+2],false,lmp);
        if (time_decay_constant_out <= 0.0){
          error->all(FLERR,"FixMLML: Illegal fix mlml hysteresis-time value: {}", time_decay_constant_out);
        }
        iarg += 3;
      } else if (strcmp(arg[iarg], "blend") == 0){
        if (blend_set) error->all(FLERR, "FixMLML: blend specified multiple times");
        blend_set = true;
        // check there is one more argument
        if (iarg+1 > narg) error->all(FLERR,"FixMLML: Illegal fix mlml on-fly command");
        // if argument is "linear" set blend type to 0
        // if argument is "cubic" set blend type to 1
        // else error with unrecognised blend type
        if (strcmp(arg[iarg+1], "linear") == 0){
          blend_type = 0;
        } else if (strcmp(arg[iarg+1], "cubic") == 0){
          blend_type = 1;
        } else {
          error->all(FLERR,"FixMLML: Illegal fix mlml blend type: {}", arg[iarg+1]);
        }
        iarg += 2;
      } else if (strcmp(arg[iarg], "pair_only") == 0){
        if (set_pair_only) error->all(FLERR, "FixMLML: pair_only specified multiple times");
        set_pair_only = true;
        // check there are no more arguments
        if (iarg+1 > narg) error->all(FLERR,"FixMLML: Illegal fix mlml pair_only command");
        // next value is the region to set pair only
        pair_only_region = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
        if (pair_only_region < 1 || pair_only_region > 2){
          error->all(FLERR,"FixMLML: Illegal fix mlml pair_only region: {}", pair_only_region);
        }
        iarg += 2;
      }else{
        error->all(FLERR,"FixMLML: Unrecognised fix mlml keyword: {}", arg[iarg]);
      }
    }
  }

  // create arrays
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntot = nlocal + nghost;
  memory->create(d2_eval_prev, ntot, 2, "FixMLML: d2_eval_prev");
  pair_forces_idx = -1;
  if (set_pair_only) pair_forces_idx = atom->add_custom("pair_only_forces", 1, 3, 1);


}

FixMLML::~FixMLML()
{
  if (copymode) return;
  if (local_qm_atom_list) memory->destroy(local_qm_atom_list);
  if (core_qm_atom_idx) memory->destroy(core_qm_atom_idx);
  atom->delete_callback(id, Atom::GROW);
  //if (d2_eval_prev) memory->destroy(d2_eval_prev);
}

/* ---------------------------------------------------------------------- */

int FixMLML::setmask()
{
  int mask = 0;
  mask |= FixConst::PRE_FORCE;
  mask |= FixConst::MIN_PRE_FORCE;
  if (set_pair_only) mask |= FixConst::PRE_REVERSE;
  if (!setup_only) mask |= FixConst::END_OF_STEP;
  if (!setup_only) mask |= FixConst::MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMLML::init()
{
  double max_cutoff = fmax(rqm, fmax(bw, rblend));
  neighbor->add_request(this, NeighConst::REQ_FULL)->set_cutoff(max_cutoff);
}

void FixMLML::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

void FixMLML::grow_arrays(int nmax)
{
  int ntot = atom->nlocal + atom->nghost;
  // grow the local_qm_atom_list array
  memory->grow(local_qm_atom_list, nmax, "FixMLML: local_qm_atom_list");
  memory->grow(core_qm_atom_idx, nmax, "FixMLML: core_qm_atom_idx");

}

void FixMLML::setup_pre_force(int){
  // called right at the start of simulation
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntot = nlocal + nghost;

  if (ntot > prev_ntot){
    memory->grow(d2_eval_prev, ntot, 2, "FixMLML: d2_eval_prev");
    // if we are using pair_only_forces, then grow the array
    if (set_pair_only){
      double** pair_only_forces = atom->darray[pair_forces_idx];
      memory->grow(pair_only_forces, ntot, 3, "FixMLML: pair_only_forces");
    }
    prev_ntot = ntot;
  }

  int **i2_potential = (int**)atom->extract("i2_potential");
  double **d2_eval = (double**)atom->extract("d2_eval");
  // set all arrays to 0
  for (int i = 0; i < ntot; i++){
    i2_potential[i][0] = 0;
    i2_potential[i][1] = 0;
    d2_eval[i][0] = 0.0;
    d2_eval[i][1] = 0.0;
  }

  // if no initial group set, then start all atoms as QM
  if (all_pot_one_flag){
    error->warning(FLERR, "FixMLML: fix_classify command does not have an initialisation group, all atoms will be evaluated with potential 1 until first fix evaluation");
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    if (i2_potential == nullptr || d2_eval == nullptr){
      error->all(FLERR, "FixMLML: both i2_potential and d2_eval must be allocated");
    }

    for (int i = 0; i < nlocal + nghost; i++){
      i2_potential[i][0] = 1;
      i2_potential[i][1] = 0;
      d2_eval[i][0] = 1.0;
      d2_eval[i][1] = 0.0;
    }

    first_set = false;
    // setting this to true here means that it will slowly decay from
    // all expensive to the fix atoms after the first fix evaluation
    initial_allocation = true;
  }else{
    // if there is an initial group, switch at the start
    // to doing it based on group. If not, then do it
    // with the fix
    if (init_flag){
      if (update->ntimestep > 0) {
        error->warning(FLERR, "FixMLML: Using initialisation group until first fix evaluation!");
        first_set = false;
        // setting this to true here means that it will slowly decay
        // from the initial group atoms to the fix atoms after the
        // first fix evaluation
        initial_allocation = true;
      }
      fflag = false;
      gflag = true;
    }

    this->allocate_regions();

    if (init_flag){
      fflag = true;
      gflag = false;
    }
  }
}

void FixMLML::min_setup(int)
{
  // ensure that fflag is FALSE
  if (fflag) error->all(FLERR, "FixMLML: minimisations must use group based classification");
  this->setup_pre_force(0);
}

void FixMLML::min_post_force(int)
{
  if (fflag) error->all(FLERR, "FixMLML: minimisations must use group based classification");
  this->end_of_step();
}

// pre_force, resize and zero the d2_pair_only_forces array
void FixMLML::pre_force(int /*vflag*/) {

  int ntot = atom->nlocal + atom->nghost;
  if (ntot > prev_ntot){
    memory->grow(d2_eval_prev, ntot, 2, "FixMLML: d2_eval_prev");
    // if we are using pair_only_forces, then grow the array
    if (set_pair_only){
      double** pair_only_forces = atom->darray[pair_forces_idx];
      memory->grow(pair_only_forces, ntot, 3, "FixMLML: pair_only_forces");
    }
    prev_ntot = ntot;
  }

  if (set_pair_only) {
    double** pair_only_forces = atom->darray[pair_forces_idx];
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;
    int ntot = nlocal + nghost;

    // zero the pair_only_forces array
    for (int i = 0; i < ntot; i++){
      pair_only_forces[i][0] = 0.0;
      pair_only_forces[i][1] = 0.0;
      pair_only_forces[i][2] = 0.0;
    }
  }
}

void FixMLML::setup_pre_reverse(int flag1, int flag2) {
  if (set_pair_only) {
    this->pre_reverse(flag1, flag2);
  }
}

void FixMLML::pre_reverse(int flag1, int flag2) {
  // if we are using pair_only_forces, then we need to zero
  // out all non pair_contributions to the forces according to
  // region supplied in input script

  // this should have been set in pair_hybrid_overlay_mlml
  double** pair_only_forces = atom->darray[pair_forces_idx];
  double** d2_eval = (double**)atom->extract("d2_eval");
  double** f = atom->f;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntot = nlocal + nghost;

  // iterate over all local and ghost atoms, subtracting off the fraction of
  // the non pair forces according to d2_eval
  for (int i = 0; i < ntot; i++){
    for (int j = 0; j < 3; j++){
      f[i][j] -= (d2_eval[i][pair_only_region-1]) *
                 (f[i][j] - pair_only_forces[i][j]);
    }
  }

}

void FixMLML::end_of_step()
{
  // at the end of the timestep this is called (every nevery timesteps)
  this->allocate_regions();
}


void FixMLML::allocate_regions(){
  int **i2_potential = (int**)atom->extract("i2_potential");
  double **d2_eval = (double**)atom->extract("d2_eval");
  int *ilist = list->ilist;
  int *jlist;
  int *num_neigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int inum = list->inum;
  bool atom_is_qm;
  bool just_qm;
  int n_core_qm = 0;
  int ntot = nlocal + nghost;

  bigint current_timestep = update->ntimestep;
  double dt = update->dt;

  //std::cout<<"ALLOCATING REGIONS!!!"<<std::endl;
  
  if (i2_potential == nullptr || d2_eval == nullptr){
    error->all(FLERR, "FixMLML: both i2_potential and d2_eval must be allocated");
  }

  // if we are tracking based on fix:
  if (fflag){
    // check if this is a timestep where global qm list is updated
    if (current_timestep != 0 && current_timestep % nfreq == 0){
      // if it is, get fix vector and update list
      first_set = true;
      int ifix = modify->find_fix(fix_id);
      if (ifix == -1){
        error->all(FLERR, "FixMLML: fix id does not exist");
      }
      classify_fix = modify->fix[ifix];
      classify_vec = classify_fix->vector_atom;
      if (classify_vec == nullptr){
        error->all(FLERR, "FixMLML: fix does not output compatible vector");
      }
      // communicate the global QM list to all processors
      update_global_QM_list();
    }

    // if we have not yet set any core QM atoms, 
    // then just return with i2_potential and d2_eval being
    // set in the setup_pre_force function
    if (!first_set){
      return;
    }
  }


  // reset all i2 and d2 values to 0
  for (int i = 0; i < nlocal + nghost; i++){
    i2_potential[i][0] = 0;
    i2_potential[i][1] = 0;
    d2_eval_prev[i][0] = d2_eval[i][0];
    d2_eval_prev[i][1] = d2_eval[i][1];
    d2_eval[i][0] = 0.0;
    d2_eval[i][1] = 0.0;
  }



  // send these to all ghost atoms on other processors
  comm->forward_comm(this);

  // iterate over all local atoms
  // if the atom is in the QM region, set it and all neighbours
  // within r_core to be evaluated with the QM potential
  for (int ii = 0; ii < inum; ii++){
    atom_is_qm = false;
    int i = ilist[ii];

    // if we are tracking based on groups, just check the group
    if (gflag){
      if (group2bit & atom->mask[i]){
        atom_is_qm = true;
      }
    }

    // if we are tracking based on fix, check the global QM list
    if (fflag){
      for (int jj = 0; jj < tot_qm; jj++){
        if (all_qm[jj] == atom->tag[i]){
          atom_is_qm = true;
          break;
        }
      }
    }

    // if the atom is QM, set it and all neighbours within r_core
    // to be evaluated with the QM potential
    if (atom_is_qm){
      // i2_potential[i][0] = 1;
      d2_eval[i][0] = 1.0;
      jlist = firstneigh[i];
      for (int jj = 0; jj < num_neigh[i]; jj++){
        int j = jlist[jj];
        j &= NEIGHMASK;
        if (check_cutoff(i, j, rqm)){
          // i2_potential[j][0] = 1;
          d2_eval[j][0] = 1.0;
        }
      }
    }
  }
  
  // communicate these core QM atoms to all other processors
  comm->reverse_comm(this);
  comm->forward_comm(this);

  // get ids of all core qm atoms
  int *core_qm_idx = new int[nlocal];
  for (int ii = 0; ii < inum; ii++){
    int i = ilist[ii];
    // if the atom is only QM, store the ID
    if (d2_eval[i][0] == 1.0){
      core_qm_idx[n_core_qm] = i;
      n_core_qm++;
    }
  }
  

  // create the blending region between the MM and QM regions.
  // do this by iterating over the neighbours of all core atoms
  // and setting d2_eval of those that are within the blending width to a linear blend
  for (int ii = 0; ii < n_core_qm; ii++){
    int i = core_qm_idx[ii];
    for (int jj = 0; jj < num_neigh[i]; jj++){
      jlist = firstneigh[i];
      int j = jlist[jj];
      j &= NEIGHMASK;
      if (check_cutoff(i, j, rblend)){
        d2_eval[j][0] = fmax(d2_eval[j][0], blend(i, j));
        // i2_potential[j][0] = 1;
      }
    }
  }

  // communicate QM blending atoms and atoms outside MM buffer
  comm->reverse_comm(this);
  comm->forward_comm(this);

  if (time_decay_hysteresis && initial_allocation){
    // iterate over d2_eval, setting d2_eval[i][0] to
    // ((d2_eval[i][0] - d2_eval_prev[i][0])/time_const)*dt*Nevery + d2_eval_prev[i][0]
    // this is a discretised exponential decay of d2_eval[i][0]
    for (int i = 0; i < nlocal + nghost; i++){
      double change = (d2_eval[i][0] - d2_eval_prev[i][0]);
      // only decay if change is negative
      // atoms should enter QM region instantly but leave slowly
      if (change < 0.0){
        double multiplier = (dt*static_cast<double>(nevery))/time_decay_constant_out;
        // if multiplier is greater than 1, set it to 1
        if (multiplier > 1.0){
          multiplier = 1.0;
        }
        d2_eval[i][0] = change*multiplier + d2_eval_prev[i][0];
      }
      if (change > 0.0){
        double multiplier = (dt*static_cast<double>(nevery))/time_decay_constant_in;
        // if multiplier is greater than 1, set it to 1
        if (multiplier > 1.0){
          multiplier = 1.0;
        }
        d2_eval[i][0] = change*multiplier + d2_eval_prev[i][0];
      }
      // if the value is less than threshold and change is decreasing set it to 0
      if (change < 0 && d2_eval[i][0] < 0.01){
        d2_eval[i][0] = 0.0;
      }
      // if the value is greater than threshold and change is increasing set it to 1
      if (change > 0 && d2_eval[i][0] > 0.99){
        d2_eval[i][0] = 1.0;
      }
    }
  }

  // communicate new d2_eval
  comm->reverse_comm(this);
  comm->forward_comm(this);

  // now we need to populate i2_potential
  
  // now create the QM buffer region
  // first get all QM atoms in core and blending region
  int *qm_and_blend_idx = new int[nlocal];
  int n_qm_and_blend_idx = 0;
  for (int ii = 0; ii < inum; ii++){
    int i = ilist[ii];
    if (d2_eval[i][0] > 0.0){
      i2_potential[i][0] = 1;
      qm_and_blend_idx[n_qm_and_blend_idx] = i;
      n_qm_and_blend_idx++;
    }
  }

  // now iterate over these atoms and add any neighbours within the buffer width
  // to the QM buffer region
  for (int ii = 0; ii < n_qm_and_blend_idx; ii++){
    int i = qm_and_blend_idx[ii];
    for (int jj = 0; jj < num_neigh[i]; jj++){
      jlist = firstneigh[i];
      int j = jlist[jj];
      j &= NEIGHMASK;
      if (check_cutoff(i, j, bw)){
        i2_potential[j][0] = 1;
      }
    }
  }

  // now form the MM buffer region
  // do this by iterating over all QM atoms
  // and checking if they are within the buffer width
  // of any MM atoms.

  for (int ii = 0; ii < n_qm_and_blend_idx; ii++){
    int i = qm_and_blend_idx[ii];
    if (d2_eval[i][0] == 1.0){
      just_qm = true;
      jlist = firstneigh[i];
      // iterate over it's neighbours
      for (int jj = 0; jj < num_neigh[i]; jj++){
        int j = jlist[jj];
        j &= NEIGHMASK;
        // if any neighbour within the buffer width has an MM component
        // atom i must be part of the MM buffer, so set just_qm to false
        if (check_cutoff(i, j, bw)){
          if (d2_eval[j][0] < 1.0){
            just_qm = false;
            break;
          }
        }
      }

      if (just_qm){
        // this atom is outside the MM buffer
        // so i2 should be set to 0
        // note that we are setting it to 1 (?!)
        // this is for communication reasons 
        // (highest value gets taken in communication - ghost atoms would overrule a 0)
        // at the end, we flip i2_potential[i][1] to 1-i2_potential[i][1]
        i2_potential[i][1] = 1;
      }
    }
  }


  // communicate i2_potential to all procs
  comm->reverse_comm(this);
  comm->forward_comm(this);

  // final thing is to make d2_eval[i][1] = 1.0-d2_eval[i][0] for all atoms
  for (int ii = 0; ii < inum; ii++){
    int i = ilist[ii];
    d2_eval[i][1] = 1.0 - d2_eval[i][0];

    // also flip i2_potential[1] to 1-i2_potential[1]
    i2_potential[i][1] = 1 - i2_potential[i][1];
  }
  // now forward comm to ensure d2_eval is correct on ghost atoms
  comm->forward_comm(this);

  // clean up memory
  delete[] core_qm_idx;
  delete[] qm_and_blend_idx;

  initial_allocation = true;
}

double FixMLML::get_x(int i, int j)
{
  return atom->x[i][j];
}

int FixMLML::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  int **i2_potential = (int**)atom->extract("i2_potential");
  double **d2_eval = (double**)atom->extract("d2_eval");

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = ubuf(i2_potential[i][0]).d;
    buf[m++] = ubuf(i2_potential[i][1]).d;
    buf[m++] = d2_eval[i][0];
    buf[m++] = d2_eval[i][1];
  }

  return m;
}

void FixMLML::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  int **i2_potential = (int**)atom->extract("i2_potential");
  double **d2_eval = (double**)atom->extract("d2_eval");

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    i2_potential[j][0] = std::max((int) ubuf(buf[m++]).i, i2_potential[j][0]);
    i2_potential[j][1] = std::max((int) ubuf(buf[m++]).i, i2_potential[j][1]);
    d2_eval[j][0] = fmax(buf[m++], d2_eval[j][0]);
    d2_eval[j][1] = fmax(buf[m++], d2_eval[j][1]);
  }
}

int FixMLML::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;
  int **i2_potential = (int**)atom->extract("i2_potential");
  double **d2_eval = (double**)atom->extract("d2_eval");

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = ubuf(i2_potential[j][0]).d;
    buf[m++] = ubuf(i2_potential[j][1]).d;
    buf[m++] = d2_eval[j][0];
    buf[m++] = d2_eval[j][1];
  }
  return m;
}

void FixMLML::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  int **i2_potential = (int**)atom->extract("i2_potential");
  double **d2_eval = (double**)atom->extract("d2_eval");
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    i2_potential[i][0] = (int) ubuf(buf[m++]).i;
    i2_potential[i][1] = (int) ubuf(buf[m++]).i;
    d2_eval[i][0] = buf[m++];
    d2_eval[i][1] = buf[m++];
  }
}

bool FixMLML::check_cutoff(int i, int j, double cutoff)
{
  double dx = this->get_x(i, 0) - this->get_x(j, 0);
  double dy = this->get_x(i, 1) - this->get_x(j, 1);
  double dz = this->get_x(i, 2) - this->get_x(j, 2);
  double rsq = dx*dx + dy*dy + dz*dz;
  if (rsq < cutoff*cutoff) return true;
  return false;
}

double FixMLML::blend(int i, int j)
{

  double delta[3];
  delta[0] = this->get_x(i, 0) - this->get_x(j, 0);
  delta[1] = this->get_x(i, 1) - this->get_x(j, 1);
  delta[2] = this->get_x(i, 2) - this->get_x(j, 2);
  double r = sqrt(delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2]);
  double x = r/rblend;

  if (blend_type==0) {
    return 1.0 - x; // linear blending
  } else if (blend_type==1) {
    return 1.0 - ((3*(x * x)) - (2*(x * x * x))); //cubic blending
  } else {
    error->all(FLERR, "FixMLML: unrecognised blend type");
    return 0.0;
  }
  
}


void FixMLML::update_global_QM_list()
{
  int nlocal = atom->nlocal;
  int *atom_ids = atom->tag;
  int tot_qm_loc;
  int world_size;
  tot_qm_loc = 0;

  // first store all local global qm atom indices
  for (int i = 0; i < nlocal; i++){
    if (classify_vec[i] >= lb && classify_vec[i] <= ub){
      local_qm_atom_list[tot_qm_loc] = atom_ids[i];
      tot_qm_loc++;
    }
  }

  // now gather all global qm atom indices
  // first gather the total number of qm atoms from all processors into an array
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  int *tot_qm_glob = new int[world_size];
  MPI_Allgather(&tot_qm_loc, 1, MPI_INT, tot_qm_glob, 1, MPI_INT, MPI_COMM_WORLD);
  tot_qm = 0;

  // get the total number of qm atoms  
  for (int i = 0; i < world_size; i++){
    tot_qm += tot_qm_glob[i];
  }

  // build an array of displacements for the allgatherv operation
  if (prev_qm_tot < tot_qm){
    memory->grow(all_qm, tot_qm, "FixMLML: all_qm");
    prev_qm_tot = tot_qm;
  }
  int *displs = new int[world_size];
  displs[0] = 0;
  for (int i = 1; i < world_size; i++){
    displs[i] = displs[i-1] + tot_qm_glob[i-1];
  }
  
  // now gather all global qm atom indices onto all processors
  MPI_Allgatherv(local_qm_atom_list, tot_qm_loc, MPI_INT, all_qm, tot_qm_glob, displs, MPI_INT, MPI_COMM_WORLD);

  delete[] tot_qm_glob;
  delete[] displs;
}
