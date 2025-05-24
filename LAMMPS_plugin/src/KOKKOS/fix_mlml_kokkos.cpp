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
#include "neigh_list_filter_kokkos.h"
#include "neighbor.h"
#include "comm.h"
#include "domain.h"
#include "group.h"
#include "modify.h"
#include "memory.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "respa.h"
#include "update.h"
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
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;
  int flag, cols;
  
  potential_1_idx = atomKK->find_custom("potential_1", flag, cols);
  if (potential_1_idx == -1){
    error->all(FLERR, "FixMLMLKokkos: d_potential_1 must be defined");
  }
  if (flag != 1 || cols != 0) {
    error->all(FLERR, "FixMLMLKokkos: d_potential_1 must be a single column double vector");
  }

  potential_2_idx = atomKK->find_custom("potential_2", flag, cols);
  if (potential_2_idx == -1){
    error->all(FLERR, "FixMLMLKokkos: d_potential_2 must be defined");
  }
  if (flag != 1 || cols != 0) {
    error->all(FLERR, "FixMLMLKokkos: d_potential_2 must be a single column double vector");
  }

  eval_1_idx = atomKK->find_custom("eval_1", flag, cols);
  if (eval_1_idx == -1){
    error->all(FLERR, "FixMLMLKokkos: d_eval_1 must be defined");
  }
  if (flag != 1 || cols != 0) {
    error->all(FLERR, "FixMLMLKokkos: d_eval_1 must be a single column double vector");
  }

  eval_2_idx = atomKK->find_custom("eval_2", flag, cols);
  if (eval_2_idx == -1){
    error->all(FLERR, "FixMLMLKokkos: d_eval_2 must be defined");
  }
  if (flag != 1 || cols != 0) {
    error->all(FLERR, "FixMLMLKokkos: d_eval_2 must be a single column double vector");
  }

}

template<class DeviceType>
FixMLMLKokkos<DeviceType>::~FixMLMLKokkos(){
  if (copymode) return;
}

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::init()
{
  FixMLML::init();
  // currently just need neigh list on host, so none of this logic is needed

  // auto request = neighbor->find_request(this);
  // request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
  //                          !std::is_same_v<DeviceType,LMPDeviceType>);
  // request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  // for now, request a neigh list on the host
  // request->set_kokkos_host(true);
  // request->set_kokkos_device(false);
  // request->enable_full();
  // std::cout<<"FixMLMLKokkos::init() finished"<<std::endl;
}

// template<class DeviceType>
// void FixMLMLKokkos<DeviceType>::init_list(int id, NeighList *ptr)
// {
//   std::cout<<"FixMLMLKokkos::init_list() called"<<std::endl;
//   FixMLML::init_list(id, ptr);
//   k_list = static_cast<NeighListKokkos<DeviceType>*>(ptr);
//   std::cout<<"FixMLMLKokkos::init_list() finished"<<std::endl;
// }

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::setup_pre_force(int){
  atomKK->sync(Host, DVECTOR_MASK);
  atomKK->modified(Host, DVECTOR_MASK);
  // called right at the start of simulation
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntot = nlocal + nghost;
  prev_ntot = ntot;
  memory->create(d_eval_prev_1, ntot, "FixMLMLKokkos: d_eval_prev_1");
  memory->create(d_eval_prev_2, ntot, "FixMLMLKokkos: d_eval_prev_2");

  // if no initial group set, then start all atoms as QM
  if (all_pot_one_flag){
    error->warning(FLERR, "FixMLMLKokkos: fix_classify command does not have an initialisation group, all atoms will be evaluated with potential 1 until first fix evaluation");
    auto h_dvector = atomKK->k_dvector.h_view;
    int nlocal = atom->nlocal;
    int nghost = atom->nghost;

    for (int i = 0; i < nlocal + nghost; i++){
      h_dvector(potential_1_idx, i) = 1.0;
      h_dvector(potential_2_idx, i) = 0.0;
      h_dvector(eval_1_idx, i) = 1.0;
      h_dvector(eval_2_idx, i) = 0.0;
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
        error->warning(FLERR, "FixMLMLKokkos: Using initialisation group until first fix evaluation!");
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

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::allocate_regions(){
  // mark what we are reading and modifying
  atomKK->sync(Host, DVECTOR_MASK | X_MASK);
  atomKK->modified(Host, DVECTOR_MASK);
  auto h_dvector = atomKK->k_dvector.h_view;

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

  auto h_x = atomKK->k_x.h_view;

  bigint current_timestep = update->ntimestep;
  double dt = update->dt;


  if (prev_ntot<ntot){
    prev_ntot = ntot;
    memory->grow(d_eval_prev_1, ntot, "FixMLMLKokkos: d_eval_prev_1");
    memory->grow(d_eval_prev_2, ntot, "FixMLMLKokkos: d_eval_prev_2");
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
    h_dvector(potential_1_idx, i) = 0.0;
    h_dvector(potential_2_idx, i) = 0.0;
    d_eval_prev_1[i] = h_dvector(eval_1_idx, i);
    d_eval_prev_2[i] = h_dvector(eval_2_idx, i);
    h_dvector(eval_1_idx, i) = 0.0;
    h_dvector(eval_2_idx, i) = 0.0;
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
      h_dvector(eval_1_idx, i) = 1.0;
      jlist = firstneigh[i];
      for (int jj = 0; jj < num_neigh[i]; jj++){
        int j = jlist[jj];
        j &= NEIGHMASK;
        if (check_cutoff(i, j, rqm)){
          // i2_potential[j][0] = 1;
          h_dvector(eval_1_idx, j) = 1.0;
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
    if (h_dvector(eval_1_idx, i) == 1.0){
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
        h_dvector(eval_1_idx, j) = fmax(h_dvector(eval_1_idx, j), blend(i, j));
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
      double change = (h_dvector(eval_1_idx, i) - d_eval_prev_1[i]);
      // only decay if change is negative
      // atoms should enter QM region instantly but leave slowly
      if (change < 0.0){
        double multiplier = (dt*static_cast<double>(nevery))/time_decay_constant_out;
        // if multiplier is greater than 1, set it to 1
        if (multiplier > 1.0){
          multiplier = 1.0;
        }
        h_dvector(eval_1_idx, i) = change*multiplier + d_eval_prev_1[i];
      }
      if (change > 0.0){
        double multiplier = (dt*static_cast<double>(nevery))/time_decay_constant_in;
        // if multiplier is greater than 1, set it to 1
        if (multiplier > 1.0){
          multiplier = 1.0;
        }
        h_dvector(eval_1_idx, i) = change*multiplier + d_eval_prev_1[i];
      }
      // if the value is less than threshold and change is decreasing set it to 0
      if (change < 0 && h_dvector(eval_1_idx, i) < 0.01){
        h_dvector(eval_1_idx, i) = 0.0;
      }
      // if the value is greater than threshold and change is increasing set it to 1
      if (change > 0 && h_dvector(eval_1_idx, i) > 0.99){
        h_dvector(eval_1_idx, i) = 1.0;
      }
    }
  }

  // communicate new d2_eval
  comm->reverse_comm(this);
  comm->forward_comm(this);
  
  // now create the QM buffer region
  // first get all QM atoms in core and blending region
  int *qm_and_blend_idx = new int[nlocal];
  int n_qm_and_blend_idx = 0;
  for (int ii = 0; ii < inum; ii++){
    int i = ilist[ii];
    if (h_dvector(eval_1_idx, i) > 0.0){
      h_dvector(potential_1_idx, i) = 1.0;
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
        h_dvector(potential_1_idx, j) = 1.0;
      }
    }
  }

  // now form the MM buffer region
  // do this by iterating over all QM atoms
  // and checking if they are within the buffer width
  // of any MM atoms.

  for (int ii = 0; ii < n_qm_and_blend_idx; ii++){
    int i = qm_and_blend_idx[ii];
    if (h_dvector(eval_1_idx, i) == 1.0){
      just_qm = true;
      jlist = firstneigh[i];
      // iterate over it's neighbours
      for (int jj = 0; jj < num_neigh[i]; jj++){
        int j = jlist[jj];
        j &= NEIGHMASK;
        // if any neighbour within the buffer width has an MM component
        // atom i must be part of the MM buffer, so set just_qm to false
        if (check_cutoff(i, j, bw)){
          if (h_dvector(eval_1_idx, j) < 1.0){
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
        h_dvector(potential_2_idx, i) = 1.0;
      }
    }
  }


  // communicate i2_potential to all procs
  comm->reverse_comm(this);
  comm->forward_comm(this);

  // final thing is to make d2_eval[i][1] = 1.0-d2_eval[i][0] for all atoms
  for (int ii = 0; ii < inum; ii++){
    int i = ilist[ii];
    h_dvector(eval_2_idx, i) = 1.0 - h_dvector(eval_1_idx, i);

    // also flip i2_potential[1] to 1-i2_potential[1]
    h_dvector(potential_2_idx, i) = 1.0 - h_dvector(potential_2_idx, i);
  }
  // now forward comm to ensure d2_eval is correct on ghost atoms
  comm->forward_comm(this);

  // clean up memory
  delete[] core_qm_idx;
  delete[] qm_and_blend_idx;

  initial_allocation = true;
}

template<class DeviceType>
int FixMLMLKokkos<DeviceType>::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  auto h_dvector = atomKK->k_dvector.h_view;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = h_dvector(potential_1_idx, i);
    buf[m++] = h_dvector(potential_2_idx, i);
    buf[m++] = h_dvector(eval_1_idx, i);
    buf[m++] = h_dvector(eval_2_idx, i);
  }

  return m;
}

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  auto h_dvector = atomKK->k_dvector.h_view;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    h_dvector(potential_1_idx, j) = fmax(buf[m++], h_dvector(potential_1_idx, j));
    h_dvector(potential_2_idx, j) = fmax(buf[m++], h_dvector(potential_2_idx, j));
    h_dvector(eval_1_idx, j) = fmax(buf[m++], h_dvector(eval_1_idx, j));
    h_dvector(eval_2_idx, j) = fmax(buf[m++], h_dvector(eval_2_idx, j));
  }
}

template<class DeviceType>
int FixMLMLKokkos<DeviceType>::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i,j,m;
  auto h_dvector = atomKK->k_dvector.h_view;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = h_dvector(potential_1_idx, j);
    buf[m++] = h_dvector(potential_2_idx, j);
    buf[m++] = h_dvector(eval_1_idx, j);
    buf[m++] = h_dvector(eval_2_idx, j);
  }
  return m;
}

template<class DeviceType>
void FixMLMLKokkos<DeviceType>::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;
  auto h_dvector = atomKK->k_dvector.h_view;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    h_dvector(potential_1_idx, i) = buf[m++];
    h_dvector(potential_2_idx, i) = buf[m++];
    h_dvector(eval_1_idx, i) = buf[m++];
    h_dvector(eval_2_idx, i) = buf[m++];
  }
}

template<class DeviceType>
double FixMLMLKokkos<DeviceType>::get_x(int i, int j)
{
  return atomKK->k_x.h_view(i, j);
}

namespace LAMMPS_NS {
template class FixMLMLKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixMLMLKokkos<LMPHostType>;
#endif
}
