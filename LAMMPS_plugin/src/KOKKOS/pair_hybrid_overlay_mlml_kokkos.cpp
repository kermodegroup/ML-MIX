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

#include "pair_hybrid_overlay_mlml_kokkos.h"

#include "atom.h"
#include "atom_kokkos.h"
#include "memory_kokkos.h"
#include "neigh_list_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair.h"
#include "respa.h"
#include "suffix.h"
#include "update.h"
#include "kokkos_type.h"
#include "kokkos.h"


#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairHybridOverlayMLMLKokkos<DeviceType>::PairHybridOverlayMLMLKokkos(LAMMPS *lmp) : PairHybridOverlayMLML(lmp) {
  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<LMPDeviceType>::space;

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
template <class DeviceType>
PairHybridOverlayMLMLKokkos<DeviceType>::~PairHybridOverlayMLMLKokkos() {
  if (copymode) return; 
  memoryKK->destroy_kokkos(k_f_copy, f_copy);
  memoryKK->destroy_kokkos(k_ilist_temp, ilist_temp);
  memoryKK->destroy_kokkos(k_ilist_copy, ilist_copy);
}

template <class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::init_style()
{
  PairHybridOverlayMLML::PairHybridOverlay::PairHybrid::init_style();
  for (int m = 0; m < nstyles; m++)
    if (styles[m]->execution_space == Host)
      lmp->kokkos->allow_overlap = 0;
}

template<class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::resize_arrays(){
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntot = nlocal + nghost;
  if (ntot > last_ntot) {
    memoryKK->grow_kokkos(k_f_copy, f_copy, ntot, 3, "pair:f_copy");
    last_ntot = ntot;
  }
  if (nlocal > last_nlocal){
    memoryKK->grow_kokkos(k_ilist_temp, ilist_temp, nlocal, "pair:ilist_temp");
    memoryKK->grow_kokkos(k_ilist_copy, ilist_copy, nlocal, "pair:ilist_copy");
    MemKK::realloc_kokkos(d_keep_flags, "mlml:keep_flags", nlocal);
    MemKK::realloc_kokkos(d_prefix_sum, "mlml:prefix_sum", nlocal);
    last_nlocal = nlocal;
  }
}

template<class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::coeff(int narg, char **arg){
  PairHybridOverlayMLML::coeff(narg, arg);
  // sync pot_eval_arr to device
  k_pot_eval_arr.modify<LMPHostType>();
  k_pot_eval_arr.sync<DeviceType>();
}


template<class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::settings(int narg, char **arg){
  // check if the first argument is "zero"
  if (strcmp(arg[0], "zero") == 0){
    error->all(FLERR, "PairHybridOverlayMLMLKokkos: zero flag is not supported");
  }
  PairHybridOverlayMLML::settings(narg, arg);
}

template<class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::allocate_mem(){
  PairHybridOverlayMLML::PairHybridOverlay::PairHybrid::allocate();
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntot = nlocal + nghost;
  memory->create(ilist_copy, nlocal, "pair:ilist_copy");
  memory->create(ilist_temp, nlocal, "pair:ilist_temp");
  memory->create(f_copy, ntot, 3, "pair:f_copy");
  memory->create(pot_eval_arr, nstyles, "pair:pot_eval_arr");
  memory->create(f_summed, 3, "pair:f_summed");
  MemKK::realloc_kokkos(d_keep_flags, "mlml:keep_flags", nlocal);
  MemKK::realloc_kokkos(d_prefix_sum, "mlml:prefix_sum", nlocal);
  memoryKK->create_kokkos(k_ilist_copy, ilist_copy, nlocal, "kokkos:ilist_copy");
  memoryKK->create_kokkos(k_ilist_temp, ilist_temp, nlocal, "kokkos:ilist_temp");
  memoryKK->create_kokkos(k_f_copy, f_copy, ntot, 3, "kokkos:f_copy");
  memoryKK->create_kokkos(k_pot_eval_arr, pot_eval_arr, nstyles, "kokkos:pot_eval_arr");
  resize_arrays();
  if (respa_enable == 1) {
    error->warning(FLERR,"PairHybridOverlayMLMLKokkos does not currently support RESPA, disabling");
    respa_enable = 0;
  }
}

template<class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::compute(int eflag, int vflag)
{
  // sync the DVECTOR_MASK and forces to device
  atomKK->sync(execution_space, DVECTOR_MASK | F_MASK);
  int i,j,m,n;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  
  // get dvector view on device
  auto d_dvector = atomKK->k_dvector.view<DeviceType>();
  auto d_fcopy = k_f_copy.view<DeviceType>();
  auto d_f = atomKK->k_f.view<DeviceType>();
  auto d_pot_eval_arr = k_pot_eval_arr.view<DeviceType>();
  int idx;
  
  resize_arrays();

  // check if no_virial_fdotr_compute is set and global component of
  //   incoming vflag = VIRIAL_FDOTR
  // if so, reset vflag as if global component were VIRIAL_PAIR
  // necessary since one or more sub-styles cannot compute virial as F dot r

  if (no_virial_fdotr_compute && (vflag & VIRIAL_FDOTR))
    vflag = VIRIAL_PAIR | (vflag & ~VIRIAL_FDOTR);

  ev_init(eflag,vflag);

  // check if global component of incoming vflag = VIRIAL_FDOTR
  // if so, reset vflag passed to substyle so VIRIAL_FDOTR is turned off
  // necessary so substyle will not invoke virial_fdotr_compute()

  int vflag_substyle;
  if (vflag & VIRIAL_FDOTR) vflag_substyle = vflag & ~VIRIAL_FDOTR;
  else vflag_substyle = vflag;

  double *saved_special = save_special();

  // check if we are running with r-RESPA using the hybrid keyword
  // if so, then crash as this is not supported with hybrid/overlay/mlml

  Respa *respa = nullptr;
  respaflag = 0;
  if (utils::strmatch(update->integrate_style,"^respa")) {
    respa = dynamic_cast<Respa *>(update->integrate);
    if (respa->nhybrid_styles > 0){
      error->all(FLERR, "pair_style hybrid/overlay/mlml: r-RESPA is not supported with hybrid/overlay/mlml/kokkos");
    }
  }

  for (m = 0; m < nstyles; m++) {
    
    set_special(m);

    // invoke compute() unless compute flag is turned off or
    // outerflag is set and sub-natomsstyle has a compute_outer() method

    modify_neighbor_list(m);

    // copy forces
    Kokkos::parallel_for("CopyForces", nlocal + nghost, KOKKOS_LAMBDA(const int i) {
      d_fcopy(i, 0) = d_f(i, 0);
      d_fcopy(i, 1) = d_f(i, 1);
      d_fcopy(i, 2) = d_f(i, 2);
    });
    
    Kokkos::fence();

    // zero original forces
    Kokkos::parallel_for("ZeroForces", nlocal + nghost, KOKKOS_LAMBDA(const int i) {
      d_f(i, 0) = 0.0;
      d_f(i, 1) = 0.0;
      d_f(i, 2) = 0.0;
    });
    
    // this ensures forces on both device and host are zeroed before
    // the compute call
    atomKK->modified(execution_space, F_MASK);
    atomKK->sync(Host, F_MASK);

    Kokkos::fence();
    
    if (styles[m]->compute_flag == 0) continue;
    atomKK->sync(styles[m]->execution_space,styles[m]->datamask_read);
    if (outerflag && styles[m]->respa_enable)
      styles[m]->compute_outer(eflag,vflag_substyle);
    else styles[m]->compute(eflag,vflag_substyle);
    atomKK->modified(styles[m]->execution_space,styles[m]->datamask_modify);

    // after the compute call, the forces are modified on device
    // but will still be 0 on host (probably)

    if ((pot_eval_arr[m] - 1) == 0) {
      idx = eval_1_idx;
    } else if ((pot_eval_arr[m] - 1) == 1) {
      idx = eval_2_idx;
    }
    Kokkos::fence();
    // for (int i = 0; i < nlocal + nghost; i++) {
    //   std::cout<<"f["<<i<<"] = ("<<atom->f[i][0]<<", "<<atom->f[i][1]<<", "<<atom->f[i][2]<<")"<<std::endl;
    // }

    // we need to combine the forces appropriately with the right eval
    // vector on device
    Kokkos::parallel_for("RestoreForces", nlocal + nghost, KOKKOS_LAMBDA(const int i) {
      d_fcopy(i, 0) += d_f(i, 0) * d_dvector(idx, i);
      d_fcopy(i, 1) += d_f(i, 1) * d_dvector(idx, i);
      d_fcopy(i, 2) += d_f(i, 2) * d_dvector(idx, i);

      d_f(i, 0) = d_fcopy(i, 0);
      d_f(i, 1) = d_fcopy(i, 1);
      d_f(i, 2) = d_fcopy(i, 2);
    });
    Kokkos::fence();

    // now we re-sync up host and device, to ensure all is smooth
    atomKK->modified(execution_space, F_MASK);
    atomKK->sync(Host, F_MASK);

    // print forces after restoring
    // for (int i = 0; i < nlocal + nghost; i++) {
    //   std::cout<<"f["<<i<<"] after restoring = ("<<atom->f[i][0]<<", "<<atom->f[i][1]<<", "<<atom->f[i][2]<<")"<<std::endl;
    // }
    // restore neigh list
    restore_neighbor_list(m);

    restore_special(saved_special);

    // global values are set to 0.0
    // atom local values are multiplied by d2_eval

    if (eflag_global) {
      eng_vdwl += 0.0; //styles[m]->eng_vdwl;
      eng_coul += 0.0; //styles[m]->eng_coul;
    }
    if (vflag_global) {
      for (n = 0; n < 6; n++) virial[n] += 0.0; //styles[m]->virial[n];
    }
    if (eflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double *eatom_substyle = styles[m]->eatom;
      for (i = 0; i < n; i++) eatom[i] += 0.0; //eatom_substyle[i]*d2_eval[i][pot_eval_arr[m]-1];
    }
    if (vflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      double **vatom_substyle = styles[m]->vatom;
      for (i = 0; i < n; i++)
        for (j = 0; j < 6; j++)
          vatom[i][j] += 0.0; // vatom_substyle[i][j]*d2_eval[i][pot_eval_arr[m]-1];
    }

    // substyles may be CENTROID_SAME or CENTROID_AVAIL

    if (cvflag_atom) {
      n = atom->nlocal;
      if (force->newton_pair) n += atom->nghost;
      if (styles[m]->centroidstressflag == CENTROID_AVAIL) {
        double **cvatom_substyle = styles[m]->cvatom;
        for (i = 0; i < n; i++)
          for (j = 0; j < 9; j++)
            cvatom[i][j] += 0.0; // cvatom_substyle[i][j]*d2_eval[i][pot_eval_arr[m]-1];
      } else {
        double **vatom_substyle = styles[m]->vatom;
        for (i = 0; i < n; i++) {
          for (j = 0; j < 6; j++) {
            cvatom[i][j] += 0.0; // vatom_substyle[i][j]*d2_eval[i][pot_eval_arr[m]-1];
          }
          for (j = 6; j < 9; j++) {
            cvatom[i][j] += 0.0; // vatom_substyle[i][j-3]*d2_eval[i][pot_eval_arr[m]-1];
          }
        }
      }
    }
  }
  atomKK->modified(execution_space, F_MASK);
  delete[] saved_special;
}

template<class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::modify_neighbor_list(int m){
  int nlocal = atom->nlocal;
  NeighList *list_m = styles[m]->list;
  NeighListKokkos<DeviceType> *k_list = static_cast<NeighListKokkos<DeviceType>*>(list_m);
  auto d_neighbors = k_list->d_neighbors;
  auto d_ilist = k_list->d_ilist;
  auto d_numneigh = k_list->d_numneigh;
  inum_copy = list_m->inum;
  int *iskip = list_m->iskip;
  auto d_type = atomKK->k_type.view<DeviceType>();
  int inumnew_local;
  inum_new = 0;
  auto keep_flags = this->d_keep_flags;
  auto prefix_sum = this->d_prefix_sum;
  auto d_ilist_temp = k_ilist_temp.view<DeviceType>();
  auto d_ilist_copy = k_ilist_copy.view<DeviceType>();
  auto d_dvector = atomKK->k_dvector.view<DeviceType>();
  int idx;

  bool skip_flag = false;
  if (iskip != nullptr) {
    skip_flag = true;
    // create a view of the iskip array on the device
    memoryKK->create_kokkos(k_iskip, iskip, atom->ntypes, "kokkos:iskip");
  } else {
    int* iskip_zeros = new int[atom->ntypes];
    for (int i = 0; i < atom->ntypes; i++) {
      iskip_zeros[i] = 0;
    }
    memoryKK->create_kokkos(k_iskip, iskip_zeros, atom->ntypes, "kokkos:iskip");
  }
  auto d_iskip = k_iskip.view<DeviceType>();


  // manually copy the data over
  Kokkos::parallel_for("CopyIlist", inum_copy, KOKKOS_LAMBDA(const int i) {
    d_ilist_copy(i) = d_ilist(i);
  });
  
  if ((pot_eval_arr[m] - 1) == 0) {
    idx = potential_1_idx;
  } else if ((pot_eval_arr[m] - 1) == 1) {
    idx = potential_2_idx;
  }
  
  Kokkos::parallel_for("zero_views", nlocal, KOKKOS_LAMBDA(const int i) {
    keep_flags(i) = 0;
    prefix_sum(i) = 0;
  });

  Kokkos::fence();
  
  Kokkos::parallel_for("mark_to_keep", nlocal, KOKKOS_LAMBDA(const int i) {
    bool keep = true;
    if (skip_flag){
      if (d_iskip[d_type(i)] == 1) keep = false;
    }
    
    if (keep){
      if (d_dvector(idx, i) > 0.999) {
        keep_flags(i) = 1;
      } else {
        keep_flags(i) = 0;
      }
    } else {
      // if the atom is skipped from the list
      // need to shift all numbers in the prefix sum (list of indices)
      // after this point back by one
      keep_flags(i) = -1;
    }
  });

  Kokkos::fence();

  Kokkos::parallel_scan("exclusive_scan", inum_copy, KOKKOS_LAMBDA(const int i, int& update, const bool final_pass) {
    int val = keep_flags(i);
    if (final_pass) prefix_sum(i) = update;
    update += val;
  });
  
  Kokkos::fence();

  Kokkos::parallel_for("write_ilist_temp", inum_copy, KOKKOS_LAMBDA(const int i) {
    if (keep_flags(i) == 1) {
      d_ilist_temp(prefix_sum(i)) = d_ilist[i];
    }
  });

  Kokkos::parallel_reduce("count_total", inum_copy, KOKKOS_LAMBDA(const int i, int& sum) {
    if (keep_flags(i) > 0) sum += keep_flags(i);
  }, inumnew_local);

  Kokkos::fence();
  list_m->inum = inumnew_local;
  inum_new = inumnew_local;
  Kokkos::parallel_for("CopyIlistTemp", inum_copy, KOKKOS_LAMBDA(const int i) {
    d_ilist(i) = d_ilist_temp(i);
  });
  if (!skip_flag) memoryKK->destroy_kokkos(k_iskip, iskip);
}

template<class DeviceType>
void PairHybridOverlayMLMLKokkos<DeviceType>::restore_neighbor_list(int m){
  NeighList *list_m = styles[m]->list;
  NeighListKokkos<DeviceType> *k_list = static_cast<NeighListKokkos<DeviceType>*>(list_m);
  list_m->inum = inum_copy;
  auto d_ilist_copy = k_ilist_copy.view<DeviceType>();
  auto d_ilist = k_list->d_ilist;

  // restore ilist
  Kokkos::parallel_for("RestoreIlist", inum_copy, KOKKOS_LAMBDA(const int i) {
    d_ilist(i) = d_ilist_copy(i);
  });
}


namespace LAMMPS_NS {
template class PairHybridOverlayMLMLKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairHybridOverlayMLMLKokkos<LMPHostType>;
#endif
}
