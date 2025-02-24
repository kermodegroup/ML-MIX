# ----------------------------------------------------------------------
# This script is part of the ML-MIX repository.
# 
# Copyright (2025) Fraser Birks
#
# This script is licensed under the MIT License. You are free to use, modify,
# and distribute it, provided that this copyright notice and license text 
# remain intact.
#
# See the LICENSE file in the repository root for full details.
# ----------------------------------------------------------------------


import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
sys.path.append('.')
import utils
import ase.io
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from lammps import lammps
from matscipy import parameter
import params
import os
from ase.io.lammpsdata import write_lammps_data

calc_cmds = parameter('cmds')
in_file_root = parameter('in_file_root')
in_file_path = parameter('in_file_path')
base_struct_name = parameter('base_struct_name')
multi_potential = parameter('multi_potential')
dump_name = parameter('dump_name')
dump_freq = parameter('dump_freq')
n_steps = parameter('n_steps')
mlml_nevery = parameter('mlml_nevery')
rqm = parameter('rqm')
bw = parameter('bw')
damping = parameter('damping')
rblend = parameter('rblend')
dump_files = parameter('dump_files')
thermo_freq = parameter('thermo_freq',100)
load_balancing = parameter('load_balancing',True)
HeW_sim = parameter('HeW_sim',False)
tiled = parameter('tiled', False)
rcb = parameter('rcb', False)
if HeW_sim:
    mass_W = parameter('mass_W')
    mass_He = parameter('mass_He')
else:
    mass = parameter('mass')

plt.ioff()
comm = MPI.COMM_WORLD
lmps = lammps()
rank = comm.Get_rank()
print(rank)

if rank == 0:
    struct = ase.io.read(base_struct_name,parallel=False)
    selected_atom = np.where(struct.arrays['selected_atoms'] == 1)[0]
    rseed = parameter('rseed')
    if not HeW_sim:
        write_lammps_data(f'{in_file_path}/{in_file_root}.lj',struct,velocities=True)
    else:
        write_lammps_data(f'{in_file_path}/{in_file_root}.lj',struct,velocities=True,specorder=['W','He'])
else:
    selected_atom = np.zeros(1,dtype=int)
    rseed = 0
selected_atom = comm.bcast(selected_atom,root=0)
rseed = comm.bcast(rseed,root=0)

input_file = f'{in_file_path}/{in_file_root}.lj'

if HeW_sim:
    mass_cmds = [f'mass 1 {mass_W}',f'mass 2 {mass_He}']
else:
    mass_cmds = [f'mass 1 {mass}']

utils.bare_bones_setup_lammps(lmps,input_file,mass_cmds,calc_cmds,
                             sim_tstep=0.001, thermo_freq=thermo_freq,multi_potential=multi_potential)

if tiled:
    lmps.command(f'comm_style tiled')

max_comm_size = np.max((rqm,bw,rblend)) + 2
lmps.command(f'comm_modify cutoff {max_comm_size}')

lmps.command(f'group bond_atoms id {selected_atom[0]+1}')
if multi_potential:
    lmps.command(f'fix mlml_fix all mlml {mlml_nevery} {rqm} {bw} {rblend} group bond_atoms')


lmps.command('fix 5 all nve')
if load_balancing:
    if not rcb:
        lmps.command(f'fix b_fix all balance 5 1.1 shift xyz 5 1.05 weight time 1.0')
    else:
        #if rcb
        lmps.command(f'fix b_fix all balance 5 1.1 rcb weight time 1.0')

if dump_files:
    if multi_potential:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom d2_eval[1] d2_eval[2] i_potential[1] i_potential[2]')
    else:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom')    



lmps.command(f'run 6') # to invoke a balance command
lmps.command(f'run {n_steps}')