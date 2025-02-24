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

sys.path.append('.')
sys.path.append(f'{script_path}/../../../utils')

import utils
import ase.io.lammpsdata
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from lammps import lammps
from matscipy import parameter
import params
import os

cmds_mix = parameter('cmds_mix')
cmds_just_expensive = parameter('cmds_just_expensive')
mass = parameter('mass')
in_file_root = parameter('in_file_root')
in_file_path = parameter('in_file_path')
base_struct_name = parameter('base_struct_name')
multi_potential = parameter('multi_potential')
dump_name = parameter('dump_name')
bond_dump_name = parameter('bond_dump_name')
T = parameter('T')
dump_freq = parameter('dump_freq')
n_steps = parameter('n_steps')
rqm = parameter('rqm')
bw = parameter('bw')
lb = parameter('lb')
ub = parameter('ub')
damping = parameter('damping')
rblend = parameter('rblend')
dump_files = parameter('dump_files')
thermo_freq = parameter('thermo_freq',100)

coord_nevery = parameter('coord_nevery')
coord_nrepeat = parameter('coord_nrepeat')
coord_nfreq = parameter('coord_nfreq')

therm_time = parameter('therm_time')

mlml_nevery = parameter('mlml_nevery')

plt.ioff()
comm = MPI.COMM_WORLD
lmps = lammps()
rank = comm.Get_rank()
print(rank)

if rank == 0: 
    rseed = parameter('rseed')
    struct = ase.io.read(f'{base_struct_name}',parallel=False)
    ase.io.lammpsdata.write_lammps_data(f'{in_file_path}/{in_file_root}.lj',struct,atom_style='atomic')
else:
    rseed = 0
rseed = comm.bcast(rseed,root=0)

input_file = f'{in_file_path}/{in_file_root}.lj'


mass_cmds = [f'mass 1 {mass}']
if multi_potential:
    calc_cmds = cmds_mix
else:
    calc_cmds = cmds_just_expensive

utils.bare_bones_setup_lammps(lmps,input_file,mass_cmds,calc_cmds,
                             sim_tstep=0.001, thermo_freq=thermo_freq,multi_potential=multi_potential)

largest_cutoff = np.max((rqm,bw,rblend))
lmps.command(f'comm_modify cutoff {largest_cutoff+2.0}')


lmps.command('fix 5 all nve')
lmps.command('compute ca all coord/atom cutoff 2.0')
lmps.command(f'fix av_ca all ave/atom {coord_nevery} {coord_nrepeat} {coord_nfreq} c_ca') #5 20
lmps.command(f'fix b_fix all balance 10 1.1 shift xyz 5 1.05 weight time 1.0')

if multi_potential:
    lmps.command(f'fix mlml_fix all mlml {mlml_nevery} {rqm} {bw} {rblend} fix_classify av_ca {coord_nfreq} {lb} {ub}')
    #add weak langevin thermostat to MM region
    
if dump_files:
    if multi_potential:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom d2_eval[1] d2_eval[2] i_potential[1] i_potential[2] f_av_ca')
    else:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom f_av_ca')    

if therm_time>0:
    # add thermostat
    lmps.command(f'fix thermostat all langevin {T} {T} {damping} {rseed} zero yes')
    lmps.command(f'run {therm_time}')
    lmps.command('unfix thermostat')

if multi_potential:
    lmps.command(f'fix mlml_langevin all langevin/mlml {T} {T} 2.0 12345 2')
    #lmps.command(f'fix thermostat all langevin {T} {T} 3.0 12345 zero yes')


lmps.command(f'run {n_steps}')