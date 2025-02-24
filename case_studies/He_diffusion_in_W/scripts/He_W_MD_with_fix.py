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
import ase.io.lammpsdata
import matplotlib.pyplot as plt
import numpy as np
from mpi4py import MPI
from lammps import lammps
from matscipy import parameter
import params

lmps = lammps()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

mass_W = parameter('mass_W')
mass_He = parameter('mass_He')
in_file_name = parameter('in_file_name')
load_file = f'{in_file_name}.xyz'
lammps_file = f'{in_file_name}.lj'
dump_name = parameter('dump_name')
T = parameter('T')
dump_freq = parameter('dump_freq')
multi_potential = parameter('multi_potential')
dump_files=parameter("dump_files")
he_dump_freq = parameter('he_dump_freq')
tstep = parameter('tstep')
burn_in_time = parameter('burn_in_time',1000)
mlml_fix = parameter('mlml_fix',True)
mlml_nevery = parameter('mlml_nevery',20)
rqm = parameter('rqm',4.0)
bw = parameter('bw',6.0)
rblend = parameter('rblend',4.0)
qm_type = parameter('qm_type',2)
weak_damping = parameter('weak_damping',5)
thermo_freq = parameter('thermo_freq',100)
therm_time = parameter('therm_time')
damping = parameter('damping')
he_dump_name = parameter('he_dump_name', 'HeDump.lammpstrj')

if multi_potential:
    cmds = parameter('cmds')
else:
    cmds = parameter('ace_cmds')

n_steps = parameter('n_steps')

# load the structure 
if rank == 0:
    rseed = parameter('rseed')
    print('Random seed:',rseed)
    struct = ase.io.read(load_file,parallel=False)
    ase.io.lammpsdata.write_lammps_data(lammps_file,struct,velocities=True,specorder=['W','He'])
else:
    rseed = 0
rseed = MPI.COMM_WORLD.bcast(rseed,root=0) 
# set up sim

mass_cmds = [f'mass 1 {mass_W}',f'mass 2 {mass_He}']

utils.bare_bones_setup_lammps(lmps,lammps_file,mass_cmds,cmds,
                             sim_tstep=0.001, thermo_freq=thermo_freq,multi_potential=multi_potential)

largest_cutoff = np.max((rqm,bw,rblend))
lmps.command(f'comm_modify cutoff {largest_cutoff+2.0}')

lmps.command('fix 5 all nve')
lmps.command('compute ca all coord/atom cutoff 2.0')
lmps.command(f'fix b_fix all balance 10 1.1 shift xyz 5 1.05 weight time 1.0')
lmps.command(f'group He_group type {qm_type}')
if multi_potential:
    lmps.command(f'fix mlml_fix all mlml {mlml_nevery} {rqm} {bw} {rblend} group He_group')
    #add weak langevin thermostat to MM region
    
if dump_files:
    if multi_potential:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom d2_eval[1] d2_eval[2] i_potential[1] i_potential[2]')
    else:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom')    

if therm_time>0:
    # add thermostat
    lmps.command(f'fix thermostat all langevin {T} {T} {damping} {rseed} zero yes')
    lmps.command(f'run {therm_time}')
    lmps.command('unfix thermostat')

if multi_potential:
    # make a group thats all minus the He group
    lmps.command(f'group notHe subtract all He_group')
    lmps.command(f'fix mlml_langevin all langevin/mlml {T} {T} {weak_damping} {rseed} 2')
    #lmps.command(f'fix thermostat notHe langevin {T} {T} {weak_damping} {rseed} zero yes')
    lmps.command(f'dump HeDump He_group custom {he_dump_freq} {he_dump_name} id type xs ys zs vx vy vz fx fy fz d2_eval[1] d2_eval[2] i_potential[1] i_potential[2]')
else:
    lmps.command(f'dump HeDump He_group custom {he_dump_freq} {he_dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3]')

lmps.command(f'run {n_steps}')

