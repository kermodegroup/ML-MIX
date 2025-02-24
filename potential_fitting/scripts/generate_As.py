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

#create a 10x10x10 supercell of tungsten
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS
from ase.calculators.lammpslib import LAMMPSlib
from ase.io.lammpsdata import write_lammps_data
import numpy as np
from mpi4py import MPI
from lammps import lammps
from matscipy import parameter
import sys
sys.path.append('.')
import params
import os
import ase
from ase.calculators.singlepoint import SinglePointCalculator


lmps = lammps()
me = MPI.COMM_WORLD.Get_rank()
comm_world = MPI.COMM_WORLD
#isolate each rank on a seperate communicator before passing in
single_comm = comm_world.Split(color=me, key=me)
##
el = parameter('el')
a0_init = parameter('a0_init')
mass = parameter('mass')
cmds = parameter('cmds')
T = parameter('T')
damping = parameter('damping')
rseed = parameter('rseed')
n = parameter('As_n')
As_prune = parameter('As_prune')
nsteps = parameter('nsteps')
dump_freq = parameter('dump_freq')

# if me == 0:
#     if not os.path.exists('As_data'):
#         os.makedirs('As_data')
#     calc = LAMMPSlib(comm=single_comm,amendments=[f"mass 1 {mass}"],lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)
#     lattice = parameter('lattice')
#     print('optimising lattice parameter')
#     unit_cell = lattice(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
#     unit_cell.calc = calc
#     ecf = ExpCellFilter(unit_cell)
#     uc_optimise = LBFGS(ecf)
#     uc_optimise.run(fmax=0.0001)
#     a0 = unit_cell.get_cell()[0,0] #get the optimised lattice parameter
#     print('optimised lattice parameter:',a0)

#     #set up a nxnxn supercell
#     full_cell = unit_cell*(n,n,n)
#     write_lammps_data('input_struct.data',full_cell,velocities=True)
    
# #start lammps simulation
# lmps.command('clear') 
# lmps.command('dimension 3')
# lmps.command('boundary p p p')
# lmps.command('atom_style atomic')
# lmps.command('units metal')
# lmps.command(f'read_data input_struct.data')
# lmps.command(f'mass 1 {mass}')
# #lmps.command(f'mass 2 4')
# lmps.commands_list(cmds)
# lmps.command('compute pe all pe/atom')
# lmps.command('fix sf all store/force')
# lmps.command('fix 5 all nve')
# lmps.command(f'fix thermostat all langevin {T} {T} {damping} {rseed} zero yes')
# lmps.command(f'dump myDump all custom {dump_freq} dump.lammpstrj id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe')
# lmps.command(f'thermo 100')
# lmps.command(f'run {nsteps}')



if me == 0:
    traj = ase.io.read('dump.lammpstrj',index=':',parallel=False)
    new_traj = []
    print('traj length:',len(traj))
    for i,t in enumerate(traj):
        if i%As_prune != 0:
            continue
        #set all symbols to el
        t.set_chemical_symbols([el]*len(t))
        fs = [t.arrays[f'f_sf[{i}]'] for i in range(1,4)]
        pe_per_atom = t.arrays['c_pe']
        pe = pe_per_atom.sum()
        force_raw = np.concatenate(fs,axis=1)
        t.calc = SinglePointCalculator(t,energy=pe,forces=force_raw)    
        #delete the f_sf arrays
        for i in range(1,4):
            del t.arrays[f'f_sf[{i}]'] #JuLIP cant cope with these
        new_traj.append(t)

    print('refined traj length:',len(new_traj))
    ase.io.write('As_data/data_for_As.xyz',new_traj,format='extxyz',parallel=False)