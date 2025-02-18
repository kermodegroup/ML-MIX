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

calc_cmds = parameter('cmds')
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
nve_sim = parameter('nve_sim',False)
mlml_nevery = parameter('mlml_nevery')
rqm = parameter('rqm')
bw = parameter('bw')
damping = parameter('damping')
rblend = parameter('rblend')
dump_files = parameter('dump_files')
thermo_freq = parameter('thermo_freq',100)

plt.ioff()
comm = MPI.COMM_WORLD
lmps = lammps()
rank = comm.Get_rank()
print(rank)

if rank == 0: 
    struct = ase.io.read(base_struct_name,parallel=False)
    rigid_atoms = np.where(struct.arrays['selected_atoms'] == 1)[0]
    rseed = parameter('rseed')
else:
    rigid_atoms = np.zeros(2,dtype=int)
    rseed = 0
rigid_atoms = comm.bcast(rigid_atoms,root=0)
rseed = comm.bcast(rseed,root=0)

input_file = f'{in_file_path}/{in_file_root}.lj'


mass_cmds = [f'mass 1 {mass}']

utils.bare_bones_setup_lammps(lmps,input_file,mass_cmds,calc_cmds,
                             sim_tstep=0.001, thermo_freq=thermo_freq,multi_potential=multi_potential)

max_comm_size = np.max((rqm,bw,rblend)) + 2
lmps.command(f'comm_modify cutoff {max_comm_size}')

lmps.command(f'group bond_atoms id {rigid_atoms[0]+1} {rigid_atoms[1]+1}')
if multi_potential:
    lmps.command(f'fix mlml_fix all mlml {mlml_nevery} {rqm} {bw} {rblend} group bond_atoms')

lmps.command(f'group non_rigid_atoms subtract all bond_atoms')

lmps.command('fix 5 non_rigid_atoms nve')

if nve_sim:
    lmps.command(f'fix 4 bond_atoms rigid single')
else:
    lmps.command(f'fix 4 bond_atoms rigid single langevin {T} {T} {damping} {rseed}')
    lmps.command(f'fix thermostat non_rigid_atoms langevin {T} {T} {damping} {rseed} zero yes')

lmps.command(f'fix b_fix all balance 10 1.1 shift xyz 5 1.05 weight time 1.0')


#dump out the bond force
lmps.command(f'dump fa_dump bond_atoms custom 1 {bond_dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3]')

if dump_files:
    if multi_potential:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom d2_eval[1] d2_eval[2] i_potential[1] i_potential[2]')
    else:
        lmps.command(f'dump myDump all custom {dump_freq} {dump_name} id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom')    



lmps.command(f'run {n_steps}')