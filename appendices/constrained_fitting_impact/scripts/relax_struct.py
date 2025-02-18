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
in_file_path = parameter('in_file_path')
in_file_names = parameter('in_file_names')
multi_potential = parameter('multi_potential')
mlml_nevery = parameter('mlml_nevery')
rqm = parameter('rqm')
bw = parameter('bw')
rblend = parameter('rblend')
thermo_freq = parameter('thermo_freq',1)

plt.ioff()
comm = MPI.COMM_WORLD
lmps = lammps()
rank = comm.Get_rank()
print(rank)

for name in in_file_names:
    print("Relaxing",name)
    if rank == 0: 
        struct = ase.io.read(f'{in_file_path}/{name}.xyz',parallel=False)
        central_atom_ids = np.where(struct.arrays['selected_atoms'] == 1)[0]
        rseed = parameter('rseed')
    else:
        central_atom_ids = None
        rseed = 0
    central_atom_ids = comm.bcast(central_atom_ids,root=0)
    rseed = comm.bcast(rseed,root=0)
    input_file = f'{in_file_path}/{name}.lj'
    mass_cmds = [f'mass 1 {mass}']

    utils.bare_bones_setup_lammps(lmps,input_file,mass_cmds,calc_cmds,
                                sim_tstep=0.001,thermo_freq=thermo_freq,multi_potential=multi_potential)

    max_comm_size = np.max((rqm,bw,rblend)) + 2
    lmps.command(f'comm_modify cutoff {max_comm_size}')

    central_atom_ids_str = ' '.join(map(str, central_atom_ids + 1))
    lmps.command(f'group central_atom id {central_atom_ids_str}')
    if multi_potential:
        lmps.command(f'fix mlml_fix all mlml {mlml_nevery} {rqm} {bw} {rblend} group central_atom')

    if multi_potential:
        lmps.command(f'dump myDump all custom 100 relax.lammpstrj id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom d2_eval[1] d2_eval[2] i_potential[1] i_potential[2] ix iy iz')
    else:
        lmps.command(f'dump myDump all custom 100 relax.lammpstrj id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom ix iy iz')  

    lmps.command(f'fix b_fix all balance 10 1.1 shift xyz 5 1.05 weight time 1.0')  
    lmps.command(f'run 11') # balance processors
    lmps.command("min_style fire")
    #lmps.command("min_modify norm inf")
    try:
        lmps.command("minimize 0.0 1.0e-3 2000 2000")
        if multi_potential:
            lmps.command(f'write_dump all custom {name}.lammpstrj id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom d2_eval[1] d2_eval[2] i_potential[1] i_potential[2] ix iy iz')
        else:
            lmps.command(f'write_dump all custom {name}.lammpstrj id type xs ys zs vx vy vz f_sf[1] f_sf[2] f_sf[3] c_pe_peratom ix iy iz')   
    except:
        continue

