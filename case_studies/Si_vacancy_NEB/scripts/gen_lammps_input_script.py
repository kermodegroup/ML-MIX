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
in_file_name = parameter('in_file_name')
final_image_name = parameter('final_image_name')
multi_potential = parameter('multi_potential')
mlml_nevery = parameter('mlml_nevery')
rqm = parameter('rqm')
bw = parameter('bw')
rblend = parameter('rblend')
thermo_freq = parameter('thermo_freq',1)
nimages = parameter('nimages')

plt.ioff()
comm = MPI.COMM_WORLD
# lmps = lammps()
class LammpsMock:
    def __init__(self, input_file=f'neb_input.in'):
        self.input_file = input_file
        with open(self.input_file, 'w') as f:
            f.write("# LAMMPS input file\n")

    def command(self, cmd):
        with open(self.input_file, 'a') as f:
            f.write(cmd + '\n')

    def commands_list(self, cmds):
        for cmd in cmds:
            self.command(cmd)
lmps = LammpsMock()
rank = comm.Get_rank()
assert comm.Get_size() == 1

print(rank)
if rank == 0: 
    struct = ase.io.read(f'{in_file_path}/{in_file_name}.xyz',parallel=False)
    central_atom_id = np.where(struct.arrays['selected_atoms'] == 1)[0]
    print(central_atom_id)
    rseed = parameter('rseed')
else:
    central_atom_id = 0
    rseed = 0
central_atom_id = comm.bcast(central_atom_id,root=0)
rseed = comm.bcast(rseed,root=0)
input_file = f'{in_file_path}/{in_file_name}.lj'
final_file = f'{in_file_path}/{final_image_name}.neb_file'
mass_cmds = [f'mass 1 {mass}']

utils.bare_bones_setup_lammps(lmps,input_file,mass_cmds,calc_cmds,
                            sim_tstep=0.001,thermo_freq=thermo_freq,multi_potential=multi_potential)

lmps.command("unfix sf")
max_comm_size = np.max((rqm,bw,rblend)) + 2
lmps.command(f'comm_modify cutoff {max_comm_size}')

central_atom_ids = ' '.join(map(str, central_atom_id + 1))
lmps.command(f'group central_atom id {central_atom_ids}')
if multi_potential:
    lmps.command(f'fix mlml_fix all mlml {mlml_nevery} {rqm} {bw} {rblend} group central_atom')

 

lmps.command(f'fix b_fix all balance 10 1.1 shift xyz 5 1.05 weight time 1.0')  

#lmps.command(f'run 11') # balance processors
lmps.command("min_style fire")
lmps.command("min_modify norm inf")
lmps.command("fix neb_fix all neb 10.0")

# set up dump file
# make universe style variable 
lmps.command(f'variable ufile universe ' + ' '.join([str(i) for i in range(nimages)]))
if multi_potential:
    lmps.command('dump myDump all custom 50 ${ufile}.lammpstrj id type xs ys zs vx vy vz fx fy fz c_pe_peratom d2_eval[1] d2_eval[2] i_potential[1] i_potential[2] ix iy iz')
else:
    lmps.command('dump myDump all custom 50 ${ufile}.lammpstrj id type xs ys zs vx vy vz fx fy fz c_pe_peratom ix iy iz') 


#lmps.command("minimize 0.0 1.0e-5 10000 10000")
lmps.command(f"neb 0.0 0.0005 1000 500 1 final {final_file}")

