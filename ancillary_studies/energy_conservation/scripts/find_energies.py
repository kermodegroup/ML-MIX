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

import json
from ase.calculators.lammpslib import LAMMPSlib
import ase.io
import numpy as np
from tqdm import tqdm

potential_files = ['2_10', 
                   '2_15', 
                   '3_15']

studies = ['abrupt_mixing', 'blended_mixing']

subfile_set_1 = ['r=4','r=6', 'r=8', 'r=10']
subfile_set_2 = ['r=4','r=5', 'r=6']
subfiles = [subfile_set_1, subfile_set_2]

study = 1

expensive_pot = '0ce218cb3e'

cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {expensive_pot}.yace Si', 
        f'pair_coeff 1 1 table {expensive_pot}_pairpot.table Si_Si']
mass = 28.0855
mass_cmds = [f'mass 1 {mass}']
calc = LAMMPSlib(lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)


for potential_file in potential_files:
    energy_dict = {}
    for subfile in subfiles[study]:
        path = f'{potential_file}/{studies[study]}/{subfile}/dump.lammpstrj'
        traj = ase.io.read(path, index=':')
        total_energies = []
        for t in tqdm(traj):
            t.calc = calc
            pe = t.get_potential_energy()
            vs = t.get_velocities()
            #get norm of velocities
            abs_vs = np.linalg.norm(vs, axis=1)
            ke_manual = 0.5 * np.sum(abs_vs**2) * mass
            total_energies.append(pe + ke_manual)
        energy_dict[subfile] = total_energies
    

    print('writing energies to file: ', f'{potential_file}/{studies[study]}/energies.json')
    with open(f'{potential_file}/{studies[study]}/energies.json', 'w') as f:
        json.dump(energy_dict, f)
