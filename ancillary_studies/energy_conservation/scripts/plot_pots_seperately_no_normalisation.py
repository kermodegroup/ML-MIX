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


import os
import sys
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import matplotlib.pyplot as plt
import global_plot_settings
global_plot_settings.bigger_text()
import json
import os
import numpy as np
# load up reference bond force and error

# ref_bond_force_pathcode = './results/reference/abrupt_mixing_change_rqm/bond_forces.json'
# ref_bond_force_error_path = './results/reference/abrupt_mixing_change_rqm/bond_force_errors.json'

# with open(ref_bond_force_path, 'r') as f:
#     ref_bond_forces = json.load(f)

# ref_bond_force = ref_bond_forces['.']
# with open(ref_bond_force_error_path, 'r') as f:
#     ref_bond_force_errors = json.load(f)
# ref_bond_error = ref_bond_force_errors['.']

x_vals = [4,6,8,10]

# iterate through the potential files and plot the bond forces and errors
potential_files = [['2_10'], 
                   ['2_15'], 
                   ['3_15']]

studies = ['abrupt_mixing_change_rqm', 'smooth_mixing_r_tot_10']

subfile_set_1 = ['r=4','r=6', 'r=8', 'r=10']
subfile_set_2 = ['r=4','r=5', 'r=6']
subfiles = [subfile_set_1, subfile_set_2]

#change this to plot smooth mixing or abrupt mixing
study = 0

if not os.path.exists('plots'):
    os.makedirs('plots')

if not os.path.exists(f'plots/{studies[study]}'):
    os.makedirs(f'plots/{studies[study]}')

for potential_set in potential_files:

    for potential_file in potential_set:
        plt.figure()
        print(f'{potential_file}')
        energies_path = f'./results/{potential_file}/{studies[study]}/energies.json'
        with open(energies_path, 'r') as f:
            energies_dict = json.load(f)

        for subfile in subfiles[study]:
            energies_dict[subfile] = np.array(energies_dict[subfile])
            energies_dict[subfile] -= energies_dict[subfile][0]
            mask = np.where(energies_dict[subfile]<10000) #remove the large values
            energies_dict[subfile] = energies_dict[subfile][mask]
            plt.plot(energies_dict[subfile], label=f'{subfile}')

    
        # Add labels and legend
        plt.xlabel('time (ps)')
        plt.ylabel('Energy Change (eV)')
        plt.title(f'Energy Change vs Time for {potential_file}')
        plt.tight_layout()
        plt.legend()
        plt.savefig(f'plots/{studies[study]}/energy_vs_time_for_{potential_set}.png')