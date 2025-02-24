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
sys.path.append('.')
import matplotlib.pyplot as plt
import global_plot_settings
global_plot_settings.bigger_text()
import json
import os
import numpy as np

# iterate through the potential files and plot the bond forces and errors
potential_files = [['2_10'], 
                   ['2_15'], 
                   ['3_15']]
subfile_set_1 = ['r=10']
subfile_set_2 = ['r=6', 'r=5',  'r=4']
blending_widths = ['r_blend=4', 'r_blend=5', 'r_blend=6']
subfiles = [subfile_set_1, subfile_set_2]

studies = ['abrupt_mixing_change_rqm', 'smooth_mixing_r_tot_10']

if not os.path.exists('plots'):
    os.makedirs('plots')

if not os.path.exists(f'plots/comparison'):
    os.makedirs(f'plots/comparison')

for potential_set in potential_files:
    plt.figure()
    for potential_file in potential_set:
        for j, study in enumerate(studies):
            print(f'{potential_file}')
            energies_path = f'./results/{potential_file}/{study}/energies.json'
            with open(energies_path, 'r') as f:
                energies_dict = json.load(f)
            for i, subfile in enumerate(subfiles[j]):
                r_energies = energies_dict[subfile]
                r_energies = np.array(r_energies)
                r_energies -= r_energies[0]
                mask = np.where(r_energies<600)
                r_energies = r_energies[mask]
                flux = (r_energies[-1] - r_energies[0])/(4*np.pi*(10**2)*len(r_energies))
                print(f'flux for {study}, {subfile}, {potential_file} is {flux}')
                # plot on bar chart
                if study == 'abrupt_mixing_change_rqm':
                    plt.bar('r_blend=0', flux)
                else:
                    plt.bar(blending_widths[i], flux, label=f'{study}')


        # Add labels and legend
        plt.xlabel(r'Blending width (\AA)')
        plt.ylabel(r'Energy Flux ($\mathrm{eV/ps\AA}^2$)')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
        #plt.title(f'Energy flux vs blending width for {potential_file}')
        plt.tight_layout()
        plt.savefig(f'plots/comparison/energy_flux_blend_width_for_{potential_file}.png')