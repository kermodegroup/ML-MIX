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

rvals = [[4,6,8,10], [10,10,10]]

# iterate through the potential files and plot the bond forces and errors
potential_files = [['2_10'], 
                   ['2_15'], 
                   ['3_15'],]

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

        for i, subfile in enumerate(subfiles[study]):
            energies_dict[subfile] = np.array(energies_dict[subfile])
            energies_dict[subfile] -= energies_dict[subfile][0]
            mask = np.where(energies_dict[subfile]<10000) #remove the large values
            energies_dict[subfile] = energies_dict[subfile][mask]
            plt.plot(energies_dict[subfile]/(4*np.pi*rvals[study][i]**2), label=f'{subfile}')

    
        # Add labels and legend
        plt.xlabel('time (ps)')
        plt.ylabel(r'Energy Change/Surface Area ($\mathrm{eV/\AA}^{2}$)')
        plt.title(f'Energy Change vs Time for {potential_file}')
        plt.tight_layout()
        plt.legend()
        plt.savefig(f'plots/{studies[study]}/energy_vs_time_for_{potential_set}_normalised.png')