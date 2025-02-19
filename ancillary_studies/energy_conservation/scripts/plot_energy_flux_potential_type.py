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
# potential_files = [['2_10_control', '2_10'], 
#                    ['2_15_control', '2_15'], 
#                    ['3_15_control', '3_15']]
potential_files = [['2_10', '2_15', '3_15']]

studies = ['abrupt_mixing_change_rqm', 'smooth_mixing_r_tot_10']

#change this to plot smooth mixing or abrupt mixing
study = 0

if not os.path.exists('plots'):
    os.makedirs('plots')

if not os.path.exists(f'plots/{studies[study]}'):
    os.makedirs(f'plots/{studies[study]}')

for potential_set in potential_files:
    plt.figure()
    for potential_file in potential_set:
        
        print(f'{potential_file}')
        energies_path = f'./results/{potential_file}/{studies[study]}/energies.json'
        with open(energies_path, 'r') as f:
            energies_dict = json.load(f)
        r_10_energies = energies_dict['r=10']
        r_10_energies = np.array(r_10_energies)
        r_10_energies -= r_10_energies[0]
        mask = np.where(r_10_energies<600)
        r_10_energies = r_10_energies[mask]
        flux = (r_10_energies[-1] - r_10_energies[0])/(4*np.pi*(10**2)*len(r_10_energies))
        print(f'flux for {potential_file} is {flux}')
        # plot on bar chart
        plt.bar(potential_file, flux, label=f'{potential_file}')


    # Add labels and legend
    plt.xlabel('Cheap potential')
    plt.ylabel(r'Energy Flux ($\mathrm{eV/ps\AA}^2$)')
    #set y axis to log
    plt.yscale('log')
    plt.title(f'Energy flux')
    plt.tight_layout()
    plt.savefig(f'plots/{studies[study]}/energy_flux_for_{potential_set}.png')