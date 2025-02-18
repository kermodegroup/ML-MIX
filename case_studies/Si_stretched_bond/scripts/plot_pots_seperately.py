import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import matplotlib.pyplot as plt
import global_plot_settings
import json
import os
global_plot_settings.normal_text()
# load up reference bond force and error

ref_bond_force_path = './results/NVE_ref/bond_forces.json'
ref_bond_force_error_path = './results/NVE_ref/bond_force_errors.json'

with open(ref_bond_force_path, 'r') as f:
    ref_bond_forces = json.load(f)

ref_bond_force = ref_bond_forces['.']
with open(ref_bond_force_error_path, 'r') as f:
    ref_bond_force_errors = json.load(f)
ref_bond_error = ref_bond_force_errors['.']

x_vals = [4,6,8,10]


potential_files = [['mixed']]

subfile_set_1 = ['r=0','r=4','r=6', 'r=8', 'r=10']
subfiles = [subfile_set_1]
study = 0
xsets = [[0,4,6,8,10]]

if not os.path.exists('plots'):
    os.makedirs('plots')

for potential_set in potential_files:
    plt.figure()

    for potential_file in potential_set:
        print(f'{potential_file}')
        forces_path = f'./results/{potential_file}/bond_forces.json'
        error_path = f'./results/{potential_file}/bond_force_errors.json'
        with open(forces_path, 'r') as f:
            forces_dict = json.load(f)
        with open(error_path, 'r') as f:
            error_dict = json.load(f)

        forces = [forces_dict[subfile] for subfile in subfiles[study]]
        errors = [error_dict[subfile] for subfile in subfiles[study]]

        plt.errorbar(xsets[study], forces, yerr=errors, fmt='o-', label=f'{potential_file}')

    # Plot the reference value as a dotted horizontal line
    plt.axhline(y=ref_bond_force, color='r', linestyle='--', label='Ref')

    # Plot the ribbon around the reference value
    ref_x_vals = [-5, 12]
    plt.fill_between(ref_x_vals, ref_bond_force - ref_bond_error, ref_bond_force + ref_bond_error, color='r', alpha=0.2)

    # Add labels and legend
    plt.xlabel(r'R ($\AA$)')
    plt.ylabel('Average Bond Force eV/A')
    plt.title(f'Average Bond Force vs Expensive Potential Radius')
    plt.tight_layout()
    plt.xlim([xsets[study][0]-1, xsets[study][-1]+1])
    #plt.ylim([-3, 8])
    plt.legend()
    plt.savefig(f'plots/bond_forces_vs_r_for_{potential_set}.png')