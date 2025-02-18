from ase.calculators.lammpslib import LAMMPSlib
import ase.io
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from ase.utils.forcecurve import fit_images
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import global_plot_settings
global_plot_settings.normal_text()
pot = '0ce218cb3e'

cmds = ['pair_style hybrid/overlay pace table spline 5401',
                        f'pair_coeff * * pace {pot}.yace Si', 
                        f'pair_coeff 1 1 table {pot}_pairpot.table Si_Si']

calc = LAMMPSlib(lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)


paths = ["results/all_expensive","results/mixed_9"]
image_num = [9, 9]
forcefits = []
for j,path in enumerate(paths):
    traj = ase.io.read(f"{path}_final_images.xyz", index=':')
    for t in traj:
        print(t.get_potential_energy())
    forcefit = fit_images(traj)
    forcefits.append(forcefit)
    print(f"Path {path} done")
    print(f"energy barrier: {np.max(forcefit.energies) - np.min(forcefit.energies)}")

plt.figure()
labels = ["Expensive","ML/ML"]
colours = ['tab:blue','tab:orange']
for i,fit in enumerate(forcefits):
    plt.plot(fit.path, fit.energies, 'o', label = labels[i], color=colours[i])
    plt.plot(fit.fit_path, fit.fit_energies, color=colours[i])
    for x,y in fit.lines:
        plt.plot(x,y,'-',color=colours[i])

plt.xlabel("Reaction Coordinate (\AA)")
plt.ylabel("Energy Change (eV)")
#plt.title(f"NEB Energy Profile {path}")
plt.legend()
plt.savefig(f"plots/neb_energy_profiles.png")


