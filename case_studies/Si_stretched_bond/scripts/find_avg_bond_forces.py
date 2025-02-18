import json
from ase.calculators.lammpslib import LAMMPSlib
import ase.io
import numpy as np
from tqdm import tqdm

folder_set = ['results/NVE_ref', 'results/mixed']

reference_subfolder = ['.']
potential_subfolders = ['r=0', 'r=4','r=6', 'r=8', 'r=10']

file_sets = [reference_subfolder, potential_subfolders]

expensive_pot = '0ce218cb3e'

cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {expensive_pot}.yace Si', 
        f'pair_coeff 1 1 table {expensive_pot}_pairpot.table Si_Si']
mass = 28.0855
mass_cmds = ['mass 1 28.0855']
calc = LAMMPSlib(lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)


for i, potential in enumerate(folder_set): #reference_file:
    print(f'{potential}')
    forces_dict = {}
    error_dict = {}
    for subfile in file_sets[i]:
        print(f'{subfile}')
        path = f'{potential}/{subfile}/stretched_bond_dump.lammpstrj'
        traj = ase.io.read(path, index=':')
        bond_forces = []
        for i,t in enumerate(traj):
            #v = t.get_velocities()
            positions = t.get_positions()
            vector_between_atoms = positions[0,:] - positions[1,:]
            unit_vector = vector_between_atoms/np.linalg.norm(vector_between_atoms)
            forces = np.hstack((t.arrays['f_sf[1]'],t.arrays['f_sf[2]'],t.arrays['f_sf[3]']))
            force_1_along_bond = -np.dot(forces[0,:],unit_vector)
            force_2_along_bond = np.dot(forces[1,:],unit_vector)
            sum_force_along_bond = force_1_along_bond + force_2_along_bond
            bond_forces.append(sum_force_along_bond)
        mean_bond_force = np.mean(bond_forces)
        err_mbf = np.std(bond_forces)/np.sqrt(len(bond_forces))
        print(f'combined force is {mean_bond_force} +/- {err_mbf}')
        forces_dict[subfile] = mean_bond_force
        error_dict[subfile] = err_mbf

    print('writing forces to file: ', f'{potential}/bond_forces.json')
    with open(f'{potential}/bond_forces.json', 'w') as f:
        json.dump(forces_dict, f)
    print('writing errors to file: ', f'{potential}/bond_force_errors.json')
    with open(f'{potential}/bond_force_errors.json', 'w') as f:
        json.dump(error_dict, f)

