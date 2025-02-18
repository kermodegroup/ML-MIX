from ase.calculators.lammpslib import LAMMPSlib
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import utils
import numpy as np
import ase
from ase.io.lammpsdata import write_lammps_data
from matscipy.neighbours import neighbour_list
mass = 28.0855
pot = '0ce218cb3e'

mass_cmds = [f'mass 1 {mass}']
cmds = ['pair_style hybrid/overlay pace table spline 5401',
                        f'pair_coeff * * pace {pot}.yace Si', 
                        f'pair_coeff 1 1 table {pot}_pairpot.table Si_Si']

calc = LAMMPSlib(lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)
si_vac, r_vac, vac_pos = utils.make_vacancy(calc, n_cell=10, rattle=False)


#get the atom that's closest to the vacancy
near_index = np.argmin(r_vac)

#get vector between nearest atom and vacancy
r_vac_near = vac_pos - si_vac.positions[near_index]

#shift the atom to the vacancy
final_image = si_vac.copy()
final_image.positions[near_index] += r_vac_near

structs = [si_vac, final_image]

mask = np.zeros(len(si_vac), dtype=bool)
# get all atoms with low coordination
for struct in structs:
    i = neighbour_list('i', struct, 4.0)
    coord = np.bincount(i)
    mask = mask | (coord == np.min(coord))

si_vac.arrays["selected_atoms"] = mask

ase.io.write('first_neb_image.xyz',si_vac)
write_lammps_data('first_neb_image.lj',si_vac,atom_style='atomic')
ase.io.write('final_neb_image.xyz',final_image)

# for the final image, just write the number of atoms, and on each line the atom id and coordinates
with open('final_neb_image.neb_file','w') as f:
    f.write(f'{len(final_image)}\n')
    for i, pos in enumerate(final_image.positions):
        f.write(f'{i+1} {pos[0]} {pos[1]} {pos[2]}\n')



