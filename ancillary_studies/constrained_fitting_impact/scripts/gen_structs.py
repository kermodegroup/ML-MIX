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

import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
sys.path.append('.')
import ase.io
import numpy as np
from ase.io.lammpsdata import write_lammps_data
from matscipy import parameter
import params


unit_cell = parameter("unit_cell")
n = parameter("n")
strains = parameter("strains")
rqm = parameter("rqm")

full_cell = unit_cell*(n,n,n)
#get id of atom closest to cell centre
cell_middle = full_cell.get_cell().diagonal()/2
atom_positions = full_cell.get_positions()
atom_distances = np.linalg.norm(atom_positions - cell_middle,axis=1)
closest_atom = np.argmin(atom_distances)
print(f'Atom closest to cell centre: {closest_atom}')
#get all atoms within 6.0 angstroms of the closest atom
close_atoms = np.where(atom_distances < 6.0)[0]


# apply strain as simple xy shear
strain_matrix = np.eye(3)
for i,strain in enumerate(strains):
    strain_matrix[0,1] = strain

    # generate the strained cell
    cell = full_cell.get_cell()
    cell = np.dot(cell,strain_matrix)
    full_cell.set_cell(cell, scale_atoms=True)        
    full_cell.arrays['selected_atoms'] = np.zeros(len(full_cell))
    full_cell.arrays['selected_atoms'][close_atoms] = 1
    ase.io.write(f'shear_{strain}.xyz',full_cell)
    write_lammps_data(f'shear_{strain}.lj',full_cell, velocities=True)


