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

import numpy as np
from ase.build import bulk
from ase.calculators.lammpslib import LAMMPSlib
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS
import ase.io
from ase.calculators.singlepoint import SinglePointCalculator
import sys
sys.path.append('.')
from matscipy import parameter
import params
import os
from ase.lattice.cubic import Diamond


def create_strained_configurations(alat=3.1855, symbol="W", max_strain=0.01, steps=17, size=5, ):

    # strain_mat: uniaxial when a=b, pure shear when a.b = 0
    def strain_mat(a, b): 
        return (np.outer(a, b)+np.outer(b, a)) / np.linalg.norm(a) / np.linalg.norm(b) / 2.0
    
    strain_types = {r"Hydrostatic" : np.diag([1, 1, 1]) / 3.,
                r"Uniaxial [100]" : strain_mat([1, 0, 0] , [1, 0, 0]),
                r"Uniaxial [111]" : strain_mat([1, 1, 1] , [1, 1, 1]),
                r"Shear $[100],[010]$" : strain_mat([1, 0, 0] , [0, 1, 0]), 
                r"Shear $[110],[001]$" : strain_mat([1, 1, 0], [0, 0, 1]),
                r"Shear $[110],[1\bar{1}2]$" : strain_mat([1, 1, 0], [1, -1, 2])}

    bulk_atoms = bulk(symbol, a=alat, cubic=True)
    bulk_atoms *= size
    bulk_cell = bulk_atoms.cell.copy()[:]
    
    strained_configs = {}
    for strain_type_label, strain_type in strain_types.items():
        strained_configs[strain_type_label] = []
        for strain_value in np.linspace(-max_strain, +max_strain, steps):
            
            strain = np.eye(3) + strain_value * strain_type
            #print(strain)
            strained_atoms = bulk_atoms.copy()
            strained_atoms.set_cell(strain @ bulk_cell, scale_atoms=True)
            if lattice == Diamond:
                #relax
                strained_atoms.calc = calc
                opt = LBFGS(strained_atoms)
                opt.run(fmax=0.0001)
            
            strained_atoms.info["strain_type"] = strain_type_label
            strained_atoms.info["strain_value"] = strain_value
            
            strained_configs[strain_type_label].append(strained_atoms)
    
    return strained_configs



el = parameter('el')
a0_init = parameter('a0_init')
mass = parameter('mass')
cmds = parameter('cmds')
lattice = parameter('lattice')
Ah_n = parameter('Ah_n', 5)
ncfgs = parameter('ncfgs', 17)
max_strain = parameter('max_strain', 0.005)

#make folder for the results if it doesn't exist
if not os.path.exists('Ah_data'):
    os.makedirs('Ah_data')


calc = LAMMPSlib(amendments=[f"mass 1 {mass}"],lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)

unit_cell = lattice(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
unit_cell.calc = calc
ecf = ExpCellFilter(unit_cell)
uc_optimise = LBFGS(ecf)
uc_optimise.run(fmax=0.0001)
a0 = unit_cell.get_cell()[0,0]

config_set = create_strained_configurations(alat=a0, symbol=el, max_strain=max_strain, steps=ncfgs, size=Ah_n)
eval_configs = []
total_configs = sum([len(config_set[strain_val]) for strain_val in config_set])
i = 0
print("Evaluating configurations...")
for strain_val in config_set:
    for config in config_set[strain_val]:
        i += 1
        print(f"{i}/{total_configs}", end='\r')
        config.calc = calc
        energy = config.get_potential_energy()
        forces = config.get_forces()
        config.calc = SinglePointCalculator(config, energy=energy, forces=forces)
        eval_configs.append(config)


ase.io.write(f'Ah_data/data_for_Ah.xyz', eval_configs, parallel=False)

# Also generate pure bulk reference
bulk_atoms = unit_cell*(Ah_n,Ah_n,Ah_n)
bulk_atoms.calc = calc
energy = bulk_atoms.get_potential_energy()
forces = bulk_atoms.get_forces()
bulk_atoms.calc = SinglePointCalculator(bulk_atoms,energy=energy,forces=forces)

ase.io.write(f'Ah_data/just_bulk.xyz',bulk_atoms,parallel=False)
