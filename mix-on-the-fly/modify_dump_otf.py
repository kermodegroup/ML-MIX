print("########Entering Python########")
import sys
import os
import yaml
import ase.io
import warnings
import numpy as np
from ase.calculators.singlepoint import SinglePointCalculator
print("Python path: ", sys.executable)
print("Generating otf-fit-data.xyz...")

config_file = sys.argv[1]
config = yaml.safe_load(open(config_file, 'r'))
# check if essential keys in config file
if 'force-keys' not in config:
    raise ValueError("force-keys not in config file")
if 'fit-elements' not in config:
    raise ValueError("fit-elements not in config file")
if 'dump-name' not in config:
    raise ValueError("dump-name not in config file")
if 'mask-key' not in config:
    config['mask-key'] = None
if 'pae-key' not in config:
    config['pae-key'] = None
if 'auto-mask' not in config:
    config['auto-mask'] = None


traj = ase.io.read(config['dump-name'], index=':')

#if traj is not a list, make it a list
if not isinstance(traj, list):
    traj = [traj]

fit_elements = config['fit-elements']

new_traj = []
for i, t_it in enumerate(traj):
    t = t_it.copy()
    t.calc = t_it.calc
    keep_arr_keys = ['positions', 'numbers', 'momenta']

    symbols = t.get_chemical_symbols()
    if i == 0:
        # check the first traj to see if all the symbols are H
        if all(symbol == 'H' for symbol in t.get_chemical_symbols()):
            warnings.warn("All atoms are H, did you forget to dump the mass?")

    if config['auto-mask'] is None and config['mask-key'] is None:
        warnings.warn("No criteria for masking provided, fitting all atoms in fit-elements")
    
    mask = np.ones(len(t), dtype=bool)
    if config['mask-key'] is not None:   
        mask = t.arrays[config['mask-key']]>0.99
        t.arrays['mask'] = mask
        keep_arr_keys.append('mask')
    
    if config['auto-mask'] is not None:
        # get all unique symbols
        unique_symbols = np.unique(symbols)
        for symb in unique_symbols:
            if symb not in fit_elements:
                atom_ids = np.where(symbols == symb)[0]
                pos = t.get_positions()
                for atom_id in atom_ids:
                    distances = np.linalg.norm(pos - pos[atom_id], axis=1)
                    mask[distances < config['auto_mask']['auto_mask_cutoff']] = False

    # if the mask is all False, cycle to next loop (do not add)
    if np.all(~mask):
        continue
    default_keys = False
    for key in ['fx', 'fy', 'fz']:
        if config['force-keys'][key] == key:
            default_keys = True
        else:
            if default_keys:
                raise ValueError("If you change one force key from fx, fy, fz, you have to change all of them")
    

    if default_keys:
        try:
            forces = t.get_forces()
        except RuntimeError:
            raise ValueError("No forces with keys fx fy fz in trajectory")
        t.arrays['fit_forces'] = forces
        
    if not default_keys:
        forces = np.zeros((len(t), 3))
        for j, key in enumerate(['fx', 'fy', 'fz']):
            forces[:, j] = t.arrays[config['force-keys'][key]].flatten()

        t.arrays['fit_forces'] = forces
    
    keep_arr_keys.append('fit_forces')
    
    if config['pae-key'] is not None:
        t.arrays['pae'] = t.arrays[config['pae-key']]
        keep_arr_keys.append('pae')
    
    # now, remove all the keys we don't want to keep
    keys_list = list(t.arrays.keys())
    for key in keys_list:
        if key not in keep_arr_keys:
            del t.arrays[key]
    
    # now, delete all atoms with symbols not in fit_elements
    el_mask = np.array([symbol in fit_elements for symbol in symbols])
    t = t[el_mask]
    new_traj.append(t)

ase.io.write("otf-fit-data.xyz", new_traj, format='extxyz', write_results=True)
print("########Exiting Python########")
