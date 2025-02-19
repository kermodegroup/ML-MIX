import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
sys.path.append('.')
import ase.io
import numpy as np
import matplotlib.pyplot as plt
from matscipy.neighbours import neighbour_list
from ase.geometry import get_distances

import global_plot_settings
global_plot_settings.large_text()

strains = np.arange(0,0.031,0.005)
in_file_names = [f'shear_{strain}' for strain in strains]
in_file_paths = 'results/'
base_folder = ['all_expensive']
folders = ['mixed_matched','mixed_unmatched']

def strain_error(at0, u_ref, u, cutoff=3.0):
    I, J,  = neighbour_list('ij', at0, cutoff)
    v = u - u_ref
    dv = np.linalg.norm(v[I, :] - v[J, :], axis=1)
    return np.linalg.norm(dv)

def get_correct_positions(struct,path):
    pos = struct.get_positions()
    ids,ix,iy,iz = read_in_image_from_dump(f'{path}')
    #get cell vectors x y and z
    cell = struct.get_cell()
    a = cell[0]
    b = cell[1]
    c = cell[2]
    pos[ids] += ix[:,None]*a + iy[:,None]*b + iz[:,None]*c
    return pos

def read_in_image_from_dump(dump_file):
    # just read in id,ix,iy,iz columns
    # read in lammps dump file
    with open(dump_file) as f:
        lines = f.readlines()
    # find the line where the data starts
    for i, line in enumerate(lines):
        if line.startswith('ITEM: ATOMS'):
            start = i
            break
    # read in the data
    data = np.loadtxt(lines[start+1:])
    ids = np.array((data[:,0] - 1), dtype=int)
    ix = data[:,-3]
    iy = data[:,-2]
    iz = data[:,-1]
    return ids,ix,iy,iz



Dus = {}
for folder in folders:
    Dus[folder] = []

for name in in_file_names:
    base_path = f'{in_file_paths}/{base_folder[0]}/{name}.lammpstrj'
    struct = ase.io.read(base_path,parallel=False)
    base_pos = get_correct_positions(struct,base_path)
    for folder in folders:
        new_path = f'{in_file_paths}/{folder}/{name}.lammpstrj'
        try:
            struct2 = ase.io.read(new_path,parallel=False)
        except:
            print(f'No file found for {new_path}')
            Dus[folder].append(-100)
            continue
        pos = get_correct_positions(struct2,new_path)
        diff = pos - base_pos
        struct.arrays[folder] = diff
        Du = strain_error(struct,base_pos,pos)
        Dus[folder].append(Du)
    ase.io.write(f'plots/{name}.xyz',struct,format = 'extxyz')

#now plot the total diff against the strain 
# prune the -100 values from strain list

mask_matched = np.array(Dus['mixed_matched']) != -100
mask_unmatched = np.array(Dus['mixed_unmatched']) != -100
print(strains)
print(strains[mask_matched])
print(strains[mask_unmatched])

plt.plot(strains[mask_matched], np.array(Dus['mixed_matched'])[mask_matched], label='Constrained')
plt.plot(strains[mask_unmatched], np.array(Dus['mixed_unmatched'])[mask_unmatched], label='Unconstrained')
plt.xlabel(r'$\epsilon_{xy}$')
plt.ylabel('Du')
plt.legend()
plt.savefig('plots/Du vs strain.png')
