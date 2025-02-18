import numpy as np
import ase.io
import mpi4py as MPI
from matscipy.neighbours import neighbour_list
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../utils')
from utils import read_lammps_dump

def get_seed_atoms(struct):
    #find 2 atoms nearest center of struct
    center = np.mean(struct.get_positions(), axis=0)
    dists = np.linalg.norm(struct.get_positions() - center, axis=1)
    idx = np.argsort(dists)
    seed_atoms = idx[:2]
    return seed_atoms

def build_regions_lammps(lmps, struct, r_core, r_blend, r_buff, pick_seed_with='group', nsteps=1, fix_nevery=10, comm=None, rank=0, path='./'):
    largest_cutoff = np.max((r_core,r_buff,r_blend))
    lmps.command(f'comm_modify cutoff {largest_cutoff+2.0}')
    seed_atoms = get_seed_atoms(struct)
    # set up dump
    lmps.command(f'dump dump1 all custom {nsteps} {path}/dump.lammpstrj id type x y z fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]')
    lmps.command(f'group seed_atoms id {" ".join([str(i+1) for i in seed_atoms])}')
    if pick_seed_with == 'group':
        lmps.command(f'fix mlml_fix all mlml 1 {r_core} {r_buff} {r_blend} group seed_atoms')
    elif pick_seed_with == 'fix':
        lmps.command('compute ca all coord/atom cutoff 4.0')
        lmps.command(f'fix av_ca all ave/atom 1 1 {fix_nevery} c_ca')
        lmps.command(f'fix mlml_fix all mlml 1 {r_core} {r_buff} {r_blend} fix_classify av_ca {fix_nevery} 45.5 inf')
    elif pick_seed_with == 'fix_and_init_group':
        lmps.command('compute ca all coord/atom cutoff 4.0')
        lmps.command(f'fix av_ca all ave/atom 1 1 {fix_nevery} c_ca')
        lmps.command(f'fix mlml_fix all mlml 1 {r_core} {r_buff} {r_blend} fix_classify av_ca {fix_nevery} 45.5 inf init_group seed_atoms')

    if pick_seed_with != 'group':    
        lmps.command(f'dump fix_dump all custom {fix_nevery} {path}/fix_dump.lammpstrj id type x y z fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2] f_av_ca')
    lmps.command(f'run {nsteps}')
    
    if rank == 0:
        out_dump = read_lammps_dump(f'{path}/dump.lammpstrj')[-1]
        i2_potential_1 = out_dump.arrays['i2_potential[1]']
        i2_potential_2 = out_dump.arrays['i2_potential[2]']
        d2_eval_1 = out_dump.arrays['d2_eval[1]']
        d2_eval_2 = out_dump.arrays['d2_eval[2]']
        i2_potential = np.zeros((len(struct),2), dtype=bool)
        i2_potential[:,0] = i2_potential_1.flatten()
        i2_potential[:,1] = i2_potential_2.flatten()
        d2_eval = np.zeros((len(struct),2), dtype=float)
        d2_eval[:,0] = d2_eval_1.flatten()
        d2_eval[:,1] = d2_eval_2.flatten()
    else:
        d2_eval = None
        i2_potential = None
    
    if comm is not None:
        i2_potential = comm.bcast(i2_potential, root=0)
        d2_eval = comm.bcast(d2_eval, root=0)
    return i2_potential, d2_eval


def build_regions_python(struct, r_core, r_blend, r_buff, pick_seed_with='group', comm=None, rank=0, path='./'):
    
    i2_potential = np.zeros((len(struct),2), dtype=bool)
    d2_eval = np.zeros((len(struct),2), dtype=float)
    core = np.zeros(len(struct), dtype=bool)

    if pick_seed_with == 'group':
        seed_atoms = get_seed_atoms(struct)
    elif pick_seed_with == 'fix':
        i = neighbour_list('i', struct, 4.0)
        coord = np.bincount(i)
        seed_atoms = np.where((coord>=46))[0]
    elif pick_seed_with == 'all_expensive':
        i2_potential[:,0] = True
        i2_potential[:,1] = False
        d2_eval[:,0] = 1.0
        d2_eval[:,1] = 0.0
        core = np.ones_like(core)
        struct.arrays['i2_potential'] = i2_potential
        struct.arrays['core'] = core
        struct.arrays['d2_eval'] = d2_eval
        if rank == 0:
            ase.io.write(f'{path}/struct_i2.xyz', struct, parallel=False)
        return i2_potential, d2_eval



    i2_potential[:,1] = True
    for idx in seed_atoms:
        dists = struct.get_distances(idx, np.arange(len(struct)), mic=True)
        core = (dists<r_core)|core

    core_atoms = np.where(core)[0]
    i2_potential[core,0] = True
    i2_potential[core,1] = False
    blend = np.zeros(len(struct), dtype=bool)

    d2_eval[core,0] = 1.0

    for idx in core_atoms:
        dists = struct.get_distances(idx, np.arange(len(struct)), mic=True)
        blend_adj = ((dists<r_blend)&(~core))
        blend_adj_idx = np.where(blend_adj)[0]
        for idx_adj in blend_adj_idx:
            dist = dists[idx_adj]
            d2_eval[idx_adj,0] = max(d2_eval[idx_adj,0], 1.0 - dist/r_blend)
        blend = blend_adj|blend
    
    i2_potential[blend,0] = True
    blend_atoms = np.where(blend)[0]
    
    buff_0 = np.zeros(len(struct), dtype=bool)
    buff_1 = np.zeros(len(struct), dtype=bool)
    for idx in blend_atoms:
        dists = struct.get_distances(idx, np.arange(len(struct)), mic=True)
        buff_0 = ((dists<r_buff)&(~(core|blend)))|buff_0
        buff_1 = ((dists<r_buff)&(core))|buff_1

    i2_potential[buff_0,0] = True
    i2_potential[buff_1,1] = True
    d2_eval[:,1] = 1.0 - d2_eval[:,0]

    struct.arrays['i2_potential'] = i2_potential
    struct.arrays['core'] = core
    struct.arrays['blend'] = blend
    struct.arrays['buff_0'] = buff_0
    struct.arrays['buff_1'] = buff_1
    struct.arrays['d2_eval'] = d2_eval
    if rank == 0:
        ase.io.write(f'{path}/struct_i2.xyz', struct, parallel=False)
    


    return i2_potential, d2_eval



