import numpy as np
import ase.io
import mpi4py as MPI


def read_lammps_dump(filename):
    traj = ase.io.read(filename,parallel=False)
    print(traj)
    snapshots = []
    current_snapshot = []
    reading_atoms = False
    headers = []
    
    with open(filename, 'r') as file:
        for line in file:
            if line.startswith("ITEM: TIMESTEP"):
                if current_snapshot:
                    current_snapshot = np.array(current_snapshot)
                    snapshots.append(current_snapshot[np.argsort(current_snapshot[:,0])])
                    current_snapshot = []
                next(file)  # Skip the timestep value
            elif line.startswith("ITEM: ATOMS"):
                headers = line.split()[2:]
                reading_atoms = True
            elif reading_atoms:
                values = line.split()
                if "i2_potential[1]" in headers:
                    id_idx = headers.index("id")
                    i2_idx_1 = headers.index("i2_potential[1]")
                    i2_idx_2 = headers.index("i2_potential[2]")
                    d2_idx_1 = headers.index("d2_eval[1]")
                    d2_idx_2 = headers.index("d2_eval[2]")
                    current_snapshot.append([int(values[id_idx]),
                                             int(values[i2_idx_1]), 
                                             int(values[i2_idx_2]), 
                                             float(values[d2_idx_1]), 
                                             float(values[d2_idx_2])])
    
    if current_snapshot:
        current_snapshot = np.array(current_snapshot)
        snapshots.append(current_snapshot[np.argsort(current_snapshot[:,0])])
    traj.arrays['i2_potential[1]'] = snapshots[0][:,1]
    traj.arrays['i2_potential[2]'] = snapshots[0][:,2]
    traj.arrays['d2_eval[1]'] = snapshots[0][:,3]
    traj.arrays['d2_eval[2]'] = snapshots[0][:,4]
    return traj

def get_seed_atoms(struct):
    #find 2 atoms nearest center of struct
    center = np.mean(struct.get_positions(), axis=0)
    dists = np.linalg.norm(struct.get_positions() - center, axis=1)
    idx = np.argsort(dists)
    seed_atoms = idx[:2]
    return seed_atoms

def build_regions_lammps(lmps, struct, r_core, r_blend, r_buff, comm=None, rank=0):
    largest_cutoff = np.max((r_core,r_buff,r_blend))
    lmps.command(f'comm_modify cutoff {largest_cutoff+2.0}')
    seed_atoms = get_seed_atoms(struct)
    # set up dump
    lmps.command('dump dump1 all custom 1 dump.lammpstrj id type x y z fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]')
    # set up fix
    lmps.command(f'group seed_atoms id {" ".join([str(i+1) for i in seed_atoms])}')
    lmps.command(f'fix mlml_fix all mlml 1 {r_core} {r_buff} {r_blend} group seed_atoms')

    lmps.command('run 0')
    
    if rank == 0:
        out_dump = read_lammps_dump('dump.lammpstrj')

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


def build_regions_python(struct, r_core, r_blend, r_buff, comm=None, rank=0):
    seed_atoms = get_seed_atoms(struct)
    pos = struct.get_positions()

    i2_potential = np.zeros((len(struct),2), dtype=bool)
    d2_eval = np.zeros((len(struct),2), dtype=float)
    core = np.zeros(len(struct), dtype=bool)

    i2_potential[:,1] = True
    for idx in seed_atoms:
        pos_atom = pos[idx,:]
        dists = np.linalg.norm(pos - pos_atom, axis=1)
        core = (dists<r_core)|core

    core_atoms = np.where(core)[0]
    i2_potential[core,0] = True
    i2_potential[core,1] = False
    blend = np.zeros(len(struct), dtype=bool)

    d2_eval[core,0] = 1.0

    for idx in core_atoms:
        pos_atom = pos[idx,:]
        dists = np.linalg.norm(pos - pos_atom, axis=1)
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
        pos_atom = pos[idx,:]
        dists = np.linalg.norm(pos - pos_atom, axis=1)
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
        ase.io.write('struct_i2.xyz', struct, parallel=False)
    


    return i2_potential, d2_eval



