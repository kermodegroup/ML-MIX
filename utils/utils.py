import os
import numpy as np
from ase.build import bulk
from ase.lattice.cubic import Diamond
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS
import os
import ase.io

def bare_bones_setup_lammps(lmps,input_file,mass_cmds,calc_commands,
                             sim_tstep=0.001, thermo_freq=100,multi_potential=False):
    
    """Set up the simulation by passing an active LAMMPS object a number of commands"""
    # ---------- Initialize Simulation --------------------- 

    utils_path = os.path.dirname(os.path.abspath(__file__))

    lmps.command('clear')
    lmps.command(f'plugin load {utils_path}/../LAMMPS_plugin/build/mlmlplugin.so')
    lmps.command(f'plugin load {utils_path}/../LAMMPS_plugin/build/langevinmlmlplugin.so')
    lmps.command(f'plugin load {utils_path}/../LAMMPS_plugin/build/hybridoverlaymlmlplugin.so')

    lmps.command('dimension 3')
    lmps.command('boundary p p p')
    lmps.command('atom_style atomic')
    lmps.command('units metal')
    #lmps.command('comm_style tiled')

    #-----------------Read atoms------------------
    if multi_potential:
        lmps.command('fix eval_pot all property/atom i2_potential 2 ghost yes')
        lmps.command('fix eval_arr all property/atom d2_eval 2 ghost yes')
    lmps.command(f'read_data {input_file}')

    #----------Define potential-------
    lmps.commands_list(mass_cmds)

    ############ set up potential ################
    #set up both potentials, ensuring that interactions between types are defined
    lmps.commands_list(calc_commands)
    lmps.command('fix sf all store/force')

    # ---------- define potential energy compute ---------
    lmps.command('compute pe_peratom all pe/atom')

    # ---------- set timestep length -------------
    lmps.command(f'timestep {sim_tstep}')

    # Specify the output frequency for thermo data
    lmps.command(f'thermo {thermo_freq}')

def make_vacancy(calc, n_cell=10, rattle=True):
    """Create a supercell of Si with a vacancy in the centre."""
    el              = 'Si'
    a0_init         = 5.43
    lattice = Diamond
    print('optimising lattice parameter')
    unit_cell = Diamond(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
    unit_cell.calc = calc
    ecf = ExpCellFilter(unit_cell)
    uc_optimise = LBFGS(ecf)
    uc_optimise.run(fmax=0.0001)
    a0 = unit_cell.get_cell()[0,0] #get the optimised lattice parameter

    si = bulk('Si', a=a0, cubic=True) * n_cell

    # pick atom closest to centre of cell
    vac_index = ((si.positions - np.diag(si.cell)/2) ** 2).sum(axis=1).argmin()
    vac_pos = si.positions[vac_index]

    # measure distance from vacancy to all other atoms
    r_vac = si.get_distances(vac_index, np.arange(0, len(si)), mic=True)
    
    del si[vac_index]
    r_vac = np.r_[r_vac[:vac_index], r_vac[vac_index+1:]]
    
    if rattle:
        si.rattle(0.05) # displace atoms from equilibrium to break symmetry
    
    print(f'Si vacancy cell with {len(si)} atoms')
    return si, r_vac, vac_pos


def read_lammps_dump(filename):
    traj = ase.io.read(filename,parallel=False, index=':')
    if len(traj) == 1:
        traj = [traj]
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
                reading_atoms=False
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

    print(f'Read {len(snapshots)} snapshots')
    print(snapshots)
    for i, t in enumerate(traj):
        t.arrays['i2_potential[1]'] = snapshots[i][:,1]
        t.arrays['i2_potential[2]'] = snapshots[i][:,2]
        t.arrays['d2_eval[1]'] = snapshots[i][:,3]
        t.arrays['d2_eval[2]'] = snapshots[i][:,4]
    return traj