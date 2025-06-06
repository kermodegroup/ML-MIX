import os

import ase.io

def set_up_lammps(lmps, struct, mass_cmds, multi_potential=True, rank=0,path='./', kokkos=False):
    
    #write input file
    input_file = f'{path}/data.in'
    if rank == 0:
        ase.io.write(input_file, struct, format='lammps-data', parallel=False)

    # ---------- Initialize Simulation --------------------- 
    lmps.command('clear') 
    plugin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../LAMMPS_plugin/build'))
    if not kokkos:
        lmps.command(f'plugin load {plugin_path}/hybridoverlaymlmlplugin.so')
        lmps.command(f'plugin load {plugin_path}/mlmlplugin.so')

    lmps.command('dimension 3')
    lmps.command('boundary p p p')
    lmps.command('atom_style atomic')
    lmps.command('units metal')
    lmps.command('atom_modify map yes')
    #lmps.command('comm_style tiled')

    #-----------------Read atoms------------------
    if multi_potential:
        if not kokkos:
            lmps.command('fix eval_pot all property/atom i2_potential 2 ghost yes')
            lmps.command('fix eval_arr all property/atom d2_eval 2 ghost yes')
        else:
            lmps.command('fix ml_mix_props all property/atom d_potential_1 d_potential_2 d_eval_1 d_eval_2 ghost yes')

    lmps.command(f'read_data {input_file}')

    #----------Define potential-------
    lmps.commands_list(mass_cmds)
    lmps.command('thermo 1')