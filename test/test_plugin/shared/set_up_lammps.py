import ase.io

def set_up_lammps(lmps, struct, mass_cmds, multi_potential=True, rank=0):
    
    #write input file
    input_file = 'data.in'
    if rank == 0:
        ase.io.write(input_file, struct, format='lammps-data', parallel=False)

    # ---------- Initialize Simulation --------------------- 
    lmps.command('clear') 
    lmps.command('plugin load ../../../LAMMPS_plugin/build/hybridoverlaymlmlplugin.so')
    lmps.command('plugin load ../../../LAMMPS_plugin/build/mlmlplugin.so')

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
    lmps.command('thermo 1')