import os

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

