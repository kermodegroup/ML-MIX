import numpy as np

cheap_pot = '../HeW'

cmds = [f'pair_style uf3 3',
            f'pair_coeff * * {cheap_pot}.uf3 W He']

nprocs=0
mass_W = 183.84
mass_He = 4.0026
in_file_root = 'W_input_16'
in_file_path = '.'
base_struct_name = '../W_input_16.xyz'
multi_potential = False
dump_name = 'dump.lammpstrj'
dump_freq = 100
n_steps = 160*nprocs
mlml_nevery = 1
rqm = 6.0
bw = 6.0
rblend = 4.0
damping = 1.0
dump_files = False
thermo_freq = 1
HeW_sim = True
rseed = 12345