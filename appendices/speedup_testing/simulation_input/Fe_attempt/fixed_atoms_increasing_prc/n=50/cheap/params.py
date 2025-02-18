import numpy as np

cheap_pot = '../model_constrained_alpha_1e-7_2_10_espresso'

cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {cheap_pot}.yace Fe', 
        f'pair_coeff 1 1 table {cheap_pot}_pairpot.table Fe_Fe']

nprocs=0
mass = 55.845
in_file_root = 'Fe_input_50'
in_file_path = '.'
base_struct_name = '../Fe_input_50.xyz'
multi_potential = False
dump_name = 'dump.lammpstrj'
dump_freq = 100
n_steps = 6*nprocs
mlml_nevery = 1
rqm = 6.0
bw = 6.0
rblend = 4.0
damping = 1.0
dump_files = False
thermo_freq = 1


rseed = 12345