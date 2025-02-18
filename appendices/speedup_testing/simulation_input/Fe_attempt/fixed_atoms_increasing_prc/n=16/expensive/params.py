import numpy as np

expensive_pot = '../jace_espresso'

cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {expensive_pot}.yace Fe', 
        f'pair_coeff 1 1 table {expensive_pot}_pairpot.table Fe_Fe']

nprocs=0
mass = 55.845
in_file_root = 'Fe_input_16'
in_file_path = '.'
base_struct_name = '../Fe_input_16.xyz'
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


rseed = 12345
