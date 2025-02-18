import numpy as np

expensive_pot = '../0ce218cb3e'

cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {expensive_pot}.yace Si', 
        f'pair_coeff 1 1 table {expensive_pot}_pairpot.table Si_Si']

nprocs=0
mass = 28.0855
in_file_root = 'Si_input_32'
in_file_path = '.'
base_struct_name = '../Si_input_32.xyz'
multi_potential = False
dump_name = 'dump.lammpstrj'
dump_freq = 100
n_steps = 5*nprocs
mlml_nevery = 1
rqm = 6.0
bw = 6.0
rblend = 0.0
damping = 1.0
dump_files = False
thermo_freq = 1


rseed = 12345