import numpy as np

expensive_pot = '../W_He_ace'

cmds = ['pair_style hybrid/overlay pace table spline 5000',
            f'pair_coeff * * pace {expensive_pot}.yace W He', 
            f'pair_coeff 1 1 table {expensive_pot}_pairpot.table W_W',
            f'pair_coeff 2 1 table {expensive_pot}_pairpot.table He_W',
            f'pair_coeff 2 2 table {expensive_pot}_pairpot.table He_He']

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
