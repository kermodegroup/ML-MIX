import numpy as np

cheap_pot = '../model_constrained_alpha_1e-7_2_10'

cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {cheap_pot}.yace Si', 
        f'pair_coeff 1 1 table {cheap_pot}_pairpot.table Si_Si']

n=0
natoms=(n**3)*8
print("natoms:", natoms)
n_steps = int((540)*(64000/natoms)) # aim for about 40 min of runtime
print("n_steps",n_steps)
mass = 28.0855
in_file_root = f'Si_input_{n}'
in_file_path = '.'
base_struct_name = f'../../Si_input_{n}.xyz'
multi_potential = False
dump_name = 'dump.lammpstrj'
dump_freq = 100
mlml_nevery = 1
rqm = 6.0
bw = 6.0
rblend = 0.0
damping = 1.0
dump_files = False
thermo_freq = 1


rseed = 12345
