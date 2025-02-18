import numpy as np

cheap_pot = '../HeW'

cmds = [f'pair_style uf3 3',
            f'pair_coeff * * {cheap_pot}.uf3 W He']
n=0
natoms=(n**3)*2
print("natoms:", natoms)
n_steps = int((360)*(64000/natoms)) # aim for about 40 min of runtime
print("n_steps",n_steps)
mass_W = 183.84
mass_He = 4.0026
in_file_root = f'W_input_{n}'
in_file_path = '.'
base_struct_name = f'../../W_input_{n}.xyz'
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

HeW_sim = True

rseed = 12345
