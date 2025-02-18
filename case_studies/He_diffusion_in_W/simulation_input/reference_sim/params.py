import numpy as np
#mass of W
mass_W = 183.84
mass_He = 4.0026
pot = '../W_He_ace'
T = 400
n = 16

i = 0
he_dump_name=f'HeDump{i}.lammpstrj'
in_file_name=f'../W_He_input_{n}_unrelaxed'
multi_potential=False
dump_name = 'dump.lammpstrj'
dump_freq = 30
he_dump_freq = 10
dump_files = False
n_steps = 60000
therm_time = 2000
damping = 0.1
relax = False

tstep=0.001

mlml_nevery=1
rqm=6.0
bw=6.0
rblend=4.0
qm_type=2
weak_damping = 2.0


rseed = np.random.randint(0,100000)

ace_cmds = ['pair_style hybrid/overlay pace table spline 5000',
            f'pair_coeff * * pace {pot}.yace W He', 
            f'pair_coeff 1 1 table {pot}_pairpot.table W_W',
            f'pair_coeff 2 1 table {pot}_pairpot.table He_W',
            f'pair_coeff 2 2 table {pot}_pairpot.table He_He']

