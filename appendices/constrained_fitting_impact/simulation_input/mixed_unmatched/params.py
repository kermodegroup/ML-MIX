import numpy as np
expensive_pot = '../0ce218cb3e'
cheap_pot = '../model_unconstrained_alpha_1e-7_2_10'
strains = np.arange(0.0,0.03,0.005)

cmds = ['pair_style hybrid/overlay/mlml pace pace table spline 5000 table spline 5000',
        f'pair_coeff * * pace 1 1 {expensive_pot}.yace Si', 
        f'pair_coeff 1 1 table 1 1 {expensive_pot}_pairpot.table Si_Si',
        f'pair_coeff * * pace 2 2 {cheap_pot}.yace Si', 
        f'pair_coeff 1 1 table 2 2 {cheap_pot}_pairpot.table Si_Si']

# do just expensive
# cmds = ['pair_style hybrid/overlay pace table spline 5000',
#         f'pair_coeff * * pace {expensive_pot}.yace Si', 
#         f'pair_coeff 1 1 table {expensive_pot}_pairpot.table Si_Si']


mass = 28.0855
in_file_names = [f'shear_{strain}' for strain in strains]
in_file_path = '../generate_structures'
multi_potential = True
dump_name = 'dump.lammpstrj'
mlml_nevery = 1
rqm = 0.0
bw = 6.0
rblend = 6.0
damping = 2.0
thermo_freq = 1
rseed = 12345