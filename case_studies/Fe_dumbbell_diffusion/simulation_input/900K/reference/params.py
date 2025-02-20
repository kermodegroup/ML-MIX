import numpy as np

val=0
expensive_pot = '../jace_espresso'
cheap_pot = '../model_constrained_alpha_1e-7_2_15'


cmds_just_expensive = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {expensive_pot}.yace Fe', 
        f'pair_coeff 1 1 table {expensive_pot}_pairpot.table Fe_Fe']


cmds_mix = ['pair_style hybrid/overlay/mlml pace pace table spline 5000 table spline 5000',
        f'pair_coeff * * pace 1 1 {expensive_pot}.yace Fe', 
        f'pair_coeff 1 1 table 1 1 {expensive_pot}_pairpot.table Fe_Fe',
        f'pair_coeff * * pace 2 2 {cheap_pot}.yace Fe', 
        f'pair_coeff 1 1 table 2 2 {cheap_pot}_pairpot.table Fe_Fe']



mass = 55.845
in_file_root = 'fe_dumbell_16'
in_file_path = '../../'
base_struct_name = '../../fe_dumbell_16.xyz'
multi_potential = False
dump_name = f'dump{val}.lammpstrj'
bond_dump_name = 'stretched_bond_dump.lammpstrj'
T = 900
dump_freq = 2000
n_steps = 1000000
rqm = 6.0
bw = 6.0
rblend = 4.0
damping = 0.01
dump_files = True
thermo_freq = 20

mlml_nevery = 1

coord_nevery = 10
coord_nrepeat = 10
coord_nfreq = 100

lb = 0.4
ub = 'inf'


rseed = np.random.randint(0, 100000)

therm_time = 2000