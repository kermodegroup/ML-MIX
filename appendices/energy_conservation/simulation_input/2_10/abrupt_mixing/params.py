expensive_pot = '/home/eng/phrffk/mix_potentials/MM-MM-retrain/case_study_1_Si_stretched/0ce218cb3e'
cheap_pot = '/home/eng/phrffk/mix_potentials/MM-MM-retrain/case_study_1_Si_stretched/2_10/model_constrained_alpha_1e-7_2_10'



cmds = ['pair_style hybrid/overlay/mlml pace pace table spline 5000 table spline 5000',
        f'pair_coeff * * pace 1 1 {expensive_pot}.yace Si', 
        f'pair_coeff 1 1 table 1 1 {expensive_pot}_pairpot.table Si_Si',
        f'pair_coeff * * pace 2 2 {cheap_pot}.yace Si', 
        f'pair_coeff 1 1 table 2 2 {cheap_pot}_pairpot.table Si_Si']



value=0
mass = 28.0855
in_file_root = '300K_bond_thermalised'
in_file_path = '/home/eng/phrffk/mix_potentials/MM-MM-retrain/case_study_1_Si_stretched/'
base_struct_name = '/home/eng/phrffk/mix_potentials/MM-MM-retrain/case_study_1_Si_stretched/sb_struct.xyz'
multi_potential = True
dump_name = 'dump.lammpstrj'
bond_dump_name = 'stretched_bond_dump.lammpstrj'
T = 300
dump_freq = 1000
n_steps = 50000
nve_sim = True
mlml_nevery = 1
rqm = value
bw = 6.0
rblend = 0.0
damping = 1.0
dump_files = True
thermo_freq = 100


rseed = 12345