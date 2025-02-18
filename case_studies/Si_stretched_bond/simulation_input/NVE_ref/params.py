expensive_pot = '../0ce218cb3e'



cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {expensive_pot}.yace Si', 
        f'pair_coeff 1 1 table {expensive_pot}_pairpot.table Si_Si']



value=0
mass = 28.0855
in_file_root = '300K_bond_thermalised'
in_file_path = '../'
base_struct_name = '../sb_struct.xyz'
multi_potential = False
dump_name = 'dump.lammpstrj'
bond_dump_name = 'stretched_bond_dump.lammpstrj'
T = 300
dump_freq = 100
n_steps = 50000
nve_sim = True
mlml_nevery = 1
rqm = value
bw = 6.0
rblend = 0.0
damping = 1.0
dump_files = False
thermo_freq = 100


rseed = 12345