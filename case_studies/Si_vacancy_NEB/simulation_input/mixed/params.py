expensive_pot = '../0ce218cb3e'
cheap_pot = '../model_constrained_alpha_1e-7_2_10'

nimages = 9

cmds = ['pair_style hybrid/overlay/mlml pace pace table spline 5000 table spline 5000',
        f'pair_coeff * * pace 1 1 {expensive_pot}.yace Si', 
        f'pair_coeff 1 1 table 1 1 {expensive_pot}_pairpot.table Si_Si',
        f'pair_coeff * * pace 2 2 {cheap_pot}.yace Si', 
        f'pair_coeff 1 1 table 2 2 {cheap_pot}_pairpot.table Si_Si']

mass = 28.0855
in_file_name = "first_neb_image"
in_file_path = '../'
final_image_name = "final_neb_image"
multi_potential = True
dump_name = 'dump.lammpstrj'
mlml_nevery = 0
rqm = 4.0
bw = 6.0
rblend = 4.0
damping = 2.0
thermo_freq = 1
rseed = 12345
