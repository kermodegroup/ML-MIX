from ase.lattice.cubic import BodyCenteredCubic

el = 'W'
a0_init = 3.165
mass = 183.84
yace = f'expensive_potential/W_He_ace'
cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {yace}.yace W', f'pair_coeff 1 1 table {yace}_pairpot.table W_W']
lattice = BodyCenteredCubic

## Ah params
Ah_n = 5
ncfgs = 17
max_strain = 0.005

## As params
nsteps = 10000
dump_freq = 100
T = 1200
damping = 0.1
rseed = 12345
As_n = 10
As_prune = 5 # take 1 in every As_prune configurations from the trajectory

## analysis params
cheap_model_type = 'uf3'

