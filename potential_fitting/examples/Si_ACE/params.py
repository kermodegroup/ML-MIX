from ase.lattice.cubic import Diamond, BodyCenteredCubic
lattice = Diamond
el = 'Si'
a0_init = 5.401
mass = 28.085
pot = './expensive_potential/0ce218cb3e'
cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {pot}.yace {el}', f'pair_coeff 1 1 table {pot}_pairpot.table {el}_{el}']

## Ah params
Ah_n = 3
ncfgs = 17
max_strain = 0.005

## As params
nsteps = 10000
dump_freq = 100
T = 500
damping = 0.1
rseed = 12345
As_n = 8
As_prune = 5 # take 1 in every As_prune configurations from the trajectory

## analysis params
cheap_model_type = 'ace'
