from ase.lattice.cubic import Diamond, BodyCenteredCubic
lattice = BodyCenteredCubic
el = 'Fe'
a0_init = 2.856
mass = 55.845
pot = './expensive_potential/jace_espresso'
cmds = ['pair_style hybrid/overlay pace table spline 5000',
        f'pair_coeff * * pace {pot}.yace Fe', f'pair_coeff 1 1 table {pot}_pairpot.table Fe_Fe']

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
cheap_model_type = 'ace'
