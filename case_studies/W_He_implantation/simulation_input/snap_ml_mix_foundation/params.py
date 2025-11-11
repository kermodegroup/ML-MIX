from ase.calculators.lammpslib import LAMMPSlib
from ase.lattice.cubic import BodyCenteredCubic
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS

parallel=False

el = 'W'
mass_W = 183.84
a0_init = 3.1805

directions = [[1,0,0], [0,1,0], [0,0,1]]
nx = 15
ny = 15
nz = 30
vacuum = 50
He_energy=5
pure_W_cmds = [
"variable zblcutinner equal 4",
"variable zblcutouter equal 4.8",
"variable zblz equal 74",
"pair_style hybrid/overlay zbl ${zblcutinner} ${zblcutouter} snap",
"pair_coeff 1 1 zbl ${zblz} ${zblz}",
"pair_coeff * * snap ../../W_2940_2017_2.snapcoeff ../../W_2940_2017_2.snapparam W"]


calc = LAMMPSlib(lmpcmds=pure_W_cmds, log_file='test.log', amendments=[f"mass 1 {mass_W}"],keep_alive=True)

print('optimising lattice parameter')
unit_cell = BodyCenteredCubic(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
unit_cell.set_calculator(calc)
print(unit_cell.get_forces())
ecf = ExpCellFilter(unit_cell)
uc_optimise = LBFGS(ecf)
uc_optimise.run(fmax=0.0001)
a0 = unit_cell.get_cell()[0,0] #get the optimised lattice parameter
print('optimised lattice parameter:',a0)

rcore=6.0
rbuff=8.0
rblend=4.0

lmps_cmds= [f'variable zblcutinner equal 4',
            f'variable zblcutouter equal 4.8',
            f'variable zblz equal 74',
            f'comm_modify cutoff {rbuff+2.0}',
            'pair_style hybrid/overlay/mlml symmetrix/mace no_mpi_message_passing snap zbl ${zblcutinner} ${zblcutouter}',
            f'pair_coeff * * symmetrix/mace 1 ../../MACE_16_w_He-2-74.json W He',
            'pair_coeff 1 1 zbl 2 ${zblz} ${zblz}',
            f'pair_coeff * * snap 2 ../../W_2940_2017_2.snapcoeff ../../W_2940_2017_2.snapparam W NULL',
            f'group He_group type 2',
            f'fix mlml_fix all mlml 1 {rcore} {rbuff} {rblend} group He_group']


nrepeats = 250
T = 1000
dump_files = False
dump_freq = 10
equilibration_time = 10000
He_place_dist = 5.0 #potential cutoff

equilibration_dump = False
gpu=True

gpu_mode="symmetrix"

second_stage_time_per_loop = 0.05
second_stage_max_time = 1

multi_potential=True
