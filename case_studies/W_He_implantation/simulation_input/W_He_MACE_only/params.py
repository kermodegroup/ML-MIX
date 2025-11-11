from ase.calculators.lammpslib import LAMMPSlib
from mpi4py import MPI
from ase.lattice.cubic import BodyCenteredCubic
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS
from mace.calculators import MACECalculator

# read in from command line
# MPI setup
comm_world = MPI.COMM_WORLD
rank = comm_world.Get_rank()
#isolate each rank on a seperate communicator before passing in
single_comm = comm_world.Split(color=rank, key=rank)
el = 'W'
model_file = '../../mpa_finetune_W_He'
mass_W = 183.84
a0_init = 3.165

directions = [[1,0,0], [0,1,0], [0,0,1]]
nx = 15
ny = 15
nz = 15
vacuum = 50
He_energy=0

calc = MACECalculator(f'{model_file}.model',device="cuda",head="new")

print('optimising lattice parameter')
unit_cell = BodyCenteredCubic(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
unit_cell.set_calculator(calc)
ecf = ExpCellFilter(unit_cell)
uc_optimise = LBFGS(ecf)
uc_optimise.run(fmax=0.0001)
a0 = unit_cell.get_cell()[0,0] #get the optimised lattice parameter
print('optimised lattice parameter:',a0)

lmps_cmds = [f'pair_style mliap unified {model_file}.model-mliap_lammps.pt 0',
            f'pair_coeff * * W He']


nrepeats = 25
T = 1000
dump_files = False
dump_freq = 10
second_stage_tstep = 0.0001
second_stage_nsteps = 10000
equilibration_time = 10000
He_place_dist = 5.0 #potential cutoff

equilibration_dump = False
check_nsteps=200
gpu=True
