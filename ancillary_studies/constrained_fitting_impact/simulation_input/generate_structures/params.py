from ase.calculators.lammpslib import LAMMPSlib
from ase.lattice.cubic import Diamond
from ase.constraints import ExpCellFilter
from ase.optimize import LBFGS
import numpy as np
n = 10
strains = np.arange(0,0.03,0.005)

mass = 28.0855 # atomic mass of Si
pot = '../0ce218cb3e'
#generate a large vacancy cell
mass_cmds = [f'mass 1 {mass}']
cmds = ['pair_style hybrid/overlay pace table spline 5401',
                        f'pair_coeff * * pace {pot}.yace Si', 
                        f'pair_coeff 1 1 table {pot}_pairpot.table Si_Si']

calc = LAMMPSlib(lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)
r_init = 8
el              = 'Si'
a0_init         = 5.43
lattice = Diamond
print('optimising lattice parameter')
unit_cell = Diamond(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
unit_cell.calc = calc
ecf = ExpCellFilter(unit_cell)
uc_optimise = LBFGS(ecf)
uc_optimise.run(fmax=0.0001)
a0 = unit_cell.get_cell()[0,0] #get the optimised lattice parameter

rqm = 6.0