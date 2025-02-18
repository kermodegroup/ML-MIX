import sys
sys.path.append('.')
import params
from matscipy import parameter
from ase.optimize import LBFGS
from ase.constraints import ExpCellFilter
from ase.lattice.cubic import Diamond, BodyCenteredCubic
import ase.io
import numpy as np

calc = parameter("calc")
ns = parameter("ns")
in_file_name = parameter("in_file_name")
el = 'W'
a0_init = 3.16
lattice = BodyCenteredCubic
print('optimising lattice parameter')
unit_cell = lattice(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
unit_cell.calc = calc
ecf = ExpCellFilter(unit_cell)
uc_optimise = LBFGS(ecf)
uc_optimise.run(fmax=0.0001)
a0 = unit_cell.get_cell().diagonal()[0]

for n in ns:
    struct = unit_cell.copy()
    struct *= (n,n,n)
    #pick an atom near the centre
    cell = struct.get_cell()
    pos = struct.get_positions()
    midpoint = cell.diagonal()/2
    #get position of atom nearest to the centre
    dists = np.linalg.norm(pos-midpoint,axis=1)
    atom = np.argmin(dists)
    position = pos[atom][0]
    He = ase.Atoms('He', positions=np.array([[position + a0/2,position,position]]))
    struct += He
    struct.arrays["selected_atoms"] = np.zeros(len(struct))
    struct.arrays["selected_atoms"][-1] = 1
    print(f'writing {in_file_name}{n}.xyz, natoms: {len(struct)}')
    ase.io.write(f'{in_file_name}{n}.xyz',struct)