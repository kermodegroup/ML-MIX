# ----------------------------------------------------------------------
# This script is part of the ML-MIX repository.
# 
# Copyright (2025) Fraser Birks
#
# This script is licensed under the MIT License. You are free to use, modify,
# and distribute it, provided that this copyright notice and license text 
# remain intact.
#
# See the LICENSE file in the repository root for full details.
# ----------------------------------------------------------------------

from ase.calculators.lammpslib import LAMMPSlib
from ase.build import bulk
from ase.constraints import ExpCellFilter
from ase.lattice.cubic import Diamond, FaceCenteredCubic, SimpleCubic, BodyCenteredCubic
from matscipy.elasticity import fit_elastic_constants
from ase.optimize import LBFGS
import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('.')
from matscipy import parameter
import params
import ase.units as units


#get command line arguments in the form of
# python test_elastic_constants.py "$folder_name/UF3_lammps" "$folder_name/elastic_constants" "$investigation_name" "$model_name1" "$model_name2" "$model_name3" "$model_name4"

lammps_pot_path = sys.argv[1]
output_path = sys.argv[2]
investigation_name = sys.argv[3]
model_names = sys.argv[4:]

el = parameter('el')
a0_init = parameter('a0_init')
mass = parameter('mass')
lattice = parameter('lattice')
cmds = parameter('cmds')
cheap_model_type = parameter('cheap_model_type')

calc = LAMMPSlib(amendments=[f"mass 1 {mass}"],lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)
unit_cell = lattice(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
unit_cell.calc = calc


ecf = ExpCellFilter(unit_cell)
uc_optimise = LBFGS(ecf)
uc_optimise.run(fmax=0.0001)
elastic_symmetry = 'cubic'
a_ref = unit_cell.get_cell()[0,0]
print("Reference lattice parameter:", a_ref)
unit_cell_ace = unit_cell.copy()


#fit elastic constants
C_ref, C_err_ref = fit_elastic_constants(unit_cell,symmetry=elastic_symmetry, verbose=False, optimizer=LBFGS)

C_ref = C_ref/units.GPa
C_err_ref = C_err_ref/units.GPa

print(f'Reference elastic constants: {C_ref[0,0]:.2f} {C_ref[0,1]:.2f} {C_ref[3,3]:.2f}')

def get_elastic_constants(model_name, model_type):
    if model_type == 'ace':
        cmds = ['pair_style hybrid/overlay pace table spline 5000',
            f'pair_coeff * * pace {lammps_pot_path}/{model_name}.yace {el}', 
            f'pair_coeff 1 1 table {lammps_pot_path}/{model_name}_pairpot.table {el}_{el}']
    elif model_type == 'uf3':
        cmds= ['pair_style uf3 3',
            f'pair_coeff * * {lammps_pot_path}/{model_name}.uf3 {el}']


    calc = LAMMPSlib(amendments=[f"mass 1 {mass}"], lmpcmds = cmds,log_file='lammps_output.log',keep_alive=True)
    elastic_symmetry = 'cubic'
    uc = lattice(size=[1,1,1],symbol=el,latticeconstant=a0_init,pbc=(1,1,1))
    uc.calc = calc
    #get a0
    ecf = ExpCellFilter(uc)
    uc_optimise = LBFGS(ecf)
    try:
        uc_optimise.run(fmax=0.0001, steps=1000)
        a = uc.get_cell()[0,0]
        print("Optimised lattice parameter:", a)
    except:
        print(f'Failed to optimise {model_name}')
        uc = unit_cell_ace.copy()
        a = 0.0

    #fit elastic constants
    C, C_err = fit_elastic_constants(uc,symmetry=elastic_symmetry, verbose=False, optimizer=LBFGS)
    C = C/units.GPa
    C_err = C_err/units.GPa

    print(f'{model_name} elastic constants: {C[0,0]:.2f} {C[0,1]:.2f} {C[3,3]:.2f}')

    return C, C_err, a


c11 = []
c12 = []
c44 = []
#B = []
c11_err = []
c12_err = []
c44_err = []

a_vals = []
for model_name in model_names:
    C, C_err, a = get_elastic_constants(model_name, cheap_model_type)
    #print(f'{model_name}: {C[0,0]:.2f} {C[0,1]:.2f} {C[3,3]:.2f}')
    #now write C and C_err to files in output_path
    np.savetxt(f'{output_path}/{model_name}_C.txt', C)
    np.savetxt(f'{output_path}/{model_name}_C_err.txt', C_err)

    #add c11, c12 and c44 to lists
    c11.append(C[0,0])
    c12.append(C[0,1])
    c44.append(C[3,3])
    #B.append((C[0,0]+2*C[0,1])/3)

    c11_err.append(C_err[0,0])
    c12_err.append(C_err[0,1])
    c44_err.append(C_err[3,3])

    a_vals.append(a)

x = np.arange(len(model_names))  # the label locations
width = 0.2 # the width of the bars

bar_positions = [x - 1.5*width, x-0.5*width, x+0.5*width,x + 1.5*width]

fig, ax = plt.subplots()

#plot bar chart with errors
rects1 = ax.bar(bar_positions[0], c11, width, label='C11', yerr=c11_err)
rects2 = ax.bar(bar_positions[1], c12, width, label='C12', yerr=c12_err)
rects3 = ax.bar(bar_positions[2], c44, width, label='C44', yerr=c44_err)
#rects4 = ax.bar(bar_positions[3], B, width, label='B')


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Elastic constants (GPa)')
#ax.set_title(f'{investigation_name} elastic constants')
labels = ['Unconstrained fit', 'Constrained fit']
ax.set_xticks(x-width/2)
ax.set_xticklabels(labels)


#draw horizontal ribbons all the way across marking out the reference values with error
ax.axhline(y=C_ref[0,0], color='tab:blue', linestyle='--', label='C11 ref')
ax.axhline(y=C_ref[0,1], color='tab:orange', linestyle='--', label='C12 ref')
ax.axhline(y=C_ref[3,3], color='tab:green', linestyle='--', label='C44 ref')

#do a fill between the error bars
ax.fill_between(x-width, C_ref[0,0]-C_err_ref[0,0], C_ref[0,0]+C_err_ref[0,0], color='tab:blue', alpha=0.3)
ax.fill_between(x, C_ref[0,1]-C_err_ref[0,1], C_ref[0,1]+C_err_ref[0,1], color='tab:orange', alpha=0.3)
ax.fill_between(x+width, C_ref[3,3]-C_err_ref[3,3], C_ref[3,3]+C_err_ref[3,3], color='tab:green', alpha=0.3)
#rotate x label
#plt.xticks(rotation=45)
plt.title('Silicon (ACE to ACE)')
legend = ax.legend(loc='upper right', bbox_to_anchor=(1, 0.95))
fig.tight_layout()
#put legend in box
plt.savefig(f'{output_path}/elastic_constants.png')

#plot lattice parameter bar chart for reference and each model

fig, ax = plt.subplots()

bar_positions = x
width = 0.35
rects1 = ax.bar(bar_positions, a_vals, width, label='a')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel(r'Lattice parameter ($\AA$)')
#ax.set_title(f'{investigation_name} lattice parameter')
ax.set_xticks(x)
ax.set_xticklabels(labels)

#draw horizontal ribbon all the way across marking out the reference value
ax.axhline(y=a_ref, color='tab:blue', linestyle='--', label='a ref')

# plt.xticks(rotation=45)
ax.legend()
fig.tight_layout()
plt.savefig(f'{output_path}/lattice_parameter.png')

# Write the elastic constants and lattice parameter (and errors) to a file
with open(f'{output_path}/elastic_constants_and_lattice_parameters.txt', 'w') as f:
    f.write(f'Reference elastic constants: {C_ref[0,0]:.2f} {C_ref[0,1]:.2f} {C_ref[3,3]:.2f}\n')
    f.write(f'Reference elastic constants errors: {C_err_ref[0,0]:.2f} {C_err_ref[0,1]:.2f} {C_err_ref[3,3]:.2f}\n')
    f.write(f'Reference lattice parameter: {a_ref:.5f}\n\n')
    
    for i, model_name in enumerate(model_names):
        f.write(f'{model_name} elastic constants: {c11[i]:.2f} {c12[i]:.2f} {c44[i]:.2f}\n')
        f.write(f'{model_name} elastic constants errors: {c11_err[i]:.2f} {c12_err[i]:.2f} {c44_err[i]:.2f}\n')
        f.write(f'{model_name} lattice parameter: {a_vals[i]:.5f}\n\n')


