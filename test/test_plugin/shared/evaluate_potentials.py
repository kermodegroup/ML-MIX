import numpy as np
import ase
from lammps import lammps
import sys
sys.path.append("../shared")
from build_struct import build_struct
from build_regions import build_regions_lammps
from set_up_lammps import set_up_lammps

def run_test(property_dict_1, property_dict_2, verbose=False, zero=False, comm=None, rank=None):
    # generate structure

    elements_1 = property_dict_1["elements"]
    elements_2 = property_dict_2["elements"]
    cmds_1 = property_dict_1["cmds"]
    cmds_2 = property_dict_2["cmds"]
    ID_1 = property_dict_1["ID"]
    ID_2 = property_dict_2["ID"]

    if ("X" in elements_1) or ("X" in elements_2):
        if ("X" in elements_1) and ("X" in elements_2):
            elements = ["H"]
        elif "X" in elements_1:
            elements = elements_2
        else:
            elements = elements_1
    else:
        #elements 1 and elements 2 must contain the same elements
        for element in elements_1:
            if element not in elements_2:
                raise ValueError("Elements in potentials must match.")
        elements = elements_1
    
    #generate mass array for LAMMPS from elements
    mass_cmds = []
    for i, element in enumerate(elements):
        atomic_num = ase.data.atomic_numbers[element]
        atomic_mass = ase.data.atomic_masses[atomic_num]
        mass_cmds.append(f"mass {i+1} {atomic_mass}")
    struct = build_struct(elements, rank=rank)
    if not verbose:
        lmps = lammps(cmdargs=["-log", "none", "-screen", "none"], comm=comm)
    else:
        lmps = lammps(comm=comm)
    set_up_lammps(lmps, struct, mass_cmds)
    
    r_core = 3.0
    r_blend = 3.0
    r_buff = 5.0
    i2_potential, d2_eval = build_regions_lammps(lmps, struct, r_core, r_blend, r_buff)
    try:
        #get forces for each potential
        lmps.command(f'pair_style {cmds_1["style_name"]} {cmds_1["style_params"]}')
        lmps.command(f'pair_coeff {cmds_1["coeff_types"]} {cmds_1["coeff_params"]}')
        lmps.command('dump myDump all custom 1 just_1.lammpstrj id type xs ys zs fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]')
        lmps.command('run 0')
        lmps.command(f'undump myDump')

        lmps.command(f'pair_style {cmds_2["style_name"]} {cmds_2["style_params"]}')
        lmps.command(f'pair_coeff {cmds_2["coeff_types"]} {cmds_2["coeff_params"]}')
        lmps.command('dump myDump all custom 1 just_2.lammpstrj id type xs ys zs fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]')
        lmps.command('run 0')
        lmps.command(f'undump myDump')

        if zero:
            lmps.command(f'pair_style hybrid/overlay/mlml zero yes {cmds_1["style_name"]} {cmds_1["style_params"]} {cmds_2["style_name"]} {cmds_2["style_params"]}')
        else:
            lmps.command(f'pair_style hybrid/overlay/mlml {cmds_1["style_name"]} {cmds_1["style_params"]} {cmds_2["style_name"]} {cmds_2["style_params"]}')
            
        if ID_1 == ID_2:
            lmps.command(f'pair_coeff {cmds_1["coeff_types"]} {cmds_1["style_name"]} 1 1 {cmds_1["coeff_params"]}')
            lmps.command(f'pair_coeff {cmds_2["coeff_types"]} {cmds_2["style_name"]} 2 2 {cmds_2["coeff_params"]}')
        else:
            lmps.command(f'pair_coeff {cmds_1["coeff_types"]} {cmds_1["style_name"]} 1 {cmds_1["coeff_params"]}')
            lmps.command(f'pair_coeff {cmds_2["coeff_types"]} {cmds_2["style_name"]} 2 {cmds_2["coeff_params"]}')

        lmps.command('dump myDump all custom 1 both.lammpstrj id type xs ys zs fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]')
        lmps.command('run 0')
        lmps.command(f'undump myDump')

    except:
        return "ERR"
    
    lmps.close()

    if rank == 0:
        just_1 = ase.io.read('just_1.lammpstrj', parallel=False)
        just_2 = ase.io.read('just_2.lammpstrj', parallel=False)
        both = ase.io.read('both.lammpstrj', parallel=False)
        
        #compare forces
        forces_1 = just_1.get_forces()
        forces_2 = just_2.get_forces()
        forces_both = both.get_forces()

        #check forces_**d2_eval[:,*] == forces_both[:,*]
        manual_forces_both = forces_1*d2_eval[:,0].reshape(-1,1) \
                            + forces_2*d2_eval[:,1].reshape(-1,1)
        
        if zero:
            #get the vector sum of manual_forces_both
            total_force = np.sum(manual_forces_both, axis=0)
            #subtract this from forces_both
            manual_forces_both -= total_force/np.shape(manual_forces_both)[0]
            
        diff = manual_forces_both - forces_both
        print(np.max(diff))
        both.arrays['fb'] = forces_both
        both.arrays['fmb'] = manual_forces_both
        both.arrays['diff'] = diff
        ase.io.write('diff.xyz', both, format='extxyz', parallel=False)
        try:
            assert np.allclose(forces_both, manual_forces_both, atol=1e-6)
        except AssertionError:
            print(f"Forces for {ID_2} do not match.")
            return "❌"
        return "✅"
    
    # if not rank 0
    return None

