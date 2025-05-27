import numpy as np
import ase
from lammps import lammps
import sys
sys.path.append("../shared")
from build_struct import build_struct
from build_regions import build_regions_lammps
from set_up_lammps import set_up_lammps

# class LammpsMock:
#     def __init__(self, input_file=f'catch_commands.in'):
#         self.input_file = input_file
#         with open(self.input_file, 'w') as f:
#             f.write("# LAMMPS input file\n")

#     def command(self, cmd):
#         with open(self.input_file, 'a') as f:
#             f.write(cmd + '\n')

#     def commands_list(self, cmds):
#         for cmd in cmds:
#             self.command(cmd)


def run_test(property_dict_1, property_dict_2, verbose=False, zero=False, comm=None, rank=None, path='./', crash_on_fail=False, kokkos=False):
    # generate structure
    if kokkos:
        dump_string = 'd_potential_1 d_potential_2 d_eval_1 d_eval_2'
    else:
        dump_string = 'i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]'

    elements_1 = property_dict_1["elements"]
    elements_2 = property_dict_2["elements"]
    cmds_1 = property_dict_1["cmds"]
    cmds_2 = property_dict_2["cmds"]
    ID_1 = property_dict_1["ID"]
    ID_2 = property_dict_2["ID"]
    cutoff_1 = float(property_dict_1["cutoff"])
    cutoff_2 = float(property_dict_2["cutoff"])

    if ("X" in elements_1) or ("X" in elements_2):
        if ("X" in elements_1) and ("X" in elements_2):
            elements = ["H"]
            a = 3.0
            structure = "fcc"
        elif "X" in elements_1:
            elements = elements_2
            a = property_dict_2["a"]
            structure = property_dict_2["structure"]
        else:
            elements = elements_1
            a = property_dict_1["a"]
            structure = property_dict_1["structure"]
    else:
        #elements 1 and elements 2 must contain the same elements
        for element in elements_1:
            if element not in elements_2:
                raise ValueError("Elements in potentials must match.")
        elements = elements_1
        a = property_dict_1["a"]
        structure = property_dict_1["structure"]
    if structure == "X":
        structure = "fcc"
    if a == "X":
        a = 3.0
    
    #generate mass array for LAMMPS from elements
    mass_cmds = []
    for i, element in enumerate(elements):
        atomic_num = ase.data.atomic_numbers[element]
        atomic_mass = ase.data.atomic_masses[atomic_num]
        mass_cmds.append(f"mass {i+1} {atomic_mass}")
    struct = build_struct(elements, rank=rank, path=path, structure=structure, a=a)
    if kokkos:
        args = ["-k", "on", "g","1",
                "-sf", "kk",
                "-pk", "kokkos", "newton", "on", "neigh", "half",
]
    else:
        args = []

    if not verbose:
        args = args + ["-log", "none", "-screen", "none"]

    lmps = lammps(cmdargs=args, comm=comm)

    set_up_lammps(lmps, struct, mass_cmds,path=path, kokkos=kokkos)
    
    r_core = 3.0
    r_blend = 3.0
    r_buff = np.max([cutoff_1, cutoff_2])
    i2_potential, d2_eval = build_regions_lammps(lmps, struct, r_core, r_blend, r_buff, path=path, kokkos=kokkos)
    #
    # Unfixing here is necessary as it turns out that if you have just fix mlml and a pair_style that
    # requests a half full neighbour list, the pair_style tries to get it from the neighbourlist
    # used in fix mlml, and for some reason it doesn't work (all forces are 0). This isn't a problem
    # if you have another pair_style defined that uses the full neighbour list. But it means that right now,
    # the LJ/LJ tests fail. This is a workaround so LJ can be used as one of the pair_styles in the
    # other tests. For now, just removing LJ test until this is fixed.
    # 
    lmps.command('unfix mlml_fix')
    try:
        #get forces for each potential
        # lmps.command('compute myforce all reduce sum fx fy fz')
        # lmps.command('variable fmag equal sqrt(c_myforce[1]^2+c_myforce[2]^2+c_myforce[3]^2)')
        # lmps.command('thermo_style custom step pe ke v_fmag')
        lmps.command(f'pair_style {cmds_1["style_name"]} {cmds_1["style_params"]}')
        lmps.command(f'pair_coeff {cmds_1["coeff_types"]} {cmds_1["coeff_params"]}')
        lmps.command('run 0')
        lmps.command(f'write_dump all custom {path}/just_1.lammpstrj id type xs ys zs fx fy fz modify format float %20.15g')

        lmps.command(f'pair_style {cmds_2["style_name"]} {cmds_2["style_params"]}')
        lmps.command(f'pair_coeff {cmds_2["coeff_types"]} {cmds_2["coeff_params"]}')
        lmps.command('run 0')
        lmps.command(f'write_dump all custom {path}/just_2.lammpstrj id type xs ys zs fx fy fz modify format float %20.15g')
        
        lmps.command(f'fix mlml_fix all mlml 1 {r_core} {r_buff} {r_blend} group seed_atoms')
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

        lmps.command('run 0')
        lmps.command(f'write_dump all custom {path}/both.lammpstrj id type xs ys zs fx fy fz {dump_string} modify format float %20.15g')


    except:
        if crash_on_fail:
            raise AssertionError("Test crashed!")
        else:
            return "ERR"
    
    lmps.close()

    if rank == 0:
        just_1 = ase.io.read(f'{path}/just_1.lammpstrj', parallel=False)
        just_2 = ase.io.read(f'{path}/just_2.lammpstrj', parallel=False)
        both = ase.io.read(f'{path}/both.lammpstrj', parallel=False)
        
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
        print(np.max(np.abs(diff)))
        both.arrays['fb'] = forces_both
        both.arrays['fmb'] = manual_forces_both
        both.arrays['diff'] = diff
        ase.io.write(f'{path}/diff.xyz', both, format='extxyz', parallel=False)
        try:
            assert np.allclose(forces_both, manual_forces_both, atol=1e-5, rtol=0)
        except AssertionError:
            print(f"Forces for {ID_1} do not match.")
            if crash_on_fail:
                raise AssertionError("Test failed!")
            return "❌"
        return "✅"
    
    # if not rank 0
    return None

