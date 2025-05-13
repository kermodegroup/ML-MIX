import sys
import ase
from lammps import lammps
sys.path.append("./shared")
sys.path.append("./potentials/")
from load_potentials import load
# Initial list of all potentials
from build_struct import build_struct
from build_regions import build_regions_python, build_regions_lammps
from set_up_lammps import set_up_lammps
from write_results import write_results
import numpy as np
import json
import os
from datetime import datetime
import re
import argparse

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if size == 1:
    test_type = "Serial"
else:
    test_type = "Parallel"



def parse_arguments():
    parser = argparse.ArgumentParser(description="Select or exclude potentials for a simulation.")
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Direct output to stdout"
    )
    parser.add_argument(
        "--crash_on_fail",
        action="store_true",
        help="Crash on the first failure"
    )
    parser.add_argument(
        "--pick_seed_with",
        type=str,
        default="group",
        help="Method to pick seed atoms"
    )
    parser.add_argument(
        "--blend_type",
        type=str,
        default="linear",
        help="Type of blending to use"
    )

    return parser.parse_args()


def region_test(verbose=False,
                crash_on_fail=True,
                direct_run=False,
                pick_seed_with='group',
                nsteps=1,
                fix_nevery=10,
                blend_type='linear'):
    data_path = "./test_fix_data"
    if rank == 0:
        if not os.path.exists(data_path):
            os.makedirs(data_path)

    elements = ['H']
    mass_cmds = ['mass 1 1.0']
    multi_potential = True
    r_core = 5.0
    r_blend = 3.0
    r_buff = 3.0
    struct = build_struct(elements, path=data_path)
    if rank == 0:
        print("running region build test...", end="")
    
    if pick_seed_with == 'fix_and_init_group':
        if nsteps > fix_nevery:
            python_pick_seed_with = 'fix'
        else:
            python_pick_seed_with = 'group'
    elif pick_seed_with == 'fix':
        if nsteps > fix_nevery:
            python_pick_seed_with = 'fix'
        else:
            python_pick_seed_with = 'all_expensive'
    elif pick_seed_with == 'group':
        python_pick_seed_with = 'group'
    
    i2_potential_python, d2_eval_python = \
        build_regions_python(struct, 
                             r_core, 
                             r_blend, 
                             r_buff, 
                             comm=comm, 
                             rank=rank, 
                             path=data_path,
                             pick_seed_with=python_pick_seed_with,
                             blend_type=blend_type,)

    if verbose:
        lmps = lammps(comm=comm)
    else:
        lmps = lammps(comm=comm, cmdargs=["-log", "none", "-screen", "none"])
    set_up_lammps(lmps, struct, mass_cmds, multi_potential, rank=rank, path=data_path)
    pot = load("./potentials/BASE/BASE.json")
    cmds_1 = pot["cmds"]
    lmps.command(f'pair_style {cmds_1["style_name"]} {cmds_1["style_params"]}')
    lmps.command(f'pair_coeff {cmds_1["coeff_types"]} {cmds_1["coeff_params"]}')
    i2_potential_lammps, d2_eval_lammps = \
        build_regions_lammps(lmps, 
                             struct, 
                             r_core, 
                             r_blend, 
                             r_buff, 
                             comm=comm, 
                             rank=rank, 
                             path=data_path,
                             pick_seed_with=pick_seed_with,
                             nsteps=nsteps,
                             fix_nevery=fix_nevery,
                             blend_type=blend_type,)


    d2_eval_diff = d2_eval_python - d2_eval_lammps
    d2_eval_diff = np.abs(d2_eval_diff)
    struct.arrays['i2_potential_python'] = i2_potential_python
    struct.arrays['i2_potential_lammps'] = i2_potential_lammps
    struct.arrays['d2_eval_python'] = d2_eval_python
    struct.arrays['d2_eval_lammps'] = d2_eval_lammps
    struct.arrays['d2_eval_diff'] = d2_eval_diff
    if rank == 0:
        ase.io.write(f"{data_path}/debug.xyz", struct, format='extxyz', parallel=False)
    
    try:
        assert np.all(i2_potential_python == i2_potential_lammps)
        assert np.allclose(d2_eval_python, d2_eval_lammps, atol=1e-6)
        err_code = "✅"
    except AssertionError:
        err_code = "❌"
        if crash_on_fail:
            raise AssertionError("Test failed!")

    if rank == 0:
        print(err_code)
    lmps.close()
    if direct_run:
        return err_code
    

def test_build_with_group():
    region_test(pick_seed_with='group')
    region_test(pick_seed_with='group',blend_type='cubic')

def test_build_with_fix():
    region_test(pick_seed_with='fix',nsteps=1,fix_nevery=10)
    region_test(pick_seed_with='fix',nsteps=5,fix_nevery=10)
    region_test(pick_seed_with='fix',nsteps=20,fix_nevery=10)
    region_test(pick_seed_with='fix',nsteps=55,fix_nevery=10)
    region_test(pick_seed_with='fix',nsteps=2,fix_nevery=1,blend_type='cubic')


def test_build_with_init_group():
    region_test(pick_seed_with='fix_and_init_group',nsteps=1,fix_nevery=10)
    region_test(pick_seed_with='fix_and_init_group',nsteps=5,fix_nevery=10)
    region_test(pick_seed_with='fix_and_init_group',nsteps=20,fix_nevery=10)
    region_test(pick_seed_with='fix_and_init_group',nsteps=55,fix_nevery=10)
    region_test(pick_seed_with='fix_and_init_group',nsteps=2,fix_nevery=1,blend_type='cubic')
    region_test(pick_seed_with='fix_and_init_group',nsteps=55,fix_nevery=10,blend_type='cubic')



    

if __name__ == "__main__":
    args = parse_arguments()
    if args.pick_seed_with not in ['group', 'fix', 'fix_and_init_group']:
        raise ValueError("pick_seed_with must be 'group', 'fix' or 'fix_and_init_group'")
    err_code = region_test(verbose=args.verbose, 
                           crash_on_fail=args.crash_on_fail,
                           direct_run=True, 
                           pick_seed_with=args.pick_seed_with, 
                           nsteps=1,
                           fix_nevery=10,
                           blend_type=args.blend_type)
    results = {"Region building": {"result":err_code, "datetime":datetime.now().strftime("%Y-%m-%d %H:%M:%S")}}
    out_json = f"../feature_test_results_{test_type}.json"
    top_readme = "../README.md"
    insert_key = "feature error table"
    possible_jsons = ["../feature_test_results_Serial.json", "../feature_test_results_Parallel.json"]
    if rank == 0:
        write_results(results, out_json, top_readme, insert_key, possible_jsons)