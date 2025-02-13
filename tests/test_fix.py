import sys
import ase
from lammps import lammps
sys.path.append("./shared")
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

    return parser.parse_args()


def test_build_regions(verbose=False,crash_on_fail=True,direct_run=False):
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
    i2_potential_python, d2_eval_python = build_regions_python(struct, r_core, r_blend, r_buff, comm=comm, rank=rank, path=data_path)

    if verbose:
        lmps = lammps(comm=comm)
    else:
        lmps = lammps(comm=comm, cmdargs=["-log", "none", "-screen", "none"])
    set_up_lammps(lmps, struct, mass_cmds, multi_potential, rank=rank, path=data_path)
    i2_potential_lammps, d2_eval_lammps = build_regions_lammps(lmps, struct, r_core, r_blend, r_buff, comm=comm, rank=rank, path=data_path)

    i2_potential_combined = i2_potential_python^i2_potential_lammps

    d2_eval_diff = d2_eval_python - d2_eval_lammps
    d2_eval_diff = np.abs(d2_eval_diff)
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


    

if __name__ == "__main__":
    args = parse_arguments()
    err_code = test_build_regions(verbose=args.verbose, crash_on_fail=args.crash_on_fail, direct_run=True)
    results = {"Region building": {"result":err_code, "datetime":datetime.now().strftime("%Y-%m-%d %H:%M:%S")}}
    out_json = f"../feature_test_results_{test_type}.json"
    top_readme = "../README.md"
    insert_key = "feature error table"
    possible_jsons = ["../feature_test_results_Serial.json", "../feature_test_results_Parallel.json"]
    if rank == 0:
        write_results(results, out_json, top_readme, insert_key, possible_jsons)