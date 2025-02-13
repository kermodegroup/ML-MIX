import argparse
import sys
sys.path.append("./potentials/")
sys.path.append("./shared")
# Initial list of all potentials
from load_potentials import load
from write_results import write_results
from evaluate_potentials import run_test

import os
from datetime import datetime

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if size == 1:
    test_type = "Serial"
else:
    test_type = "Parallel"


def test_force_zero(verbose=False,crash_on_fail=True,direct_run=False):
    data_path = "./test_force_zero_data"
    if rank == 0:
        if not os.path.exists(data_path):
            os.makedirs(data_path)
    property_dict_1 = load("./potentials/BASE/BASE.json")
    directory = f"./potentials/LJ"
    json_files = [f for f in os.listdir(directory) if f.endswith(".json")]
    assert len(json_files) == 1, "Only one json file should be in the BASE directory"
    json_file = json_files[0]
    pot_name = json_file[0].split(".")[0]
    property_dict_2 = load(f'{directory}/{json_file}')
    if rank == 0:
        print(f"Running {property_dict_2['ID']}.....", end="", flush=True)
    if verbose:
        err_code = run_test(property_dict_1, 
                            property_dict_2, 
                            verbose=True,
                            zero=True, 
                            comm=comm, 
                            rank=rank, 
                            path=data_path, 
                            crash_on_fail=crash_on_fail)
    else:
        with open(os.devnull, 'w') as devnull:
            sys.stdout = devnull
            err_code = run_test(property_dict_1, 
                                property_dict_2, 
                                zero=True, 
                                comm=comm, 
                                rank=rank, 
                                path=data_path,
                                crash_on_fail=crash_on_fail)
            sys.stdout = sys.__stdout__

    if rank == 0:
        print(err_code)
    
    if direct_run:
        return err_code
    

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

# given the elements, set up the simulation
if __name__ == "__main__":
    args = parse_arguments()
    err_code = test_force_zero(verbose=args.verbose, crash_on_fail=args.crash_on_fail, direct_run=True)
    results = {}
    results[f"Force zeroing"] = {"result":err_code, "datetime":datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
    out_json = f"../feature_test_results_{test_type}.json"
    top_readme = "../README.md"
    insert_key = "feature error table"
    possible_jsons = ["../feature_test_results_Serial.json", "../feature_test_results_Parallel.json"]
    if rank == 0:
        write_results(results, out_json, top_readme, insert_key, possible_jsons)

