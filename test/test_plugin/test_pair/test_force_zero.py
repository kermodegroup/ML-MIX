import argparse
import sys
sys.path.append("./potentials/")
sys.path.append("../shared")
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


def parse_arguments():
    parser = argparse.ArgumentParser(description="Select or exclude potentials for a simulation.")

    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Direct output to stdout"
    )

    return parser.parse_args()

# given the elements, set up the simulation
if __name__ == "__main__":
    args = parse_arguments()
    results = {}
    property_dict_1 = load("./potentials/BASE/BASE.json")
    directory = f"./potentials/LJ"
    json_files = [f for f in os.listdir(directory) if f.endswith(".json")]
    for json_file in json_files:
        pot_name = json_file[0].split(".")[0]
        property_dict_2 = load(f'{directory}/{json_file}')
        if rank == 0:
            print(f"Running {property_dict_2['ID']}.....", end="", flush=True)
        if args.verbose:
            err_code = run_test(property_dict_1, property_dict_2, verbose=True, zero=True, comm=comm, rank=rank)
        else:
            with open(os.devnull, 'w') as devnull:
                sys.stdout = devnull
                err_code = run_test(property_dict_1, property_dict_2, zero=True, comm=comm, rank=rank)
                sys.stdout = sys.__stdout__

        if rank == 0:
            print(err_code)
        results[f"Force zeroing"] = {"result":err_code, "datetime":datetime.now().strftime("%Y-%m-%d %H:%M:%S")}


    out_json = f"../feature_test_results_{test_type}.json"
    top_readme = "../../../README.md"
    insert_key = "feature error table"
    possible_jsons = ["../feature_test_results_Serial.json", "../feature_test_results_Parallel.json"]
    if rank == 0:
        write_results(results, out_json, top_readme, insert_key, possible_jsons)

