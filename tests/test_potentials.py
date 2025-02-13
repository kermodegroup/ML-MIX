import argparse
import sys
from lammps import lammps
sys.path.append("./potentials/")
sys.path.append("./shared")
# Initial list of all potentials
from load_potentials import load
from write_results import write_results
from evaluate_potentials import run_test

import os
import numpy as np
from datetime import datetime

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if size == 1:
    test_type = "Serial"
else:
    test_type = "Parallel"

#get names of all potentials from folder names in potential folder
# as long as they contain a json file
all_potentials = [f for f in os.listdir("./potentials/") if os.path.isdir(f"./potentials/{f}")]
all_potentials = [potential for potential in all_potentials if len([f for f in os.listdir(f"./potentials/{potential}") if f.endswith(".json")]) > 0]


def parse_arguments():
    parser = argparse.ArgumentParser(description="Select or exclude potentials for a simulation.")

    # Add mutually exclusive arguments for including or excluding potentials
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--potentials",
        "-p",
        nargs="+", 
        help="List of potentials to include (e.g., --potentials LJ ACE)"
    )
    group.add_argument(
        "--exclude_pot", 
        nargs="+", 
        help="List of potentials to exclude (e.g., --exclude_pot UF3 SW)"
    )
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

def process_potentials(args):
    if args.potentials:
        # Use specified potentials
        selected_potentials = args.potentials
    elif args.exclude_pot:
        # Exclude specified potentials
        selected_potentials = [pot for pot in all_potentials if pot not in args.exclude_pot]
    else:
        selected_potentials = all_potentials

    # delete BASE
    if "BASE" in selected_potentials:
        selected_potentials.remove("BASE")

    return selected_potentials


def pot_test(selected_potentials, verbose=False, crash_on_fail=True, direct_run=False):
    data_path = "./test_potentials_data"
    if rank == 0:
        print("Selected potentials:", selected_potentials)
        if not os.path.exists(data_path):
            os.makedirs(data_path)
    property_dict_1 = load("./potentials/BASE/BASE.json")
    results = {}
    #for each potential in selected_potentials, load json files in the potential folder
    for potential in selected_potentials:
        directory = f"./potentials/{potential}"
        json_files = [f for f in os.listdir(directory) if f.endswith(".json")]
        for json_file in json_files:
            pot_name = json_file[0].split(".")[0]
            property_dict_2 = load(f'{directory}/{json_file}')
            #print running property_dict_2["ID"].....
            if rank == 0:
                print(f"Running {property_dict_2['ID']}.....", end="", flush=True)
            if verbose:
                err_code = run_test(property_dict_1,
                                    property_dict_2, 
                                    verbose=True, 
                                    comm=comm, 
                                    rank=rank, 
                                    path=data_path,
                                    crash_on_fail=crash_on_fail)
            else:
                with open(os.devnull, 'w') as devnull:
                    sys.stdout = devnull
                    err_code = run_test(property_dict_1, 
                                        property_dict_2, 
                                        comm=comm, 
                                        rank=rank, 
                                        path=data_path,
                                        crash_on_fail=crash_on_fail)
                    sys.stdout = sys.__stdout__

            if rank == 0:
                print(err_code)
            
            results[f"{property_dict_2['ID']}"] = {"result":err_code, "datetime":datetime.now().strftime("%Y-%m-%d %H:%M:%S")}
        
    if direct_run:
        return results


# Pytest specific tests for each potential in turn
def test_LJ():
    selected_potentials = ["LJ"]
    pot_test(selected_potentials)

def test_ACE():
    selected_potentials = ["ACE"]
    pot_test(selected_potentials)

def test_UF3():
    selected_potentials = ["UF3"]
    pot_test(selected_potentials)

def test_table():
    selected_potentials = ["TABLE"]
    pot_test(selected_potentials)

def test_EAM():
    selected_potentials = ["EAM"]
    pot_test(selected_potentials)


# given the elements, set up the simulation
if __name__ == "__main__":
    args = parse_arguments()
    selected_potentials = process_potentials(args)
    results = pot_test(selected_potentials, verbose=args.verbose, crash_on_fail=args.crash_on_fail, direct_run=True, data_path=data_path)
    out_json = f"test_results_{test_type}.json"
    top_readme = "../README.md"
    insert_key = "error table"
    possible_jsons = ["test_results_Serial.json", "test_results_Parallel.json"]
    if rank == 0:
        write_results(results, out_json, top_readme, insert_key, possible_jsons)


