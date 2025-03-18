import argparse
import sys
from lammps import lammps
sys.path.append("./potentials/")
sys.path.append("./shared")
# Initial list of all potentials
from load_potentials import load
from write_results import write_results
from evaluate_potentials import run_test

from set_up_lammps import set_up_lammps
from build_regions import build_regions_lammps

import os
import numpy as np
from datetime import datetime

from mpi4py import MPI
import ase.io


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if size == 1:
    test_type = "Serial"
else:
    test_type = "Parallel"

if rank == 0:
    # make on-fly-data directory if it does not exist
    if not os.path.exists("./on-fly-data"):
        os.makedirs("./on-fly-data")
path = "./on-fly-data"

lmps = lammps()

def single_eval(struct, mass_cmds, pair_coeff_cmds, dump_name):
    set_up_lammps(lmps, struct, mass_cmds, multi_potential=False, rank=0, path='./')
    lmps.command(f'pair_style hybrid/overlay pace table spline 5000')
    lmps.commands_list(pair_coeff_cmds)
    lmps.command('run 1')
    lmps.command(f'write_dump all custom {path}/{dump_name}.lammpstrj id type xs ys zs fx fy fz modify format float %20.15g')

def run_and_fit(struct, mass_cmds, config_path, pair_coeff_cmds):
    #get id of atom closest to center
    set_up_lammps(lmps, struct, mass_cmds, multi_potential=True, rank=0, path='./')
    i2_potential, d2_eval = build_regions_lammps(lmps, struct, 6.0, 4.0, 5.51, pick_seed_with='group')
    lmps.command(f'pair_style hybrid/overlay/mlml on-fly 5 ./../ {config_path} pace pace table spline 5000 table spline 5000')
    lmps.commands_list(pair_coeff_cmds)
    lmps.command('run 1')
    lmps.command(f'write_dump all custom {path}/pre_fit.lammpstrj id type xs ys zs fx fy fz d2_eval[1] d2_eval[2] i_potential[1] i_potential[2] modify format float %20.15g')
    lmps.command('run 10')
    lmps.command(f'write_dump all custom {path}/post_fit.lammpstrj id type xs ys zs fx fy fz d2_eval[1] d2_eval[2] i_potential[1] i_potential[2] modify format float %20.15g')

    return i2_potential, d2_eval


def test_on_fly_fit():
    struct = ase.io.read("./on-fly/input_structs/fe_dumbell_16.xyz")
    struct.rattle(0.01)
    mass_cmds = ['mass 1 55.845']
    config_path = "./on-fly/config_files/Fe_otf_config_1.yaml"
    pair_coeff_cmds = [
        'pair_coeff * * pace 1 1 on-fly/expensive_pots/jace_espresso.yace Fe',
        'pair_coeff 1 1 table 1 1 on-fly/expensive_pots/jace_espresso_pairpot.table Fe_Fe',
        'pair_coeff * * pace 2 2 otf_pot.yace Fe',
        'pair_coeff 1 1 table 2 2 otf_pot_pairpot.table Fe_Fe'
    ]
    i2_potential, d2_eval = run_and_fit(struct, mass_cmds, config_path, pair_coeff_cmds)

    pair_coeff_cmds = [
        'pair_coeff * * pace on-fly/expensive_pots/jace_espresso.yace Fe',
        'pair_coeff 1 1 table on-fly/expensive_pots/jace_espresso_pairpot.table Fe_Fe',
    ]
    single_eval(struct, mass_cmds, pair_coeff_cmds, "just_1")

    pair_coeff_cmds = [
        'pair_coeff * * pace otf_pot.yace Fe',
        'pair_coeff 1 1 table otf_pot_pairpot.table Fe_Fe'
    ]
    single_eval(struct, mass_cmds, pair_coeff_cmds, "just_2")

    just_1 = ase.io.read(f'{path}/just_1.lammpstrj', parallel=False)
    just_2 = ase.io.read(f'{path}/just_2.lammpstrj', parallel=False)
    pre_fit = ase.io.read(f'{path}/pre_fit.lammpstrj', parallel=False)
    post_fit = ase.io.read(f'{path}/post_fit.lammpstrj', parallel=False)

    #compare forces
    forces_1 = just_1.get_forces()
    forces_2 = just_2.get_forces()
    forces_pre_fit = pre_fit.get_forces()
    forces_post_fit = post_fit.get_forces()
    diff = forces_pre_fit - forces_1
    diff_over_fit = forces_post_fit - forces_pre_fit

    assert np.allclose(forces_pre_fit, forces_1, atol=1e-5, rtol=0)
    manual_forces_post = forces_1*d2_eval[:,0].reshape(-1,1) \
                        + forces_2*d2_eval[:,1].reshape(-1,1)

    struct.arrays['diff'] = manual_forces_post - forces_post_fit
    ase.io.write(f'{path}/diff.xyz', struct, format='extxyz', parallel=False)

    assert np.allclose(forces_post_fit, manual_forces_post, atol=1e-5, rtol=0)


# test_on_fly_fit()
# if rank == 0:
#     print("✅")
# lmps.close()