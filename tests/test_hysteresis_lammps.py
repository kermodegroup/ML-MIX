
import sys
sys.path.append("./shared")
from set_up_lammps import set_up_lammps
from build_regions import build_regions_lammps, get_seed_atoms, build_regions_python
from lammps import lammps
import ase.io
from mpi4py import MPI
import numpy as np
sys.path.append("../utils")
from utils import read_lammps_dump
import os

def predict_d2_eval(d2_eval_prev, d2_eval_target, nevery, dt, tau_in, tau_out):
    change = ((d2_eval_target - d2_eval_prev))

    in_multiplier = dt*nevery/tau_in
    out_multiplier = dt*nevery/tau_out

    if in_multiplier > 1.0:
        in_multiplier = 1.0
    if out_multiplier > 1.0:
        out_multiplier = 1.0
    
    change[change>0] = change[change>0] * in_multiplier
    change[change<0] = change[change<0] * out_multiplier

    new_d2_eval = change + d2_eval_prev

    change_less_0_mask = change<0
    change_greater_0_mask = change>0
    new_d2_eval_less_thresh_mask = new_d2_eval<0.01
    new_d2_eval_greater_thresh_mask = new_d2_eval>0.99
    
    new_d2_eval[change_less_0_mask & new_d2_eval_less_thresh_mask] = 0.0
    new_d2_eval[change_greater_0_mask & new_d2_eval_greater_thresh_mask] = 1.0
    return new_d2_eval


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
if size == 1:
    test_type = "Serial"
else:
    test_type = "Parallel"


def hysteresis_test(verbose=False,
                    crash_on_fail=True,
                    direct_run=False,
                    pick_seed_with='group',
                    nsteps=1,
                    fix_nevery=10,
                    hysteresis_time_in=0.001,
                    hysteresis_time_out=0.01,
                    nevery=1):
    
    data_path = "./test_hysteresis_data"
    if rank == 0:
        if not os.path.exists(data_path):
            os.makedirs(data_path)

    mass_W = 183.84
    mass_He = 4.002602
    mass_cmds = [f'mass 1 {mass_W}',f'mass 2 {mass_He}']
    r_core = 5.0
    r_blend = 5.0
    r_buff = 5.0
    dt = 0.001

    if rank == 0:
        struct = ase.io.read("W_He_input_16_unrelaxed.xyz", parallel=False)
    else:
        struct = None

    struct = comm.bcast(struct, root=0)

    if verbose:
        lmps = lammps(comm=comm)
    else:
        lmps = lammps(comm=comm, cmdargs=["-log", "none", "-screen", "none"])
    set_up_lammps(lmps, struct, mass_cmds, multi_potential=True, rank=rank, path=data_path)

    lmps.command(f"pair_style zero 5.0")
    lmps.command(f"pair_coeff * *")

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
                                fix_lb=21,
                                hysteresis=True,
                                hysteresis_time_in=hysteresis_time_in,
                                hysteresis_time_out=hysteresis_time_out,
                                nevery=nevery)

    if pick_seed_with == 'group':
        lmps.command('delete_atoms group seed_atoms compress no')
        seed_atoms = get_seed_atoms(struct)
        d2_eval_lammps = np.delete(d2_eval_lammps, seed_atoms, axis=0)

        
    elif pick_seed_with == 'fix':
        lmps.command('group He_atoms type 1')
        lmps.command('delete_atoms group He_atoms compress no')
        d2_eval_lammps = np.delete(d2_eval_lammps, -1, axis=0)

    elif pick_seed_with == 'fix_and_init_group':
        lmps.command('group He_atoms type 1')
        # lmps.command('delete_atoms group seed_atoms compress no')
        lmps.command('delete_atoms group He_atoms compress no')
        seed_atoms = get_seed_atoms(struct)
        # d2_eval_lammps = np.delete(d2_eval_lammps, seed_atoms, axis=0)
        d2_eval_lammps = np.delete(d2_eval_lammps, -1, axis=0)


    lmps.command(f'dump decaydump all custom 1 {data_path}/dump_{pick_seed_with}.lammpstrj id type x y z fx fy fz i2_potential[1] i2_potential[2] d2_eval[1] d2_eval[2]')
    lmps.command('run 20')

    if rank == 0:
        dump = read_lammps_dump(f'{data_path}/dump_{pick_seed_with}.lammpstrj')

        try:
            # first check the initial step matches
            if pick_seed_with == 'group':
                d2_eval_prev = d2_eval_lammps[:,0]
                d2_eval_target = np.zeros_like(d2_eval_prev)
                d2_eval_predicted = predict_d2_eval(d2_eval_prev, d2_eval_target, nevery, dt, hysteresis_time_in, hysteresis_time_out)
                next_d2_eval = dump[0].arrays['d2_eval[1]']
                assert np.allclose(d2_eval_predicted, next_d2_eval), f"Prediction failed at step 0"
            elif pick_seed_with == 'fix':
                assert np.allclose(dump[0].arrays['d2_eval[1]'], np.ones_like(d2_eval_lammps[:,0])), f"Prediction failed at step 0"



            for i in range(len(dump)-1):
                if pick_seed_with != 'group':
                    if i < fix_nevery:
                        continue
                if (nsteps+i+1) % nevery != 0:
                    continue
                image = dump[i]
                d2_eval_prev = image.arrays['d2_eval[1]']
                d2_eval_target = np.zeros_like(d2_eval_prev)
                d2_eval_predicted = predict_d2_eval(d2_eval_prev, d2_eval_target, nevery, dt, hysteresis_time_in, hysteresis_time_out)
                next_d2_eval = dump[i+1].arrays['d2_eval[1]']
                # write all three cols to a file
                assert np.allclose(d2_eval_predicted, next_d2_eval), f"Prediction failed at step {i+1}"
        except AssertionError:
            err_code = "❌"
            if crash_on_fail:
                raise AssertionError("Test failed!")
            

    #finally check to make sure that when atoms are added to group, it flickers in according to hysteresis time in

    if pick_seed_with == "group":
        seed_id = 1000
        struct = dump[-1]
        i2_potential_python, d2_eval_python = build_regions_python(struct, r_core, r_blend, r_buff, pick_seed_with='group', comm=comm, rank=rank, path=data_path, seed_id=seed_id)
        lmps.command(f'group seed_atoms delete')
        lmps.command(f'group seed_atoms id {seed_id+1}')
        lmps.command(f'run 10')
        if rank == 0:
            try:
                dump = read_lammps_dump(f'{data_path}/dump_{pick_seed_with}.lammpstrj')[-10:]
                for i in range(len(dump)-1):
                    # print("curr step", nsteps+20+i, "out of", nsteps+20+len(dump)-1)
                    # print("next step", nsteps+20+i+1)
                    if (nsteps+20+i) % nevery != 0:
                        # print("skipping")
                        continue
                    image = dump[i]
                    d2_eval_prev = image.arrays['d2_eval[1]']
                    d2_eval_predicted = predict_d2_eval(d2_eval_prev, d2_eval_python[:,0], nevery, dt, hysteresis_time_in, hysteresis_time_out)
                    next_d2_eval = dump[i+1].arrays['d2_eval[1]']
                    # #save to numpy txt file
                    # full_arr = np.column_stack((d2_eval_python[:,0], d2_eval_predicted))
                    # diff = np.abs(d2_eval_predicted - next_d2_eval)
                    # sum_diff = np.sum(diff)
                    # print(f"max diff: {np.max(diff)}")
                    # sum_prev = np.sum(np.abs(next_d2_eval-d2_eval_prev))
                    # print(f"sum prev: {sum_prev}")
                    # struct.arrays['d2_eval_predicted'] = d2_eval_predicted
                    # struct.arrays['d2_eval_next'] = next_d2_eval
                    # struct.arrays['d2_eval_target'] = d2_eval_python[:,0]
                    # struct.arrays['d2_eval_diff'] = diff
                    # ase.io.write(f"{data_path}/debug_{i}.xyz", struct, format='extxyz', parallel=False)
                    # print(f"Sum of differences: {sum_diff}")
                    # np.savetxt(f"{data_path}/predicted_{i}.txt", full_arr, header="id d2_eval_python d2_eval_predicted")
                    assert np.allclose(d2_eval_predicted, next_d2_eval, atol=1e-4), f"Prediction failed at step {i+1}"
            except AssertionError:
                err_code = "❌"
                if crash_on_fail:
                    raise AssertionError("Test failed!")
    
    lmps.close()

    if direct_run:
        return err_code


def test_hysteresis_group():
    hysteresis_test(pick_seed_with='group', hysteresis_time_out=0.01, nevery=1)
    hysteresis_test(pick_seed_with='group', hysteresis_time_out=0.002, nevery=1)
    hysteresis_test(pick_seed_with='group', hysteresis_time_out=0.001, nevery=5)
    hysteresis_test(pick_seed_with='group', hysteresis_time_out=0.05, nevery=10)
    hysteresis_test(pick_seed_with='group', hysteresis_time_out=0.01, nevery=1, hysteresis_time_in=0.01)

def test_hysteresis_fix():
    hysteresis_test(pick_seed_with='fix', hysteresis_time_out=0.01, nevery=1, nsteps=1, fix_nevery=5)
    hysteresis_test(pick_seed_with='fix', hysteresis_time_out=0.001, nevery=1, nsteps=10, fix_nevery=5)
    hysteresis_test(pick_seed_with='fix', hysteresis_time_out=0.001, nevery=5, nsteps=10, fix_nevery=10)
    hysteresis_test(pick_seed_with='fix', hysteresis_time_out=0.05, nevery=10, nsteps=10, fix_nevery=10)

def test_hysteresis_fix_and_init_group():
    hysteresis_test(pick_seed_with='fix_and_init_group', hysteresis_time_out=0.01, nevery=1, nsteps=1, fix_nevery=5)
    hysteresis_test(pick_seed_with='fix_and_init_group', hysteresis_time_out=0.001, nevery=1, nsteps=10, fix_nevery=5)
    hysteresis_test(pick_seed_with='fix_and_init_group', hysteresis_time_out=0.001, nevery=5, nsteps=10, fix_nevery=10)
    hysteresis_test(pick_seed_with='fix_and_init_group', hysteresis_time_out=0.05, nevery=10, nsteps=10, fix_nevery=10)



if __name__ == "__main__":
    test_hysteresis_group()
    test_hysteresis_fix()
    test_hysteresis_fix_and_init_group()