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


import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C, WhiteKernel
import ase.io
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
from utils import read_lammps_dump
import global_plot_settings
global_plot_settings.normal_text()
import argparse


# get path to results files with argparse

parser = argparse.ArgumentParser(description='Plot speedup')
parser.add_argument('--in_path', type=str, help='Path to results files')
parser.add_argument('--out_path', type=str, help='Path to save plots')
parser.add_argument('--n_cells','-n', type=int, help='Number of cells in each dimension')
parser.add_argument('--all_algs', action='store_true', help='Plot all algorithms')
args = parser.parse_args()

path = args.in_path
n = args.n_cells
cheap_times = np.loadtxt(f'{path}/cheap/results.txt')
expensive_times = np.loadtxt(f'{path}/expensive/results.txt')
mixed_times_non_balanced = np.loadtxt(f'{path}/mixed_non_load_balanced/results.txt')
if not args.all_algs:
    try:
        mixed_times_balanced = np.loadtxt(f'{path}/mixed_load_balanced/results.txt')
    except FileNotFoundError:
        mixed_times_balanced = np.loadtxt(f'{path}/mixed_load_balanced/brick/results.txt')
else:
    mixed_times_brick = np.loadtxt(f'{path}/mixed_load_balanced/brick/results.txt')
    mixed_times_tiled_shift = np.loadtxt(f'{path}/mixed_load_balanced/tiled_shift/results.txt')
    mixed_times_tiled_rcb = np.loadtxt(f'{path}/mixed_load_balanced/tiled_rcb/results.txt')

struct = read_lammps_dump(f'{path}/dump{n}.lammpstrj')[-1]
N = len(struct)
i2_potential_1 = struct.arrays['i2_potential[1]']
i2_potential_2 = struct.arrays['i2_potential[2]']
N_1 = len(np.where(i2_potential_1==1)[0])
N_2 = len(np.where(i2_potential_2==1)[0])

# sort each by first column
expensive_times = expensive_times[expensive_times[:,0].argsort()]
mixed_times_non_balanced = mixed_times_non_balanced[mixed_times_non_balanced[:,0].argsort()]
cheap_times = cheap_times[cheap_times[:,0].argsort()]
if not args.all_algs:
    mixed_times_balanced = mixed_times_balanced[mixed_times_balanced[:,0].argsort()]
else:
    mixed_times_brick = mixed_times_brick[mixed_times_brick[:,0].argsort()]
    mixed_times_tiled_shift = mixed_times_tiled_shift[mixed_times_tiled_shift[:,0].argsort()]
    mixed_times_tiled_rcb = mixed_times_tiled_rcb[mixed_times_tiled_rcb[:,0].argsort()]

#for just cheap and expensive, divide by number of timesteps and plot speedup
n_procs = [i for i in range(1,49)]
n_steps = np.array(n_procs)*160

time_per_step_cheap = cheap_times[:,1] / n_steps
time_per_step_expensive = expensive_times[:,1] / n_steps

speedup_cheap = time_per_step_cheap[0] / time_per_step_cheap
speedup_expensive = time_per_step_expensive[0] / time_per_step_expensive

plt.plot(n_procs, speedup_cheap, label='Speedup (cheap)')
plt.plot(n_procs, speedup_expensive, label='Speedup (expensive)')
plt.plot(n_procs, n_procs, linestyle='--', color='r', label='Theoretical speedup')
plt.xlabel('Processor number')
plt.ylabel('Speedup vs serial')
plt.legend()
plt.savefig(f'{args.out_path}/speedup_cheap_vs_expensive.png')

cheap_vs_expensive_speedup = time_per_step_expensive / time_per_step_cheap
theo_speedups = (N * cheap_vs_expensive_speedup) / ((N_1*cheap_vs_expensive_speedup) + N_2)

if not args.all_algs:
    speedup_balanced = expensive_times[:,1] / mixed_times_balanced[:,1]
else:
    speedup_balanced_brick = expensive_times[:,1] / mixed_times_brick[:,1]
    speedup_balanced_tiled_shift = expensive_times[:,1] / mixed_times_tiled_shift[:,1]
    speedup_balanced_tiled_rcb = expensive_times[:,1] / mixed_times_tiled_rcb[:,1]

speedup_non_balanced = expensive_times[:,1] / mixed_times_non_balanced[:,1]
plt.figure()

kernel = C(1.0, (1e-3, 1e3)) * Matern(length_scale=1.0, length_scale_bounds=(1e-2, 1e2), nu=1.5) + WhiteKernel(noise_level=1e-5)

def fit_gp(x, y):
    gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-5, n_restarts_optimizer=10)
    gp.fit(x.reshape(-1, 1), y)
    x_pred = np.linspace(x.min(), x.max(), 100).reshape(-1, 1)
    y_pred, sigma = gp.predict(x_pred, return_std=True)
    return x_pred, y_pred, sigma

if not args.all_algs:
    x_pred_balanced, y_pred_balanced, sigma_balanced = fit_gp(expensive_times[:,0], speedup_balanced)
    plt.plot(x_pred_balanced, y_pred_balanced, color='tab:blue')
    plt.plot(expensive_times[:,0], speedup_balanced, 'o', label='Speedup (load balanced)',color='tab:blue')
    plt.fill_between(x_pred_balanced.ravel(), y_pred_balanced - 2 * sigma_balanced, y_pred_balanced + 2 * sigma_balanced, alpha=0.2,color='tab:blue')
else:
    x_pred_brick, y_pred_brick, sigma_brick = fit_gp(expensive_times[:,0], speedup_balanced_brick)
    x_pred_shift, y_pred_shift, sigma_shift = fit_gp(expensive_times[:,0], speedup_balanced_tiled_shift)
    x_pred_rcb, y_pred_rcb, sigma_rcb = fit_gp(expensive_times[:,0], speedup_balanced_tiled_rcb)
    plt.plot(x_pred_brick, y_pred_brick, label='GP Speedup (load balanced brick)', color='tab:blue')
    plt.fill_between(x_pred_brick.ravel(), y_pred_brick - 2 * sigma_brick, y_pred_brick + 2 * sigma_brick, alpha=0.2, color='tab:blue')
    plt.plot(x_pred_shift, y_pred_shift, label='GP Speedup (load balanced tiled_shift)', color='tab:purple')
    plt.fill_between(x_pred_shift.ravel(), y_pred_shift - 2 * sigma_shift, y_pred_shift + 2 * sigma_shift, alpha=0.2, color='tab:purple')
    plt.plot(x_pred_rcb, y_pred_rcb, label='GP Speedup (load balanced tiled_rcb)', color='tab:cyan')
    plt.fill_between(x_pred_rcb.ravel(), y_pred_rcb - 2 * sigma_rcb, y_pred_rcb + 2 * sigma_rcb, alpha=0.2, color='tab:cyan')

x_pred_non_balanced, y_pred_non_balanced, sigma_non_balanced = fit_gp(expensive_times[:,0], speedup_non_balanced)
plt.plot(x_pred_non_balanced, y_pred_non_balanced, color='tab:orange')
plt.plot(expensive_times[:,0], speedup_non_balanced, 's', label='Speedup (non load balanced)',color='tab:orange')
plt.fill_between(x_pred_non_balanced.ravel(), y_pred_non_balanced - 2 * sigma_non_balanced, y_pred_non_balanced + 2 * sigma_non_balanced, alpha=0.2,color='tab:orange')

x_pred_theo, y_pred_theo, sigma_theo = fit_gp(expensive_times[:,0], theo_speedups)
plt.plot(x_pred_theo, y_pred_theo, linestyle='--', color='tab:green')
plt.plot(expensive_times[:,0], theo_speedups, '^', label='Measured upper bound', color='tab:green')
plt.fill_between(x_pred_theo.ravel(), y_pred_theo - 2 * sigma_theo, y_pred_theo + 2 * sigma_theo, alpha=0.2, color='tab:green')
# plt.plot(expensive_times[:,0], theo_speedups, linestyle='--', color='r', label='Measured upper bound')

plt.xlabel('Processor number')
plt.ylabel('Speedup vs all expensive simulation')
plt.legend()
plt.title(f"{N} atoms")
if not args.all_algs:
    plt.savefig(f'{args.out_path}/speedup_Si_{n}_fixed_atoms.png')
else:
    plt.savefig(f'{args.out_path}/speedup_Si_{n}_fixed_atoms_all_algs.png')