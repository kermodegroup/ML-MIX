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
from tqdm import tqdm
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import global_plot_settings
global_plot_settings.normal_text()
import argparse
from utils import read_lammps_dump

parser = argparse.ArgumentParser(description='Plot speedup')
parser.add_argument('--in_path', type=str, help='Path to results files')
parser.add_argument('--out_path', type=str, help='Path to save plots')
args = parser.parse_args()

in_path = args.in_path
out_path = args.out_path
# Load data
expensive_times = np.loadtxt(f'{in_path}/expensive/results.txt')
mixed_times_balanced = np.loadtxt(f'{in_path}/mixed_load_balanced/results.txt')
mixed_times_non_balanced = np.loadtxt(f'{in_path}/mixed_non_load_balanced/results.txt')
cheap_times = np.loadtxt(f'{in_path}/cheap/results.txt')

# Sort by first column
expensive_times = expensive_times[expensive_times[:,0].argsort()]
mixed_times_balanced = mixed_times_balanced[mixed_times_balanced[:,0].argsort()]
mixed_times_non_balanced = mixed_times_non_balanced[mixed_times_non_balanced[:,0].argsort()]
cheap_times = cheap_times[cheap_times[:,0].argsort()]

# Generate theoretical speedup
n_vals = cheap_times[:,0].astype(int)
theo_speedups = []
for i, n in tqdm(enumerate(n_vals)):
    struct = read_lammps_dump(f'{in_path}/get_theo_speedup/dump{n}.lammpstrj')[-1]
    N = len(struct)
    N_1 = np.count_nonzero(struct.arrays['i2_potential[1]'] == 1)
    N_2 = np.count_nonzero(struct.arrays['i2_potential[2]'] == 1)
    expensive_cheap_speedup = expensive_times[i,1] / cheap_times[i,1]
    theo_speedup = (N * expensive_cheap_speedup) / ((N_1 * expensive_cheap_speedup) + N_2)
    theo_speedups.append(theo_speedup)

theo_speedups = np.array(theo_speedups)

# Compute speedups
speedup_balanced = expensive_times[:,1] / mixed_times_balanced[:,1]
speedup_non_balanced = expensive_times[:,1] / mixed_times_non_balanced[:,1]
natoms = (expensive_times[:,0]**3) * 8  # Number of atoms

# Define GP kernel (Matern for flexibility)
kernel = C(1.0, (1e-3, 1e3)) * Matern(length_scale=1.0, nu=1.5) + WhiteKernel(noise_level=1e-5)

# Fit GP function
def fit_gp(x, y):
    x_log, y_log = np.log10(x), np.log10(y)  # Log transform for power-law behavior
    gp = GaussianProcessRegressor(kernel=kernel, alpha=1e-5, n_restarts_optimizer=10)
    gp.fit(x_log.reshape(-1, 1), y_log)
    x_pred_log = np.linspace(x_log.min(), x_log.max(), 100).reshape(-1, 1)
    y_pred_log, sigma_log = gp.predict(x_pred_log, return_std=True)
    return 10**x_pred_log, 10**y_pred_log, 10**(y_pred_log - 2 * sigma_log), 10**(y_pred_log + 2 * sigma_log)

# Fit GPs to data
x_pred_bal, y_pred_bal, y_lower_bal, y_upper_bal = fit_gp(natoms, speedup_balanced)
x_pred_non_bal, y_pred_non_bal, y_lower_non_bal, y_upper_non_bal = fit_gp(natoms, speedup_non_balanced)
x_pred_theo, y_pred_theo, y_lower_theo, y_upper_theo = fit_gp(natoms, theo_speedups)

# Plot results
plt.figure()
plt.plot(natoms, speedup_balanced, 'o', label='Speedup (load balanced)',color='tab:blue')
plt.plot(natoms, speedup_non_balanced, 's', label='Speedup (non load balanced)',color='tab:orange')
plt.plot(natoms, theo_speedups, '^', label='Measured upper bound', color='tab:green')

# GP fits with uncertainty
plt.plot(x_pred_bal, y_pred_bal,color='tab:blue')
plt.fill_between(x_pred_bal.ravel(), y_lower_bal, y_upper_bal, alpha=0.2,color='tab:blue')

plt.plot(x_pred_non_bal, y_pred_non_bal,color='tab:orange')
plt.fill_between(x_pred_non_bal.ravel(), y_lower_non_bal, y_upper_non_bal, alpha=0.2,color='tab:orange')

plt.plot(x_pred_theo, y_pred_theo, linestyle='--', color='tab:green')
plt.fill_between(x_pred_theo.ravel(), y_lower_theo, y_upper_theo, alpha=0.2,color='tab:green')

# Log-log plot settings
plt.xscale('log')
plt.xlabel('Number of atoms')
plt.ylabel('Speedup vs all expensive simulation')
plt.legend()
plt.savefig(f'{out_path}/speedup_Si_fixed_proc.png')
