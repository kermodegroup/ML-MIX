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
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import global_plot_settings
global_plot_settings.normal_text()

data_names = ['W_He_attempt', 'Si_attempt', 'Fe_attempt']
labels = ['W-He', 'Si', 'Fe']
# Initialize plot
plt.figure(figsize=(8, 5))

# Define colors for lines and points
line_colours = ['tab:green', 'tab:orange', 'tab:blue']
point_colours_8K = ['#2ca02c', '#ff7f0e', '#1f77b4']  # Original colors for 8K points
point_colours_250K = ['#237a1b', '#cc6600', '#1a5e8a']  # Slightly darker shades for 250K points

# Loop through each data name and perform GP regression
for i, data_name in enumerate(data_names):
    # Load data from file (assuming first column: x values, second column: y values, third column: std deviations)
    data1 = np.loadtxt(f"all_data/{data_name}_8K_speedups.txt") 
    data2 = np.loadtxt(f"all_data/{data_name}_250K_speedups.txt")
    full_data = np.vstack((data1, data2))

    X = full_data[:,0].reshape(-1,1)  # Feature values from 1 to 48
    y = full_data[:, 1]  # Target values
    y_err = full_data[:, 2]  # Standard deviation of each point

    # Define the GP kernel: Matérn kernel with length scale 1.0 + a constant + noise term
    kernel = C(1.0, (1e-3, 1e3)) * Matern(length_scale=1.0, length_scale_bounds=(1e-2, 1e2), nu=1.5) + WhiteKernel(noise_level=1e-5)

    # Fit the Gaussian Process
    gp = GaussianProcessRegressor(kernel=kernel, alpha=y_err**2, n_restarts_optimizer=10)
    gp.fit(X, y)

    # Make predictions on a finer grid
    X_pred = np.linspace(X.min(), X.max(), 100).reshape(-1, 1)
    y_pred, sigma = gp.predict(X_pred, return_std=True)

    # Plot results
    plt.errorbar(data1[:, 0], data1[:, 1], yerr=data1[:, 2], fmt='o', label=f"{labels[i]}, 8k atoms", markersize=5, capsize=3, color=point_colours_8K[i], ecolor=point_colours_8K[i])
    plt.errorbar(data2[:, 0], data2[:, 1], yerr=data2[:, 2], fmt='s', label=f"{labels[i]}, 250k atoms", markersize=5, capsize=3, color=point_colours_250K[i], ecolor=point_colours_250K[i])
    plt.plot(X_pred, y_pred, color=line_colours[i])
    plt.fill_between(X_pred.ravel(), y_pred - 2 * sigma, y_pred + 2 * sigma, color=line_colours[i], alpha=0.2)

plt.xlabel('Processor number')
plt.ylabel('Speedup factor')
plt.title('Speedup of all cheap vs all expensive simulation')
plt.legend()
plt.ylim(1,90)
plt.savefig(f'../plots/gp_regression_final.png')
