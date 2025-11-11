# 
# Copyright (2025) Fraser Birks
#
# This script is licensed under the MIT License. You are free to use, modify,
# and distribute it, provided that this copyright notice and license text 
# remain intact.
#
# See the LICENSE file in the repository root for full details.
# ----------------------------------------------------------------------

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import global_plot_settings
from collections import defaultdict
from scipy.stats import beta

global_plot_settings.normal_text()

ref_label = "Expensive (mpa-0-medium-finetuned)" #-finetuned
plot_label = "Mixed (mpa-0-medium-finetuned, W SNAP)"

print("Usage: python plot_graphs.py <res_path> <label> <plot_ref> <ref_path>")
res_path = sys.argv[1]
pot_name = sys.argv[2]
if len(sys.argv) > 3:
    plot_ref = sys.argv[3]
    if plot_ref not in ['True', 'False']:
        print("plot_ref must be either 'True' or 'False'")
        sys.exit(1)
    plot_ref = (plot_ref == 'True')
else:
    plot_ref = False

if plot_ref:
    ref_stem = sys.argv[4]
### PLOT EXPERIMENTAL DATA ###

# Read data from file
data = np.genfromtxt(f'{script_path}/../results/experimental_res.csv', delimiter=',', dtype=None, encoding=None)

# Organize data by group (3rd column)
groups = defaultdict(dict)
for x, y, group, label in data:
    if label == 'mean':
        groups[group]['x'] = x
        groups[group]['y'] = y
    elif label == 'uppery':
        groups[group]['yhigh'] = y
    elif label == 'lowery':
        groups[group]['ylow'] = y
    elif label == 'upperx':
        groups[group]['xhigh'] = x
    elif label == 'lowerx':
        groups[group]['xlow'] = x

# Prepare arrays for plotting
x = []
y = []
xerr = [[], []]  # [lower, upper]
yerr = [[], []]

for i in range(len(groups)-1):
    g = list(sorted(groups))[i]
    entry = groups[g]
    x_val = entry['x']
    y_val = entry['y']
    
    x_lower = x_val - entry.get('xlow', x_val)
    x_upper = entry.get('xhigh', x_val) - x_val
    y_lower = y_val - entry.get('ylow', y_val)
    y_upper = entry.get('yhigh', y_val) - y_val
    
    x.append(x_val)
    y.append(y_val)
    xerr[0].append(x_lower)
    xerr[1].append(x_upper)
    yerr[0].append(y_lower)
    yerr[1].append(y_upper)

yerr_reflection = [[], []]
yerr_reflection[0] = yerr[1]
yerr_reflection[1] = yerr[0]
####### OLD MD DATA ######
old_md_data = np.loadtxt(f"{script_path}/../results/old_md_data.txt")

####### PLOT SIMULATION DATA #######

results_path = f"{res_path}/results.txt"
data = np.loadtxt(results_path)
#count columns
if data.shape[1] == 5:
    # add a column of 100s to the left hand side of data
    data = np.hstack((np.full((data.shape[0], 1), 100), data))
totals = data[:, 0]
energy = data[:, 1]
frac_pass_through = data[:, 2]
frac_reflected = data[:, 3]
mean_implantation_depth = data[:, 4]
mean_implantation_depth_error = data[:, 5]
min_depths = data[:,6]
max_depths = data[:,7]

if plot_ref:
    ref_path = f"{ref_stem}/results.txt"
    ref_data = np.loadtxt(ref_path)
    #count columns
    if ref_data.shape[1] == 5:
        # add a column of 100s to the left hand side of data
        ref_data = np.hstack((np.full((ref_data.shape[0], 1), 100), ref_data))
    ref_totals = ref_data[:, 0]
    ref_energy = ref_data[:, 1]
    ref_frac_pass_through = ref_data[:, 2]
    ref_frac_reflected = ref_data[:, 3]
    ref_mean_implantation_depth = ref_data[:, 4]
    ref_mean_implantation_depth_error = ref_data[:, 5]
    ref_min_depth = ref_data[:,6]
    ref_max_depth = ref_data[:,7]



fig, ax = plt.subplots(1, 1, figsize=(8, 6))

if plot_ref:
    # Reference dataset
    ax.errorbar(ref_energy, ref_mean_implantation_depth, 
                yerr=ref_mean_implantation_depth_error, 
                fmt='o', color='tab:blue', label=f'{ref_label}')
    
    # Add min-max range for reference
    # for xvals, y_min, y_max in zip(ref_energy, ref_min_depth, ref_max_depth):
    #     ax.vlines(xvals, y_min, y_max, color='k', lw=1)
    #     ax.plot(xvals, y_min, marker='x', color='k', ms=6)
    #     ax.plot(xvals, y_max, marker='x', color='k', ms=6)

    # Your dataset
    plot_color = "tab:orange"

else:
    plot_color = "tab:blue"

ax.errorbar(energy, mean_implantation_depth, 
            yerr=mean_implantation_depth_error, 
            fmt='o', color=plot_color, label=f"{plot_label}")

# for xvals, y_min, y_max in zip(energy, min_depths, max_depths):
#     ax.vlines(xvals, y_min, y_max, color='k', lw=1)
#     ax.plot(xvals, y_min, marker='x', color='k', ms=6)
#     ax.plot(xvals, y_max, marker='x', color='k', ms=6)


ax.legend()
ax.set_xlabel('He energy (eV)')
ax.set_ylabel('Mean implantation depth (A)')
# ax.set_title(f"{pot_name} - Mean implantation depth")

plt.savefig(f"./{pot_name}_implantation_depth.png", dpi=300)

plt.close(fig)
fig, ax = plt.subplots(1, 1, figsize=(8, 6))
if plot_ref:
    ax.plot(ref_energy, ref_frac_pass_through, 'o', color='tab:blue', label=f'{ref_label}')
    ax.plot(energy, frac_pass_through, 'o', color='tab:orange', label=f"{plot_label}")
else:
    ax.plot(energy, frac_pass_through, 'o', color='tab:blue', label=f"{ref_label}")
ax.legend()
ax.set_xlabel('He energy (eV)')
ax.set_ylabel('Fraction of He atoms passing through')
# ax.set_title(f"{pot_name} - Fraction of He atoms passing through")
plt.savefig(f"./{pot_name}_frac_pass_through.png", dpi=300)

plt.close(fig)

fig, ax = plt.subplots(1, 1, figsize=(8, 6))

if plot_ref:
    a_ref = ((ref_frac_reflected * ref_totals) + 1)
    b_ref = (ref_totals - (ref_frac_reflected * ref_totals) + 1)
    mean_ref = a_ref / (a_ref + b_ref)
    ci_lower_ref = beta.ppf(0.025, a_ref, b_ref)
    ci_upper_ref = beta.ppf(0.975, a_ref, b_ref)
    ax.errorbar(ref_energy, mean_ref, yerr=[mean_ref - ci_lower_ref, ci_upper_ref - mean_ref], fmt='o-', color='tab:blue', label=f'{ref_label} (95\% CI)')
a = ((frac_reflected * totals) + 1)
b = (totals - (frac_reflected * totals) + 1)
mean = a / (a + b)
print(a, b)
ci_lower = beta.ppf(0.025, a, b)
ci_upper = beta.ppf(0.975, a, b)
if plot_ref:
    ax.errorbar(energy, mean, yerr=[mean - ci_lower, ci_upper - mean], fmt='o-', color='tab:orange', label=f'{plot_label} (95\% CI)')
else:
    ax.errorbar(energy, mean, yerr=[mean - ci_lower, ci_upper - mean], fmt='o-', color='tab:blue', label=f'{ref_label} (95\% CI)')
ax.plot(old_md_data[:, 0], old_md_data[:, 1], 'o-', color='tab:purple', label='W EAM + pair (Borovikov et al)')
print(x)
print(y)
print(xerr)
ax.errorbar(x, 1-np.array(y), xerr=xerr, yerr=yerr_reflection, fmt='o-', capsize=5, color='tab:green', label='Experimental data (van Gorkum)')
ax.set_xlabel('He energy (eV)')
ax.set_ylabel('Fraction of He atoms reflected')
ax.set_ylim(0.2,1.05)
# ax.set_title(f"{pot_name} - Fraction of He atoms reflected")
ax.legend()
plt.tight_layout()
plt.savefig(f"./{pot_name}_frac_reflected.png", dpi=300)
plt.close(fig)
