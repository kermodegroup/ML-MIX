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
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import matplotlib.pyplot as plt
import global_plot_settings
global_plot_settings.normal_text()
pot_type = 'espresso_3'
import os

#build figure
fig, ax = plt.subplots(1, 1, figsize=(8, 6))


investigations = ['reference', '2_10_ref', '2_10_mixed']
labels = ['ACE 3_20 only', 'ACE 2_10 only', 'Mixed']
colors = ['tab:blue', 'tab:orange', 'tab:green']
for i, investigation in enumerate(investigations):
    temps = [800, 900, 1000, 1100] #600
    input_paths = [f'plots/{pot_type}_many_single_Ds/{temp}/{investigation}/mean_D.txt' for temp in temps]

    D_vals = []
    errs_D = []
    log_10_D_vals = []
    log_10_errs_D = []

    for path in input_paths:
        D_vals.append(np.loadtxt(path)[0])
        errs_D.append(np.loadtxt(path)[1])
    
    D_vals = np.array(D_vals)
    errs_D = np.array(errs_D)
    frac_errs_D = errs_D/D_vals
    log_10_D_vals = np.log10(D_vals)
    log_10_errs_D = frac_errs_D/np.log(10)
    recip_temps = 1000/np.array(temps)
    

    # plot log D vs 1000/T
    
    # plot data
    ax.errorbar(recip_temps, log_10_D_vals, yerr=log_10_errs_D, fmt='o', color=colors[i], label=labels[i])
    # fit line
    popt, pcov = np.polyfit(recip_temps, log_10_D_vals, 1, w=1/log_10_errs_D, cov=True)
    m, c = popt
    m_err = np.sqrt(pcov[0, 0])
    c_err = np.sqrt(pcov[1, 1])
    # plot fit
    recip_extra_temps = np.array([np.min(recip_temps)-0.1, np.max(recip_temps)+0.1])
    ax.plot(recip_extra_temps, m*recip_extra_temps + c, color=colors[i])

#save fig
plt.xlabel('1000/T(K)')
plt.ylabel(r'log$_{10}$ D ($\mathrm{\AA}^2$/ps)')
plt.title('Diffusion Coefficient of Fe dumbbell')
plt.legend()
plt.tight_layout()

if not os.path.exists(f'plots/{pot_type}_plots'):
    os.makedirs(f'plots/{pot_type}_plots')



plt.savefig(f'plots/{pot_type}_plots/D_vals_diff_temps.png')