import numpy as np
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))
sys.path.append(f'{script_path}/../../../utils')
import matplotlib.pyplot as plt
import global_plot_settings
global_plot_settings.normal_text()

#build figure
fig, ax = plt.subplots(1, 1, figsize=(8, 6))


investigations = ['ref_D_coeff', 'mix_D_coeff']
labels = ['ACE only', 'ACE/UF3 mixed']
colors = ['tab:blue', 'tab:orange']
for i, investigation in enumerate(investigations):
    temps = [400, 600, 800]
    input_paths = [f'plots/{investigation}/{temp}K/average_diffusion_coefficient.txt' for temp in temps]

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
    recip_extra_temps = np.array([1.0, 2.75])
    ax.plot(recip_extra_temps, m*recip_extra_temps + c, color=colors[i])

#save fig
plt.xlabel('1000/T(K)')
plt.ylabel(r'log$_{10}$ D ($\mathrm{\AA}^2$/ps)')
plt.title('Diffusion Coefficient of He in W')
plt.legend()
plt.tight_layout()
plt.savefig(f'plots/D_vals_diff_temps.png')