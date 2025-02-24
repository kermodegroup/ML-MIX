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

import ase.io
import sys
import os
script_path = os.path.dirname(os.path.abspath(__file__))

sys.path.append('.')
sys.path.append(f'{script_path}/../../../utils')

from matscipy import parameter
import params
from matscipy.neighbours import neighbour_list
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import global_plot_settings

global_plot_settings.normal_text()


dump_name_list = parameter('dump_name_list')
dump_step = parameter('dump_freq')
min_sample_number = parameter('min_sample_number')
cell_size = parameter('n')
timestep = parameter('tstep',0.001) #ps
min_sample_rate = parameter('min_sample_rate')
output_path = parameter('output_path', '.')


D_vals = []
#first read in initial structure and get ids of all heliums
for p, d_name in enumerate(dump_name_list):
    sqr_displacements = {}
    # first, get the location of the defect every single dump
    traj = ase.io.read(d_name, index=':')
    helium_id = np.where(traj[0].arrays['numbers']==2)[0] #He is type 2

    #now, take the rolling average over 3 columns at a time
    #also ensure that defect position never crosses PBC by shifting it to
    #cell centre
    defect_pos = np.zeros((len(traj),3))
    V = np.zeros(3)
    cell_centre = np.array([traj[0].cell[0,0]/2, traj[0].cell[1,1]/2, traj[0].cell[2,2]/2])
    shifted_ts = []
    for i in range(1, len(traj)-1):
        traj[i].set_pbc([True, True, True])
        traj[i].positions += V
        traj[i].wrap()
        positions = traj[i].get_positions()
        pos = positions[helium_id] - V
        defect_pos[i,:] = pos
        V = cell_centre - defect_pos[i,:]

    # #get rid of the first and last columns of defect_pos
    defect_pos = defect_pos[1:-1,:]

    print('Min sample rate = ', min_sample_rate)
    n_defect_pos = len(traj)-2
    n_samples = int(n_defect_pos/min_sample_number)+1
    #now, build arrays of the square displacement for each n
    for n in range(min_sample_rate,n_samples):
        #if n is not a key in the dictionary, add it
        if n not in sqr_displacements.keys():
            sqr_displacements[n] = np.array([])
        #sum over all time lags
        print('arr_1_index',np.arange(n_defect_pos)[n::n])
        print('arr_2_index',np.arange(n_defect_pos)[:n_defect_pos-n:n])
        shifted_defect_pos = defect_pos[n::n,:]
        reduced_defect_pos = defect_pos[:n_defect_pos-n:n,:]
        sqr_dist = np.linalg.norm(shifted_defect_pos - reduced_defect_pos, axis=1)**2
        sqr_displacements[n] = np.concatenate((sqr_displacements[n],sqr_dist))


    #now plot out the positions
    minlim = min(defect_pos.min(), 0)
    maxlim = max(defect_pos.max(), traj[0].cell[0,0])
    plt.figure()
    plt.plot(defect_pos[:,0],defect_pos[:,1])
    plt.xlabel('xpos')
    plt.ylabel('ypos')
    plt.xlim(minlim, maxlim)
    plt.ylim(minlim,maxlim)
    plt.savefig(f'{output_path}/defect_positions_xy_{p}.png')

    plt.figure()
    plt.plot(defect_pos[:,0],defect_pos[:,2])
    plt.xlabel('xpos')
    plt.ylabel('zpos')
    plt.xlim(minlim, maxlim)
    plt.ylim(minlim,maxlim)
    plt.savefig(f'{output_path}/defect_positions_xz_{p}.png')

    plt.figure()
    plt.plot(defect_pos[:,1],defect_pos[:,2])
    plt.xlabel('ypos')
    plt.ylabel('zpos')
    plt.xlim(minlim, maxlim)
    plt.ylim(minlim,maxlim)
    plt.savefig(f'{output_path}/defect_positions_yz_{p}.png')



    #now, get the mean square displacement for each n
    mean_squared_displacement = {}
    errs = {}
    arr_for_file_writing = np.zeros((len(sqr_displacements.keys()),3))
    for i,n in enumerate(sqr_displacements.keys()):
        #prune off any 0s from the front of the array
        mean_squared_displacement[n] = np.mean(sqr_displacements[n])
        errs[n] = np.std(sqr_displacements[n])/np.sqrt(len(sqr_displacements[n]))
        print('stepsize = ', n, 'mean squared displacement = ', mean_squared_displacement[n], 'error = ', errs[n])
        arr_for_file_writing[i,0] = n
        arr_for_file_writing[i,1] = mean_squared_displacement[n]
        arr_for_file_writing[i,2] = errs[n]

    np.savetxt(f'{output_path}/mean_squared_displacement_{p}.txt', arr_for_file_writing)

    plt.figure()
    t_vals = np.array(list(mean_squared_displacement.keys()))*(dump_step/1000)
    sorted_ts = np.argsort(t_vals)
    t_vals = t_vals[sorted_ts]
    mean_squared_displacement_vals = np.array(list(mean_squared_displacement.values()))[sorted_ts]
    errs_vals = np.array(list(errs.values()))[sorted_ts]

    #now perform a fit to the data

    # Define the linear function
    def linear_func(x, m, c):
        return m * x + c

    # Perform the curve fit
    popt, pcov = opt.curve_fit(linear_func, t_vals, mean_squared_displacement_vals, sigma=errs_vals, absolute_sigma=True)

    # popt contains the best-fit parameters [m, c]
    # pcov is the covariance matrix
    m, c = popt
    m_err = np.sqrt(pcov[0, 0])
    c_err = np.sqrt(pcov[1, 1])

    print(f"Slope (m) = {m:.4f} ± {m_err:.4f}")
    print(f"Intercept (c) = {c:.4f} ± {c_err:.4f}")

    # The diffusion coefficient is 1/6 of the slope
    D = m / 6
    D_err = m_err / 6
    print(f"Diffusion coefficient = {D:.2f} ± {D_err:.2f} Å²/ps")

    plt.errorbar(t_vals, mean_squared_displacement_vals, yerr=errs_vals, fmt='o-',capsize=3)
    # Plot the best-fit line
    plt.plot(t_vals, linear_func(t_vals, *popt), label=f"y = {m:.2f}x + {c:.2f}")
    plt.legend()
    #also write D on the plot in the top left corner
    plt.text(0.5, 0.8, f'D = {D:.4f} ± {D_err:.4f} Å²/ps', fontsize=12,  transform=plt.gca().transAxes)
    plt.xlabel('Time-lag (ps)')
    plt.ylabel('Mean squared displacement')
    plt.title(f'MSD vs timelag for {cell_size}x{cell_size}x{cell_size} 800K W-He')
    #plt.grid()
    plt.tight_layout()
    plt.savefig(f'{output_path}/mean_squared_displacement_{p}.png')

    D_vals.append(D)

# at the end, get the average diffusion coefficient and error in mean
D_vals = np.array(D_vals)
D_mean = np.mean(D_vals)
D_err = np.std(D_vals)/np.sqrt(len(D_vals))
print(f'Average diffusion coefficient = {D_mean:.4f} ± {D_err:.4f} Å²/ps')
#save the diffusion coefficients to a file
np.savetxt(f'{output_path}/diffusion_coefficients.txt', D_vals)

#save the average diffusion coefficient and error to a file
np.savetxt(f'{output_path}/average_diffusion_coefficient.txt', [D_mean, D_err])






