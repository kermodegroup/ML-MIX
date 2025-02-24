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


dump_name_list = parameter('dump_name_list')
dump_step = parameter('dump_step')
min_sample_number = parameter('min_sample_number')
min_sample_rate = parameter('min_sample_rate')
cell_size = parameter('n')
output_path = parameter('output_path', '.')

D_vals = []
D_errs = []
for p, d_name in enumerate(dump_name_list):
    sqr_displacements = {}

    # first, get the location of the defect every single dump

    #get the coordination and use it to assign probabilities of atoms being
    #in a defect
    traj = ase.io.read(d_name, index=':')
    big_arr = np.zeros((len(traj[0]), len(traj)))
    for j,t in enumerate(traj):
        t.pbc = [True, True, True]
        i = neighbour_list('i', t,  3.3)
        coord = np.bincount(i)
        big_arr[:, j][(coord==11) | (coord==12)] = 1
        big_arr[:,j][coord==13] = 1/2


    #now, take the rolling average over 3 columns at a time
    #also ensure that defect position never crosses PBC by shifting it to
    #cell centre
    rolling_avg = np.zeros_like(big_arr)
    defect_pos = np.zeros((len(traj),3))
    V = np.zeros(3)
    cell_centre = np.array([t.cell[0,0]/2, t.cell[1,1]/2, t.cell[2,2]/2])
    shifted_ts = []
    for i in range(1, len(traj)-1):
        traj[i].set_pbc([True, True, True])
        traj[i].positions += V
        traj[i].wrap()
        rolling_avg[:,i] = np.mean(big_arr[:,i-1:i+2], axis=1)
        #take the top two values in each
        top_two = np.argsort(rolling_avg[:,i])[-2:]
        positions = traj[i].get_positions()
        pos1 = positions[top_two[0]] - V
        pos2 = positions[top_two[1]] - V
        defect_pos[i,:] = (pos1 + pos2)/2
        #now add the difference between the position
        V = cell_centre - defect_pos[i,:] 

    #get rid of the first and last columns of defect_pos
    defect_pos = defect_pos[1:-1,:]


    #now get the min sample rate as 2x the average jump spacing
    diffs = defect_pos[1:] - defect_pos[:-1]
    norm_diff = np.linalg.norm(diffs, axis=1)
    #find indices indicating jump positions
    jump_pos = np.where(norm_diff > 1)[0]
    #get average wait time between jumps
    print('Total number of jumps,', len(jump_pos))
    jump_diff = np.diff(jump_pos)
    mean_jump_wait = np.mean(jump_diff)
    print('Mean jump wait = ', mean_jump_wait, 'frames, which is', (mean_jump_wait*dump_step)/1000, 'ps')
    #min_sample_rate = int(2*mean_jump_wait)
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
    maxlim = max(defect_pos.max(), t.cell[0,0])
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

    np.savetxt(f'{output_path}/mean_squared_displacement.txt', arr_for_file_writing)

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
    plt.title(f'MSD vs timelag for {cell_size}x{cell_size}x{cell_size} 800K Fe Dumbbell')
    plt.tight_layout()
    plt.savefig(f'{output_path}/mean_squared_displacement.png')


    D_vals.append(D)
    D_errs.append(D_err)

#save to file

np.savetxt(f'{output_path}/D_vals.txt', np.array([D_vals, D_errs]).T, header='D ± error')

#save mean and std in mean to file
mean_D = np.mean(D_vals)
#standard error from std
err_D = np.std(D_vals)/np.sqrt(len(D_vals))
np.savetxt(f'{output_path}/mean_D.txt', np.array([mean_D, err_D]), header='mean D ± sigma')