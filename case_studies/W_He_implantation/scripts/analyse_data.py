# 
# Copyright (2025) Fraser Birks
#
# This script is licensed under the MIT License. You are free to use, modify,
# and distribute it, provided that this copyright notice and license text 
# remain intact.
#
# See the LICENSE file in the repository root for full details.
# ----------------------------------------------------------------------

import os
import numpy as np
import sys
def find_implantation_files(root_dir):
    """Find all files named 'implantation_depths.txt' in subdirectories of root_dir."""
    files = []
    for dirpath, _, filenames in os.walk(root_dir):
        if 'implantation_depths.txt' in filenames:
            files.append(os.path.join(dirpath, 'implantation_depths.txt'))
    return files

def load_data(file_path):
    """Load data from a single file using NumPy. Each row is [reflection_flag, depth]."""
    try:
        data = np.loadtxt(file_path)
        if data.ndim == 1 and data.size == 3:
            data = data.reshape(1, 3)  # Handle single-line files
        elif data.ndim == 1 and data.size == 2:
            data = data.reshape(1, 2)
        return data
    except Exception as e:
        print(f"Warning: Failed to read {file_path}: {e}")
        return np.empty((0, 3))

def main(root_dir,energy):
    print(f"Analysing {root_dir}")
    files = find_implantation_files(root_dir)
    if not files:
        print("No 'implantation_depths.txt' files found.")
        return


    for i,file in enumerate(files):
        data = load_data(file)
        if i > 0:
            all_data = np.vstack((all_data, data))
        else:
            all_data = data

    if all_data.shape[0] == 0:
        print("No valid data found.")
        return

    mask = ~((all_data[:, 0] == 0.0) & (all_data[:, 1] == 0.0))
    reflections = all_data[mask, 0]
    depths = all_data[mask, 1]
    total = len(reflections)
    pass_through_mask = (depths == -10000)
    frac_pass_through = np.sum(pass_through_mask)/total
    print(f"Num passed through: {np.sum(pass_through_mask)}")
    reflections = reflections[~pass_through_mask]
    depths = depths[~pass_through_mask]
    reflected = np.sum(reflections == 1)
    embedded_depths = depths[reflections == 0]

    frac_reflected = reflected / total
    print(f"Total entries: {total}")
    print(f"Percentage reflected: {100*frac_reflected:.2f}%")

    if embedded_depths.size > 0:
        mean_depth = np.mean(embedded_depths)
        stderr_depth = np.std(embedded_depths, ddof=1) / np.sqrt(embedded_depths.size)
        print(f"Average implantation depth (of embedded): {mean_depth:.3f} ± {stderr_depth:.3f}")
        min_depth = embedded_depths[np.argmin(np.abs(embedded_depths))]
        max_depth = embedded_depths[np.argmax(np.abs(embedded_depths))]
    else:
        print("No embedded particles to compute depth statistics.")
        mean_depth = 0.0
        stderr_depth = 0.0
        min_depth = 0.0
        max_depth = 0.0

    res_arr = np.array([total, float(energy), frac_pass_through, frac_reflected, mean_depth, stderr_depth,min_depth,max_depth])
    return res_arr

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python script.py energy1 energy2 etc")
    else:
        full_arr = []
        for i in range(1,len(sys.argv)):
            energy = sys.argv[i]
            path = f"./{energy}eV"
            res_arr = main(path, energy)
            full_arr.append(res_arr)

        full_arr = np.vstack(full_arr)
        np.savetxt("results.txt", full_arr,header="total, energy (eV), fraction passed through, fraction reflected, mean depth (A), std err depth (A), min depth (A), max depth (A)")
