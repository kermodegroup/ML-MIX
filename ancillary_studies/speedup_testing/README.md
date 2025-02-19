# Ancillary Study: Speedup Profiling Tests

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[Zenodo DOI Link]

## Overview
This study profiles the following:

- All-cheap and all-expensive simulations in Si, Fe and W on between 1 and 48 cores. Three sizes of atomic domains are used for each. For Fe and W, 16x16x16, 32x32x32 and 50x50x50 unit cell domains are used. For Si, 10x10x10, 20x20x20 and 32x32x32 domains are used. For each processor number on each domain, there are three trials. 

- All-cheap and all-expensive simulations in Fe, Si and W on 27 cores for many different domain sizes, up to $10^{6}$ atoms. 

- ML/ML simulations in Si for different numbers of processors on 10x10x10, 20x20x20 and 32x32x32 domains, with different load balancing strategies. Load balancing strategies tested are: non-load balanced, brick (shift), tiled shift and tiled rcb.

- ML/ML simulations in Si on 27 processors but an increasing number of atoms, up to $10^{6}$, both load-balanced (using brick (shift)) and non-load balanced.


This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/.../params.py`.

## Running the Case Study (**Zenodo download only**)
To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/.../
mpirun -np 40 python .../scripts/speed_test_MD.py
```

Replace `/.../` with the appropriate directory name.

Before running the simulations in fixed_proc_increasing_atoms, run 
```bash
python build_structs_[el].py
```
in the corresponding `[el]_attempt/fixed_proc_increasing_atoms` folder. This is to build all of the different sizes of domains used in the simulation.

The LAMMPS output log files should be saved to `.out` files in each directory. The `extract_results.py` script can extract the relevant timings from these output fies. 

To collect and plot data from all of the all-cheap and all-expensive speedup simulations given in `/results`, run:
```bash
python /scripts/collect_all_data.py
python scripts/fit_gp_to_speedup_data.py
```
in the top directory.

Finally, to plot all the figures for each ML/ML profiling investigation, run the following:
```bash
python scripts/plot_speedup_with_GP_Si_fixed_atoms.py --in_path results/Si_attempt/fixed_atoms_increasing_prc/n\=10/trial_1/ --out_path plots/ -n 10 --all_algs
```
```bash
python scripts/plot_speedup_with_GP_Si_fixed_proc.py --in_path results/Si_attempt/fixed_proc_increasing_atoms/ --out_path plots/
```
Directories are passed in via --in_path and --out_path, and for the fixed atoms, increasing processor investigations the number of unit cell should also be passed in via -n. Specifying the --all_algs option allows a comparison between different load-balancing algorithms.

## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
