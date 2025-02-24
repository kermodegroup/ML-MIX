# Case Study: Fe Dumbbell Diffusion

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **plots and potentials** from Zenodo at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14920348.svg)](https://doi.org/10.5281/zenodo.14920348)

The full trajectory data is ~30 GB, and is therefore too large for Zenodo. This is available on request, email: fraser.birks@warwick.ac.uk.

## Overview
This study compares three things:
-  The diffusion coefficient of an Fe interstitial Dumbbell in bulk Fe measured in a simulation which uses only an expensive Fe ACE potential
-  The same diffusion coefficient measured with only a cheap Fe ACE potential
-  The same diffusion coefficient measured in an ML/ML simulation where a cheap Fe ACE potential is used for most of the W atoms and the expensive ACE potential is confined only to a small region near the dumbbell.

This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/**/params.py`.

## Running the Case Study (**Zenodo download only**)
To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/[temp]/[simulation_name]
mpirun -np 40 python ../../../scripts/fe_dumbbell_diffusion_with_fix.py
```

Replace `[temp]/[simulation_name]` with the appropriate directory name.

Note that `2_10` refers to the ML/ML simulation, `2_10_only_reference` is the all-cheap simulation and `reference` is the all-expensive simulation.

To analyse the diffusion coefficient data (**if full trajectories have been obtained**) navigate to the top directory and run
```bash
python scripts/find_dumbbell_D_many_single_runs.py T sim_name
```
Where T and sim_name are replaced a given temperature and simulation name.

Finally, plot all simulation results with:
```bash
python scripts/plot_dvals.py
```

## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
