# Ancillary Study: Impact of Constrained Fitting

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14920348.svg)](https://doi.org/10.5281/zenodo.14920348)



## Overview
This study compares the strain error in sheared and relaxed blocks of Si in ML/ML simulations compared to an all-expensive potential reference. In each case, a central sphere of atoms is simulated using an expensive Si ACE potential. For the ML/ML simulations, two different cheap potentials are used for the other atoms in the domain:
-  A cheap potential fit only to high temperature MD data (unconstrained).
-  A cheap potential fit to high temperature MD data and also constrained to have matching elastic constants (constrained).


This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/params.py`.

## Running the Case Study (**Zenodo download only**)
To generate the unrelaxed strained Si structures used in each case, run the following:
```bash
cd simulation_input/generate_structures
python ../../scripts/gen_structs.py
```

To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/[simulation_name]/
mpirun -np 40 python ../../scripts/relax_struct.py
```

Replace `[simulation_name]` with the appropriate directory name.

To analyse the data given in `results/`, navigate to the top directory and run
```bash
python scripts/plot_strain_error.py
```

## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
