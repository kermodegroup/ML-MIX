# Case Study: He implantation into W

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14920348.svg)](https://doi.org/10.5281/zenodo.14920348)

## Overview
This study compares two things:
-  The reflection coefficient for He normally deposited into {100} W surfaces at 1000 K.
-  The average implantation depth for the same conditions.

It contains code to run the study for both the 'all-expensive' case (using a fine-tuned W-He MACE-mpa-0 potential) and the 'ML/ML' case where a W SNAP potential is used for regions away from the implanting He. Note that the scripts are set up to run on GPU, meaning that ML-MIX-KOKKOS must be used, and LAMMPS must be compiled to work with both the MLIAP and symmetrix packages for MACE evaluation. 

This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/params.py`.

## Running the Case Study (**Zenodo download only**)
To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/[simulation_name]/
mpirun -np 40 python ../../scripts/implant_He.py
```

Replace `[simulation_name]` with the appropriate directory name.

To analyse the resulting output data, navigate to the results directory and run
```bash
python ../../scripts/analyse_data.py
```

Finally, plot all simulation in the plots/ directory with
```bash
python ../scripts/plot_graphs.py ../simulation_input/snap_ml_mix_foundation/ compare_plot True ../simulation_input/W_He_MACE_only/
```

## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
