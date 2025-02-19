# Case Study: Si Stretched Bond Average Force

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[Zenodo DOI Link]

## Overview
This study compares:
-  The average force in a stretched and rigidly fixed bond in a block of Si in a simulation which uses only an expensive Si ACE potential
-  The average force in a stretched and rigidly fixed bond in a block of Si in ML/ML simulations where the size of the expensive potential region around the stretched bond is varied.

This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/params.py`.

## Running the Case Study (**Zenodo download only**)
To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/[simulation_name]/
mpirun -np 40 python ../../scripts/fixed_bond_MD.py
```

Replace `[simulation_name]` with the appropriate directory name.

To analyse the stretched bond data, navigate to the top directory and run
```bash
python scripts/find_avg_bond_forces.py
```
Finally, plot all simulation results with:
```bash
python scripts/plot_pots_seperately.py
```

## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
