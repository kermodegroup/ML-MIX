# Case Study: He Diffusion in W

## ⚠️ Important Notice
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[Zenodo DOI Link]

## Overview
This study compares the diffusion coefficient measured from an simulation using only an expensive W-He ACE potential and the value measured from an ML/ML simulation where a pure W uf3 potential is used for most of the W atoms and the ACE potential is confined only to a small region near the He as it diffuses through the structure.

This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/params.py`.

## Running the Case Study
To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/[simulation_name]/
mpirun -np 40 python ../../scripts/He_W_MD_with_fix.py
```

Replace `[simulation_name]` with the appropriate directory name.

**Zenodo download only** To analyse the diffusion coefficient data, navigate to the top directory and run
```bash
python scripts/find_He_W_D_many_single_runs.py
```
To select which set of simulation results to analyse, edit the following lines of the `params.py` file in the top directory:
```python
investigation = 'ref_D_coeff'
T = 800
```
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
