# Ancillary Study: Energy Conservation Testing

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[Zenodo DOI Link]

## Overview
This study measures the energy drift over time of an NVE simulation of a block of silicon containing a stretched and rigidly fixed bond. Atoms immediately around the central bond are modelled using an expensive Si ACE potential, whilst atoms further from the bond are modelled using cheaper Si ACE potentials.

Three different cheap potentials are tested, and for each there are two different studies - one in which abrupt force mixing is used (no blending region) and the size of the core region is varied, and one where the core + blending radius is fixed at 10 Angstrom, and the blending radius is varied. 

This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/**/params.py`.

## Running the Case Study (**Zenodo download only**)
To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/[cheap_potential]/[simulation_name]/
mpirun -np 40 python ../../../scripts/fixed_bond_MD.py
```

Replace `[cheap_potential]/[simulation_name]` with the appropriate directory name.

To find the energies of each snapshot in each simulation, run the following in the `simulation_input/` directory:

```bash
python ../scripts/find_energies.py
```
Note that this script assumes the results of each simulation are saved in subfolders called `r=4`, `r=6` (matching the file structure in `results/`).

The simulation results (given in `results/`) can be plotted by running the following in the top level directory

```bash
python scripts/plot_energy_flux_potential_type.py
python scripts/plot_energy_flux_smooth_mixing.py
python scripts/plot_pots_seperately_no_normalisation.py
python scripts/plot_pots_seperately_normalised.py
```

## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
