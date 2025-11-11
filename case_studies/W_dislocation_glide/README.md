# Case Study: He Diffusion in W

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14920348.svg)](https://doi.org/10.5281/zenodo.14920348)

## Overview
This study was carried out by Matthew Nutter (m.nutter@warwick.ac.uk). Please direct all queries to him. This study compares two things:
-  The velocity of screw dislocations W measured for a variety of temperatures and shear stresses modelled with a 3, 21 ACE potential.
-  The same velocities at the same temperatures and pressures in an ML/ML simulation where the expensive 3, 21 ACE potential is confined to a region around the gliding dislocation and the remainder of the domain is simulated with a 3, 14 ACE potential.

Simple LAMMPS input scripts for running this case study can be found in `simulation_input/`.

## Running the Case Study (**Zenodo download only**)
To run a specific simulation, navigate to the relevant folder and execute the script using `mpirun`. For example:

```bash
cd simulation_input/[simulation_name]/
mpirun -np 40 lmp -in run_dislocation_glide.lammps -var temperature $1 -var stress $2 -var seed $3 -var timestep $4 -var n_timesteps $5
```

Replace `[simulation_name]` with the appropriate directory name.

Dislocation velocities can be read as one of the columns in the LAMMPS output.

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
