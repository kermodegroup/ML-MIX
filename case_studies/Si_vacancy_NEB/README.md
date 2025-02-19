# Case Study: Si Vacancy NEB

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[Zenodo DOI Link]

## Overview
This study compares:
-  The energy barrier for vacancy migration in Si found using a NEB using only an expensive Si ACE potential. 
-  The energy barrier for vacancy migration in Si found using an ML/ML NEB  where the majority of atoms in the simulation are modelled using a cheaper Si ACE potential.

This case study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/params.py`.

## Running the Case Study **Zenodo download only**
To generate the initial and final images, run the following:
```bash
cd simulation_input
python ../scripts/create_images.py 
```
To then generate an input script to run the LAMMPS NEB, run:
```bash
cd [simulation_name]
python ../../scripts/gen_lammps_input_script.py
```
Replace `[simulation_name]` with the appropriate directory name.

Finally, to run the LAMMPS NEB simulation, e.g with 9 images on 4 cores each, run:

```bash
mpirun -np 36 lmp -partition 9x4 -in neb_input.in
```

To get the final energies of each image in the NEB relaxations, run the following script in the `/simulation_input` directory:

```bash
python ../script/eval_final_images.py
```

Note that this script will need to be altered if the number of images used isn't 9.

To analyse the results provided in `results/`, navigate to this top directory and run

```bash
python scripts/plot_neb.py
```


## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
