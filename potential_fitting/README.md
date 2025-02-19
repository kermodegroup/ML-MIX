# Constrained Cheap Potential Fitting Examples

## ⚠️ Important Notice ⚠️
The version of this on GitHub is **incomplete**. It is missing some of the required expensive potentials. You can download the **complete version** from Zenodo at:

[Zenodo DOI Link]

## Overview
The examples in this folder walk through the process of generating constrained cheap potentials to match an expensive potential in a domain of interest. In the Fe and Si examples, the cheap models are small ACE potentials. In the W example, the cheap model is an ultra-fast UF3 potential. 

In each case, the optimal model parameters $\mathbf{c}$ are found as a solution to

$$
\min_{\mathbf{c}}\left(||\mathbf{y}_ {\mathrm{S}} - \mathbf{A}_ {\mathrm{S}}\mathbf{c}||^{2} + \lambda ||\mathbf{y}_ {\mathrm{H}} - \mathbf{A}_ {\mathrm{H}} \mathbf{c}||^{2}\right).
$$

Where $\mathbf{A}_ \mathrm{H} \in \mathbb{R}^{N_ {\mathrm{H}} \times N_ {\mathrm{D}}}$ and $\mathbf{A}_ \mathrm{S} \in \mathbb{R}^{N_ {\mathrm{S}} \times N_ {\mathrm{D}}}$ are the design matrices for the hard and soft constraints respectively, and $\mathbf{y}_ \mathrm{H} \in \mathbb{R}^{N_ {\mathrm{H}}}$ and $\mathbf{y}_ \mathrm{S} \in \mathbb{R}^{N_ {\mathrm{S}}}$ are the corresponding expensive potential target vectors. Obtaining a solution to the problem requires finding the minimum $\lambda$ such that the solution to this equation lies on the surface of the constraint subspace

$$
||\mathbf{A}_ {\mathrm{H}} \mathbf{c} - \mathbf{y}_ {\mathrm{H}}||^{2} = \alpha.
$$

In each example given here, $\mathbf{A}_\mathrm{H}$ data corresponds to constrained homogeneous lattice deformations, and $\mathbf{A}_\mathrm{S}$ data is from high temperature bulk MD, but this need not necessarily be the case.

This study follows the **matscipy format**:
- **Scripts for running case studies and analysis** are located in `scripts/`.
- **Input parameters for individual simulations** are located in `simulation_input/**/params.py`.

## Peforming the constrained fitting (**Zenodo download only**)
To generate the $\mathbf{A}_\mathrm{H}$ and $\mathbf{A}_\mathrm{S}$ data for each specific example, run e.g:

```bash
cd examples/[simulation_name]/
python ../../scripts/generate_Ah.py
mpirun -np 40 python ../../scripts/generate_As.py
```
Replace `[simulation_name]` with the appropriate directory name. Note that the generate_As.py script runs LAMMPS MD, so should be run MPI parallel for faster evaluation. 


Once this is done, to fit each individual constrained and unconstrained cheap potential (and evaluate its elastic constants), run the following in each directory:
```bash
chmod u+x fit_potential.sh
./fit_potential.sh
```

## Dependencies
Ensure that you have the required dependencies installed before running the simulations. Install python dependencies via:

```bash
pip install -r requirements.txt
```

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
