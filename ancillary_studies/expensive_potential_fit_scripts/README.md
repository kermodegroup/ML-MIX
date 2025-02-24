# Ancillary Study: Expensive Potential Fitting Scripts

## ⚠️ Important Notice ⚠️
The version of this case study on GitHub is **incomplete**. It is missing required potentials, result data, and plots. You can download the **full study** from Zenodo at:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14920348.svg)](https://doi.org/10.5281/zenodo.14920348)

## Overview
These folders contain the ACEpotentials.jl fitting scripts and datasets for the fitting of the Fe and W expensive potentials. Scripts are provided courtesy of Lakshmi Shenoy (Fe) and Matthew Nutter (W-He). The Fe dataset is from [Zhang et al 2023](https://doi.org/10.1038/s41524-023-01174-6) and the W-He dataset is a relabelled version of the dataset given by [Nutter et al 2024](https://doi.org/10.48550/arXiv.2406.08368). For the Si expensive potential and fitting script we refer the reader to the [ACEpotentials.jl workflow repository](https://github.com/ACEsuit/ACEworkflows/tree/main/Silicon_2023).

## Running the Scripts (**Zenodo download only**)
To fit a specific potential, run

```bash
julia acefit.jl
```

in the corresponding file. To change the number of parallel processors used for the fit, edit the number in the following line of the script:

```python
addprocs(30, exeflags="--project=$(Base.active_project())")
```

## Dependencies
Ensure that you have set up the julia environment correctly, following the instructions in repository README file.

## Citation
If you use this case study, please cite the associated publication and Zenodo dataset.
