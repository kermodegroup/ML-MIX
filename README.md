# ML-MIX ![CI](https://github.com/kermodegroup/ML-MIX/actions/workflows/ci.yml/badge.svg) ![DOCS](https://github.com/kermodegroup/ML-MIX/actions/workflows/docs.yml/badge.svg) 
A LAMMPS plugin for efficiently mixing together forces from different machine-learned interatomic potentials inspired by quantum mechanics/molecular mechanics (QM/MM) methods. Use expensive, complex potentials only where they're needed! 


<div style="display: flex; gap: 10px;">
  <img src="docs/images/blue_background_title_looped.gif" width="500">
  <img src="docs/images/Fe_dumbbell.gif" width="300">
</div>

## Key Information

> **ML-MIX is under active development**  
> Major updates and collections of bugfixes will be periodically merged from the `develop` branch into `main` and released as new versions.
> If you encounter any bugs or unexpected behaviour, please [open an issue](https://github.com/kermodegroup/ML-MIX/issues) or contact us directly via email at [fraser.birks@warwick.ac.uk](fraser.birks@warwick.ac.uk).

Full documentation can be found here: https://kermodegroup.github.io/ML-MIX/

Detailed setup instructions can be found below.

For a user wanting to start quickly, simple examples showcasing and explaining features of ML-MIX can be found in `simple_examples/`

The accompanying paper on ML-MIX can be found here: [arXiv:2502.19081](https://arxiv.org/abs/2502.19081)

All the code associated with this paper is included in this repository, in `potential_fitting/`, `case_studies/` and `ancillary_studies/`, however the GitHub version is incomplete, as it is missing potentials, data and plots. The full version of this repository can be found on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14920348.svg)](https://doi.org/10.5281/zenodo.14920348)

More information can be found in the `README` files inside each subdirectory.

## List of tested compatible pair_styles
ML-MIX is designed to be a wrapper that is compatible with any generic LAMMPS pair_style (so long as that pair_style is compatible with pair_hybrid). However, for one reason or another (and particularly if a pair_style does something that is very 'non-standard') incompatibilities can arise (which may take the form of segmentation faults). 

A common reason that a pair_style may be incompatible is if it involves an MPI communication *during compute()* (e.g, `EAM`). This MPI communication acts as a synchronisation point - meaning that all processors must halt and wait. This can *severely* bottle-neck runtimes and make doing an ML/ML simulation essentially pointless.

### CPU pair_styles

List of pair_styles *tested* and found compatible (non kokkos version):
- `lj/cut`
- `pace`
- `uf3`
- `snap`
- `symmetrix/mace` (in `no_mpi_message_passing` mode)
- `table`

List of pair_styles *tested* and currently incompatible:
- `EAM` and all derivatives
- `symmetrix/mace` (in `mpi_message_passing` and `no_domain_decomposition` mode)
- `allegro`
- `mliappy mace`

### GPU pair_styles

List of pair_styles *tested* and found compatible (kokkos versions):
- `lj/cut/kk`
- `symmetrix/mace/kk` (in `no_mpi_message_passing` mode)
- `snap/kk`

List of pair_styles *tested* and currently incompatible:
- `uf3/kk`
- `symmetrix/mace` (in `mpi_message_passing` and `no_domain_decomposition` mode)


There are likely to be many other compatible pair_styles which have not yet been tested. If you find them, please feel free to get in touch with me (fraser.birks@warwick.ac.uk) or make a PR to update the list! It is also possible that the *currently incompatible* pair_styles can be made compatible with minor tweaks (for example, a `no_mpi` mode could be made for EAM, which would just involve adding additional ghost atoms). If you wish to do this (because you desperately need to use one of the currently incompatible pair_styles), then please get in contact with me or raise a github issue. 

If you believe that a pair_style should be compatible (e.g, because there are no MPI communications involved in *compute*) and think the issue may lie with ML-MIX, then please make a Github issue.


## ⚠️ Bug Warning ⚠️ - **Fixed**
Update: This issue has been resolved since the [17/03/2025 LAMMPS Stable release](https://github.com/lammps/lammps/releases/tag/stable_29Aug2024_update2).

In versions of LAMMPS before 17/03/2025, if a pair_style tried to build a half neighborlist by pruning a full neighborlist that was built for a fix, it lead to all the forces computed by that pair_style being 0 (atoms had no neighbors). This is a problem when using ML-MIX if:
- fix mlml is defined (which needs a full neighborlist)
- *ALL* other pair_styles defined need half neighborlists (if there is even one pair_style defined which needs a full neighborlist, i.e, ACE, UF3, this isn't a problem as then half neighborlists are constructed correctly.)
- To resolve this, please update to the latest version of LAMMPS.


## Setup

### Cloning this repository

To clone this repository (and it's submodule dependencies), run
```
git clone --recurse-submodules https://github.com/kermodegroup/ML-MIX.git
```

### Building ML-MIX

There are two ways to build ML-MIX
1. In-source (required for KOKKOS integration)
2. As a plugin

Instructions for both types of builds can be found below.

### Cloning LAMMPS
The first step to building ML-MIX (either in-source or as a plugin) is cloning LAMMPS. LAMMPS can be cloned directly with
```
git clone -b release https://github.com/lammps/lammps.git lammps
```

### Building ML-MIX in source

ML-MIX can be built directly in the LAMMPS source, i.e by copying the `fix` and `pair` files into the `lammps/src` directory before starting the build. This *must* be done to use ml-mix-kokkos.

Run the `install.sh` script located in the `LAMMPS_plugin/` directory, passing the path to LAMMPS. Use the following command:

```
./install.sh /path/to/lammps/
```

#### Build *without* KOKKOS

An example cmake LAMMPS build script to compile and build LAMMPS without KOKKOS is provided in `build_scripts/example_lammps_build_script.sh`. This script includes all required packages.

Instructions to build LAMMPS can be found on the [official website](https://docs.lammps.org/Build.html). To use ML-MIX as a plugin you must build with the [PLUGIN package](https://docs.lammps.org/Packages_details.html#pkg-plugin) enabled. 

Additionally, to run the case-studies, the following packages are necessary:
- ML-PACE, which can be added using the instructions found [here](https://acesuit.github.io/ACEpotentials.jl/v0.6/tutorials/lammps/).
- ML-UF3
- RIGID (for the Si stretched bond case study)
- REPLICA (for nudged-elastic-band calculations)

Whilst it is not necessary for core functionality of the ML-MIX plugin, the case studies also depend on passing commands to the LAMMPS shared library through Python. To build LAMMPS in such a way that this is possible, please follow the instructions [here](https://docs.lammps.org/Python_install.html). Installing into a virtual environment is strongly encouraged.

#### Build *with* KOKKOS

> **Warning**
>
> Currently, ML-MIX-kokkos is in beta; when using it you are likely to encounter bugs, and only a limited number of KOKKOS pair_styles have been tested.

A kokkos-enabled build must be done *in-source* - as detailed above, the first step is therefore to run `/install.sh /path/to/lammps/` from the `LAMMPS_Plugin/` directory.

A minimal example `kokkos` build script is given in `build_scripts/example_lammps_build_script_with_kokkos.sh`. The build is more complicated, but is not made extra complicated by `ML-MIX`, and no extra compilation flags are required beyond those detailed in the [LAMMPS KOKKOS documentation](https://docs.lammps.org/Build_extras.html#kokkos). Note that the compilation *will be significantly slower* - this is due to a lot of Kokkos code being templated c++.

### Building ML-MIX as a plugin
> **Note**
>
> If ML-MIX was built in-source, this section can be skipped. You cannot do this if building with kokkos.

The ML-MIX plugin can be compiled straightforwardly with cmake in the `LAMMPS_plugin/` folder. The path to the `src` dir of LAMMPS must be provided, e.g.
```
mkdir build
cd build
cmake .. -D LAMMPS_SOURCE_DIR=/path/to/lammps/src \
cmake --build . -j 1
```
This code is also provided in `build_scripts/example_plugin_build_script.sh`.

### Loading the ML-MIX plugin
To use the ML-MIX plugin, both the `fix` and `pair_style` need to be loaded at the start of a LAMMPS input script:
```
plugin load path/to/LAMMPS_plugin/build/hybridoverlaymlmlplugin.so
plugin load path/to/LAMMPS_plugin/build/mlmlplugin.so
```

If the plugin has loaded correctly, the following lines should print
```
Loading plugin: MLML hybrid overlay pair style v0.1 by Fraser Birks (fraser.birks@warwick.ac.uk)
Loading plugin: MLML fix style v0.1 by Fraser Birks (fraser.birks@warwick.ac.uk)
```

### Installing Python packages

To install all the necessary python packages, please run (in the same Python environment which LAMMPS is installed) 
```
pip install -r requirements.txt
```

### Installing ACEpotentials.jl for constrained linear ACE fitting

We fit constrained linear ACE potentials using the [ACEpotentials.jl package](https://github.com/ACEsuit/ACEpotentials.jl). For this, Julia is needed, which can be installed following the instructions [here](https://docs.julialang.org/en/v1/manual/installation/). Install the environment contained in the `Manifest.toml` and `Project.toml` files with
```
julia --project=. -e 'using Pkg; Pkg.Registry.add("General"); Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry")); Pkg.instantiate()'
```

### Installing UF3 for constrained UF3 potential fitting

We fit constrained linear UF3 potentials using a slightly modified version of the [UF3 package](https://github.com/uf3/uf3). The necessary fork of this package is included as a submodule in this repository. We utilise the `UltraFastFeaturization` branch of UF3, which greatly speeds up generating the potential features for fitting. Building this requires the HDF5 library is installed, as detailed [here](https://github.com/uf3/uf3/tree/UltraFastFeaturization/UltraFastFeaturization). Before installation, set the environment variables

```
export HDF5_INCLUDE_DIR=/Path/to/HDF5/include
export HDF5_LIB_DIR=/Path/to/HDF5/lib
export ULTRA_FAST_FEATURIZER=True
```


and then proceed with the installation:

```
cd external/uf3
pip install .
```

## Test results for each feature
<!-- feature error table start -->
| Potential | Serial | Date and time | Parallel | Date and time |
| --- | --- | --- | --- | --- |
| Force zeroing | ✅ | 2025-03-13 15:54:08 | N/A | N/A |
| Region building | ✅ | 2025-05-13 15:38:53 | ✅ | 2025-05-13 15:23:57 |

<!-- feature error table end -->


## Tests results for different potentials
<!-- error table start -->
| Potential | Serial | Date and time | Parallel | Date and time |
| --- | --- | --- | --- | --- |
| LJ | ✅ | 2025-05-29 11:16:43 | ✅ | 2025-02-27 17:29:28 |
| table | ✅ | 2025-03-12 11:51:38 | ✅ | 2025-02-27 22:10:00 |
| ACE | ✅ | 2025-03-13 15:48:44 | ✅ | 2025-05-01 14:38:54 |
| SNAP | ✅ | 2025-06-03 17:51:59 | N/A | N/A |

<!-- error table end -->
