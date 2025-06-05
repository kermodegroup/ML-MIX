# ML-MIX ![CI](https://github.com/kermodegroup/ML-MIX/actions/workflows/ci.yml/badge.svg) ![DOCS](https://github.com/kermodegroup/ML-MIX/actions/workflows/docs.yml/badge.svg) 
A LAMMPS plugin for efficiently mixing together forces from different machine-learned interatomic potentials inspired by quantum mechanics/molecular mechanics (QM/MM) methods. Use expensive, complex potentials only where they're needed! 


<div style="display: flex; gap: 10px;">
  <img src="docs/images/blue_background_title_looped.gif" width="500">
  <img src="docs/images/Fe_dumbbell.gif" width="300">
</div>

## Key Information

> **ML-MIX is under active development**  
> Major updates and collections of bugfixes will be periodically merged from the `develop` branch into `main` and released as new versions.
> If you encounter any bugs or unexpected behaviour, please [open an issue](https://github.com/yourusername/ML-MIX/issues) or contact us directly via email at [fraser.birks@warwick.ac.uk](fraser.birks@warwick.ac.uk).

Full documentation can be found here: https://kermodegroup.github.io/ML-MIX/

Detailed setup instructions can be found below.

For a user wanting to start quickly, simple examples showcasing and explaining features of ML-MIX can be found in `simple_examples/`

The accompanying paper on ML-MIX can be found here: [arXiv:2502.19081](https://arxiv.org/abs/2502.19081)

All the code associated with this paper is included in this repository, in `potential_fitting/`, `case_studies/` and `ancillary_studies/`, however the GitHub version is incomplete, as it is missing potentials, data and plots. The full version of this repository can be found on Zenodo [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14920348.svg)](https://doi.org/10.5281/zenodo.14920348)

More information can be found in the `README` files inside each subdirectory.

## ⚠️ Bug Warning ⚠️
In the current version of LAMMPS, if a pair_style tries to build a half neighborlist by pruning a full neighborlist that is built for a fix, it seems to lead to all the forces computed by that pair_style being 0 (atoms have no neighbors). This will be a problem when using ML-MIX if:
- fix mlml is defined (which needs a full neighborlist)
- *ALL* other pair_styles defined need half neighborlists (if there is even one pair_style defined which needs a full neighborlist, i.e, ACE, UF3, this isn't a problem as then half neighborlists are constructed correctly.)

Update: This issue has been resolved in the [17/03/2025 LAMMPS Stable release](https://github.com/lammps/lammps/releases/tag/stable_29Aug2024_update2).


## Setup

### Cloning this repository

To clone this repository (and it's submodule dependencies), run
```
git clone --recurse-submodules https://github.com/kermodegroup/ML-MIX.git
```

### Installing LAMMPS
To compile the ML-MIX LAMMPS plugin, you first need to clone and build LAMMPS. LAMMPS can be cloned directly with
```
git clone -b release https://github.com/lammps/lammps.git mylammps
```

An example cmake LAMMPS build script is provided in `build_scripts/example_lammps_build_script.sh` which includes all required plugins.

Instructions to build LAMMPS can be found on the [official website](https://docs.lammps.org/Build.html). To use ML-MIX, you must build with the [PLUGIN package](https://docs.lammps.org/Packages_details.html#pkg-plugin) enabled. 

Additionally, to run the case-studies, the following packages are necessary:
- ML-PACE, which can be added using the instructions found [here](https://acesuit.github.io/ACEpotentials.jl/v0.6/tutorials/lammps/).
- ML-UF3
- RIGID (for the Si stretched bond case study)
- REPLICA (for nudged-elastic-band calculations)

Whilst it is not necessary for core functionality of the ML-MIX plugin, the case studies also depend on passing commands to the LAMMPS shared library through Python. To build LAMMPS in such a way that this is possible, please follow the instructions [here](https://docs.lammps.org/Python_install.html). Installing into a virtual environment is strongly encouraged.

### Building the ML-MIX plugin
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
| Region building | ✅ | 2025-04-15 14:15:56 | N/A | N/A |
| Force zeroing | ✅ | 2025-03-13 15:54:08 | N/A | N/A |

<!-- feature error table end -->


## Tests results for different potentials
<!-- error table start -->
| Potential | Serial | Date and time | Parallel | Date and time |
| --- | --- | --- | --- | --- |
| table | ✅ | 2025-03-12 11:51:38 | ✅ | 2025-02-27 22:10:00 |
| UF3 | ✅ | 2025-04-15 14:40:22 | ✅ | 2025-03-12 12:12:02 |
| ACE | ✅ | 2025-03-13 15:48:44 | ✅ | 2025-03-12 12:12:01 |

<!-- error table end -->
