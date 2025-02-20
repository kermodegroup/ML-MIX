# ML-MIX ![CI](https://github.com/kermodegroup/ML-MIX/actions/workflows/ci.yml/badge.svg)
A LAMMPS plugin for efficiently mixing together forces from different machine-learned interatomic potentials inspired by quantum mechanics/molecular mechanics (QM/MM) methods. Use expensive, complex potentials only where they're needed! 


<div style="display: flex; gap: 10px;">
  <img src="docs/images/blue_background_title_looped.gif" width="500">
  <img src="docs/images/Fe_dumbbell.gif" width="300">
</div>

## Key Information

Full documentation can be found here: [documentation link]

Detailed setup instructions can be found below.

For a user wanting to start quickly, simple examples showcasing and explaining features of ML-MIX can be found in `simple_examples/`

The accompanying paper on ML-MIX can be found here: [paper link]

All the code associated with this paper is included in this repository, in `potential_fitting/`, `case_studies/` and `ancillary_studies/`, however the GitHub version is incomplete, as it is missing potentials, data and plots. The full version of this repository can be found on Zenodo [zenodo link].

Scripts and examples for constrained fitting of cheap potentials can be found in `constrained_fitting/`.

Code for running all of the main studies in the paper can be found in `case_studies/`. 

Code for running all the ancillary studies in the paper can be found in `ancillary_studies/`.

More information can be found in the `README` files inside each subdirectory.

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
Instructions to build LAMMPS can be found on the [official website](https://docs.lammps.org/Build.html). To use ML-MIX, you must build with the [PLUGIN package](https://docs.lammps.org/Packages_details.html#pkg-plugin) enabled. 

Additionally, to run the case-studies, the following packages are necessary:
- ML-PACE, which can be added using the instructions found [here](https://acesuit.github.io/ACEpotentials.jl/v0.6/tutorials/lammps/).
- ML-UF3
- RIGID (for the Si stretched bond case study)
- REPLICA (for nudged-elastic-band calculations)

Whilst it is not necessary for core functionality of the ML-MIX plugin, the case studies also depend on passing commands to the LAMMPS shared library through Python. To build LAMMPS in such a way that this is possible, please follow the instructions [here](https://docs.lammps.org/Python_install.html). Installing into a virtual environment is strongly encouraged.

To aid the user, an example cmake LAMMPS build script is provided in `build_scripts/example_lammps_build_script.sh`.

### Building the ML-MIX plugin
The ML-MIX plugin can be compiled straightforwardly with cmake in the LAMMPS_plugin folder. The path to the `src` dir of LAMMPS must be provided, e.g.
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
| Region building | ✅ | 2025-02-18 16:43:39 | ✅ | 2025-02-17 18:57:25 |

<!-- feature error table end -->


## Tests results for different potentials
<!-- error table start -->
| Potential | Serial | Date and time | Parallel | Date and time |
| --- | --- | --- | --- | --- |
| EAM/fs | ✅ | 2025-02-12 16:31:51 | ✅ | 2025-02-12 16:33:36 |
| EAM/cd | ✅ | 2025-02-12 16:31:49 | ✅ | 2025-02-12 16:33:35 |
| LJ | ✅ | 2025-02-12 16:32:23 | ✅ | 2025-02-12 16:33:58 |
| UF3 | ✅ | 2025-02-12 16:31:47 | ✅ | 2025-02-12 16:33:34 |
| ACE | ✅ | 2025-02-12 16:32:22 | ✅ | 2025-02-12 16:33:57 |
| EAM/alloy | ✅ | 2025-02-12 16:31:48 | ✅ | 2025-02-12 16:33:35 |
| table | ✅ | 2025-02-12 16:31:48 | ✅ | 2025-02-12 16:33:34 |
| EAM | ✅ | 2025-02-12 16:31:50 | ✅ | 2025-02-12 16:33:35 |
| EAM/he | ✅ | 2025-02-12 16:31:51 | ✅ | 2025-02-12 16:33:36 |

<!-- error table end -->
