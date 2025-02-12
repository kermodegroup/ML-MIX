# ML-MIX
Package for ML/ML Force Mixing with Machine Learning Potentials in LAMMPS

## Setup

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

Whilst it is not necessary for the ML-MIX plugin, the examples also depend on passing commands to the LAMMPS shared library through Python. To build LAMMPS in such a way that this is possible, please follow the instructions [here](https://docs.lammps.org/Python_install.html). Installing into a virtual environment is strongly encouraged.

To aid the user, an example cmake LAMMPS build script is provided in `examples/example_lammps_build_script.sh`.

### Building the ML-MIX plugin
The ML-MIX plugin can be compiled straightforwardly with cmake in the LAMMPS_plugin folder. The path to the `src` dir of LAMMPS must be provided, e.g.
```
-D LAMMPS_SOURCE_DIR=/path/to/lammps/src
```
For an example cmake script, please see `examples/example_plugin_build_script.sh`.

### Loading the ML-MIX plugin
To use the ML-MIX plugin, both the `fix` and `pair_style` need to be loaded at the start of an input script, e.g. with
```
plugin load path/to/LAMMPS_plugin/build/hybridoverlaymlmlplugin.so
plugin load path/to/LAMMPS_plugin/build/mlmlplugin.so
```

If the plugin has loaded correctly, the following lines should print
```
Loading plugin: MLML hybrid overlay pair style v0.1 by Fraser Birks (fraser.birks@warwick.ac.uk)```
Loading plugin: MLML fix style v0.1 by Fraser Birks (fraser.birks@warwick.ac.uk)
```

### Installing Python packages

To install all the necessary python packages, please run (in the same Python environment which LAMMPS is installed) 
```
pip install -r requirements.txt
```

### Installing ACEpotentials.jl for constrained linear ACE fitting

We fit constrained linear ACE potentials using the [ACEpotentials.jl package](https://github.com/ACEsuit/ACEpotentials.jl). For this, Julia is needed, which can be installed following the instructions [here](https://docs.julialang.org/en/v1/manual/installation/). Once this has been installed, open the Julia REPL with 
```
julia
```
and install the environment according to the `Manifest.toml` and `Project.toml` files with(".")```

```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

### Installing UF3 for constrained UF3 potential fitting

We fit constrained linear UF3 potentials using a slightly modified version of the [UF3 package](https://github.com/uf3/uf3). The necessary fork of this package is included as a submodule in this repository. Specifically, we are utilising the `UltraFastFeaturization` branch of UF3, which greatly speeds up generating the potential features for fitting. Building this requires the HDF5 library is installed, as detailed [here](https://github.com/uf3/uf3/tree/UltraFastFeaturization/UltraFastFeaturization). Set the environment variables

```
export HDF5_INCLUDE_DIR=/Path/to/HDF5/include
export HDF5_LIB_DIR=/Path/to/HDF5/lib
export ULTRA_FAST_FEATURIZER=True
```


Once these are set, navigate to the uf3 submodule and install it:

```
cd external/uf3
pip install .
```

