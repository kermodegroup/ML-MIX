Setup
=====

Cloning this Repository
-----------------------
To clone this repository (and its submodule dependencies), run:

.. code-block:: bash

    git clone --recurse-submodules https://github.com/kermodegroup/ML-MIX.git

Installing LAMMPS
-----------------
To compile the ML-MIX LAMMPS plugin, first clone and build LAMMPS:

.. code-block:: bash

    git clone -b release https://github.com/lammps/lammps.git mylammps

Follow the `official build instructions <https://docs.lammps.org/Build.html>`_.

Required LAMMPS packages:

- ML-PACE
- ML-UF3
- RIGID (for Si stretched bond case study)
- REPLICA (for nudged-elastic-band calculations)

The ML-MIX case studies (given in `case_studies/`) use Python-LAMMPS integration. 
To install LAMMPS with this enabled, follow the `Python installation guide <https://docs.lammps.org/Python_install.html>`_.
Installation into a virtual environment is highly recommended.

**Example LAMMPS build script:** `build_scripts/example_lammps_build_script.sh`.

Building the ML-MIX Plugin
--------------------------
To build the ML-MIX plugin, navigate to the `LAMMPS_plugin` directory and run:

.. code-block:: bash

    mkdir build
    cd build
    cmake .. -D LAMMPS_SOURCE_DIR=/path/to/lammps/src
    cmake --build . -j 1

Ensuring to replace `/path/to/lammps/src` with the path to the cloned LAMMPS source directory.

Example script: `build_scripts/example_plugin_build_script.sh`.

Loading the ML-MIX Plugin
-------------------------
To load the ML-MIX plugin files, add the following lines to your LAMMPS input script:

.. code-block:: bash

    plugin load path/to/LAMMPS_plugin/build/hybridoverlaymlmlplugin.so
    plugin load path/to/LAMMPS_plugin/build/langevinmlmlplugin.so
    plugin load path/to/LAMMPS_plugin/build/mlmlplugin.so

If the plugins are loaded successfully, you should see outputs similar to the following for each:

.. code-block:: text

    Loading plugin: MLML hybrid overlay pair style v0.1 by Fraser Birks (fraser.birks@warwick.ac.uk)

Installing Python Packages
--------------------------
To install the Python packages required for running provided case studies, run the following:

.. code-block:: bash

    pip install -r requirements.txt

Ensure that the python enivronment is the one containing LAMMPS.

Installing ACEpotentials.jl for ACE fitting
-------------------------------------------
To install `ACEpotentials.jl <https://github.com/ACEsuit/ACEpotentials.jl>`_, run:

.. code-block:: bash

    julia --project=. -e 'using Pkg; Pkg.Registry.add("General"); Pkg.Registry.add(RegistrySpec(url="https://github.com/ACEsuit/ACEregistry")); Pkg.instantiate()'

Installing UF3 for UF3 fitting
------------------------------
We fit constrained linear UF3 potentials using a slightly modified version of the `UF3 package <https://github.com/uf3/uf3>`_. The necessary fork of this package is included as a submodule in this repository. We utilise the UltraFastFeaturization branch of UF3, which greatly speeds up generating the potential features for fitting. Building this requires the HDF5 library is installed, as detailed here. Before installation, the environment variables HDF5_INCLUDE_DIR and HDF5_LIB_DIR must be set to the HDF5 include and lib directories respectively.

.. code-block:: bash

    export HDF5_INCLUDE_DIR=/Path/to/HDF5/include
    export HDF5_LIB_DIR=/Path/to/HDF5/lib
    export ULTRA_FAST_FEATURIZER=True
    cd external/uf3
    pip install .