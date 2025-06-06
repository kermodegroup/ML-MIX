Setup
=====

Cloning this Repository
-----------------------
To clone this repository (and its submodule dependencies), run:

.. code-block:: bash

    git clone --recurse-submodules https://github.com/kermodegroup/ML-MIX.git

Cloning LAMMPS
--------------
Before building ML-MIX, you must clone the LAMMPS source code:

.. code-block:: bash

    git clone -b release https://github.com/lammps/lammps.git lammps

Building ML-MIX
===============

There are two ways to build ML-MIX:

1. In-source (required for KOKKOS integration)
2. As a plugin

Instructions for both are given below.

In-Source Build (Required for KOKKOS)
-------------------------------------
This method copies the `fix` and `pair` files into the `lammps/src/` directory. Run the installation script:

.. code-block:: bash

    ./install.sh /path/to/lammps/

**Without KOKKOS**

To build without KOKKOS, a reference build script is provided:

- `build_scripts/example_lammps_build_script.sh`

Required LAMMPS packages:

- PLUGIN
- ML-PACE (see `ACEpotentials.jl tutorial <https://acesuit.github.io/ACEpotentials.jl/v0.6/tutorials/lammps/>`_)
- ML-UF3
- RIGID (required for the Si stretched bond case study)
- REPLICA (required for nudged-elastic-band calculations)

The ML-MIX case studies require Python access to LAMMPS via its shared library interface. To enable this, follow the `LAMMPS Python installation guide <https://docs.lammps.org/Python_install.html>`_.
Installing into a virtual environment is strongly recommended.

**With KOKKOS**

.. warning::

    ML-MIX-kokkos is currently in beta. Bugs are expected and only a limited number of KOKKOS pair styles are supported.

An in-source build is required for KOKKOS. Use the `install.sh` script as above. Then follow the example script:

- `build_scripts/example_lammps_build_script_with_kokkos.sh`

Refer to the `LAMMPS KOKKOS build documentation <https://docs.lammps.org/Build_extras.html#kokkos>`_ for guidance. Note that KOKKOS builds can be significantly slower to compile due to C++ templating.

Plugin Build (Without KOKKOS)
-----------------------------
This method builds ML-MIX separately using the plugin infrastructure. Run:

.. code-block:: bash

    mkdir build
    cd build
    cmake .. -D LAMMPS_SOURCE_DIR=/path/to/lammps/src
    cmake --build . -j 1

See also: `build_scripts/example_plugin_build_script.sh`.

Loading the ML-MIX Plugin
-------------------------
To use ML-MIX via the plugin interface, load both the `fix` and `pair_style` libraries in your LAMMPS input script:

.. code-block:: bash

    plugin load path/to/LAMMPS_plugin/build/hybridoverlaymlmlplugin.so
    plugin load path/to/LAMMPS_plugin/build/mlmlplugin.so

On successful load, you should see:

.. code-block:: text

    Loading plugin: MLML hybrid overlay pair style v0.1 by Fraser Birks (fraser.birks@warwick.ac.uk)
    Loading plugin: MLML fix style v0.1 by Fraser Birks (fraser.birks@warwick.ac.uk)

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