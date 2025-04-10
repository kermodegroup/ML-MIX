On-the-fly Potential Fitting
============================

In ML-MIX, it is possible to perform on-the-fly linear ACE potential fitting. This process involves fitting a cheap potential to an expensive potential *during a simulation*, directly to the atoms which will become the cheap domain in the subsequent mixed simulation. This allows smaller, highly targeted fits to be performed.

Config YAML file
----------------

The on-the-fly fitting process involves calling python and julia scripts present in the `mix_on_the_fly` folder of the ML-MIX directory during a LAMMPS simulation. These scripts read in a user-defined configuration file in YAML format which should contain all the information necessary to process the LAMMPS dump data and fit the cheap potential.

The format of this configuration file as well as all the possible keywords are described below. An example configuration file for the on-the-fly fitting of an Fe ACE potential is provided in the `simple_examples/mix_on_the_fly` folder of the ML-MIX directory.

Steps to perform on-the-fly fitting
-----------------------------------

1. **Add the on-fly fitting keyword, timestep and relevant paths**. In the settings line of the pair hybrid/overlay/mlml command, add the following:

```
pair_style hybrid/overlay/mlml on-fly <fit_tstep> <path_to_ML_MIX> <path_to_cfg_file> ...
```
This will trigger an on-the-fly fit on timestep <fit_tstep>.

2. **Set up the mixed pair_style**. Set up a mixed simulation as it will be when the fit is complete. Label your expensive potential as MLML region 1 (this will be used for all atoms up until the fit) and label the on-the-fly cheap potential to be region 2. As we are fitting a linear ACE potential, ensure to add both the many-body and two-body terms to the pair_style. The name of the cheap potential must match the name specified in the config file.

```
pair_coeff * * pace 1 1 jace_espresso.yace Fe
pair_coeff 1 1 table 1 1 jace_espresso_pairpot.table Fe_Fe
pair_coeff * * pace 2 2 otf_pot.yace Fe
pair_coeff 1 1 table 2 2 otf_pot_pairpot.table Fe_Fe
```

3. **Set up fix mlml**. Aside from the normal keywords, we also need to set up fix mlml to output a mask that we can use to only perform the fit to the environments that will be modelled by the cheap potential. This is done by appending the output-mask keyword to the end of the fix mlml command, e.g:

```
fix mlml_fix all mlml 1 6.0 6.0 4.0 fix_classify av_ca 100 2.0 inf output-mask
```

This makes the fix output per-atom vector (accessible as `f_mlml_fix` in the dump command) an integer mask which is 1 for atoms that will be in the cheap region and 0 for atoms that will be in the expensive region.

4. **Optional: Store forces**. If you want to fit to atoms which have modified forces due to the action of a thermostat, you can define a fix to store the raw forces. This is done by adding the following line to the input script:

```
fix sf all store/force
```
Now you can dump the raw forces (just pair_style) to a file by accessing the `f_sf[1]`, `f_sf[2]` and `f_sf[3]` per-atom vectors in the dump command.

5. **Optional: Define per-atom-energies'** If you want to fit to per-atom-energies, you can define a compute to calculate these. This is done by adding the following lines to the input script:

```
compute pea all pe/atom
```

Now you can dump the per-atom energies to a file by accessing the `c_pea` per-atom vector in the dump command.

4. **Define a dump command**. This is the command to build the data that will be used for the fit. An example dump command (including stored forces, per-atom energies and the mlml mask) is shown below:

```
dump myDump all custom 100 fit_dump.lammpstrj id type mass xs ys zs f_sf[1] f_sf[2] f_sf[3] f_mlml_fix
```

You must include `mass`, `xs`, `ys` and `zs` and forces. Dumping per-atom-energies and the mlml mask is optional. Other fields will be ignored during the fitting. The dump should be plain text.

5. **Run the simulation**. Run the simulation as normal. A place-holder ACE potential will be generated initially at the start of the simulation, and the on-the-fly fitting will be triggered later at the timestep specified in the pair_style command. Note that it will only be triggered once in a simulation run, even if the timestep is reset. 


Automatic masking
-----------------

If you want to fit a potential which contains fewer elements than are present in the full expensive simulation, e.g fitting a pure W ACE potential to a mixed W-He simulation, you can specify to only fit to W in the `fit-elements` section of the configuration file, and then use the `auto-mask` keyword. This will automatically mask out any elements not included in the `fit-elements` list. The `auto-mask-cutoff` keyword specifies the cutoff distance for this masking. 

During the data processing step, any elements which are not present in `fit-elements` will be deleted. To avoid fitting errors due to fitting to expensive potential atoms with force components from the deleted species, the setting `auto-mask-cutoff` value should be set to twice the expensive potential cutoff. In practice, a value smaller than this can be used. 

###########################
Configuration File Reference
###########################

This configuration file is used to set up on-the-fly simulations in ML-MIX. Below is a detailed explanation of all available fields.

=========================
General Data Settings
=========================

**dump-name** (str)
    The name of the LAMMPS trajectory dump file where simulation data is stored.

**fit-elements** (list)
    A list of elements to include in the fitting process. Any element not listed will be automatically masked out and deleted before the fit if auto-mask is specified.

**auto-mask** (dict, optional)
    Automatic masking settings for elements not included in `fit-elements`:

    - **auto-mask-cutoff** (float)
        The cutoff distance for automatic masking (recommended to be twice the expensive potential cutoff).

**mask-key** (str, optional)
    The key used in the LAMMPS dump file for fix output vector data.

**pae-key** (str, optional)
    The key used in the LAMMPS dump file for per-atom energy data.

**force-keys** (dict)
    The header keys for the force components in the LAMMPS dump file:

    - **fx** (str): X-component of force.
    - **fy** (str): Y-component of force.
    - **fz** (str): Z-component of force.

=========================
ACEfit Settings
=========================

**e_refs** (dict)
    Reference energy values for each fitted element:

    - **Element Symbol** (float): Reference energy (e.g., `Fe: -3456.00113622806`).

**order** (int)
    The correlation order of the ACE basis.

**totaldegree** (int)
    The maximum total degree of the ACE basis.

**rcut** (float)
    Cutoff radius for neighbor interactions in the ACE fit.

**smoothness_prior_strength** (float)
    The strength of the smoothness prior, which controls regularization.

=========================
Output Settings
=========================

**output_name** (str)
    The name of the output file where the trained ACE potential will be saved.

=========================
Example Configuration File
=========================

Below is an example YAML configuration file for the on-the-fly fitting of an Fe ACE potential:

.. code-block:: yaml

    dump-name: fit_dump.lammpstrj
    fit-elements:
      - Fe
    mask-key: f_mlml_fix
    force-keys:
      fx: f_sf[1]
      fy: f_sf[2]
      fz: f_sf[3]
    e_refs:
      Fe: -3456.00113622806
    order: 2
    totaldegree: 10
    rcut: 6.0
    smoothness_prior_strength: 4
    output_name: otf_pot

An on-the-fly fitting example that uses this configuration file can be found in the `simple_examples/mix_on_the_fly` folder of the ML-MIX directory.
