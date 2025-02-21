.. index:: fix langevin_mlml

fix langevin/mlml command
====================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID langevin/mlml Tstart Tstop damp seed mlmix_region keyword values ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* langevin = style name of this fix command
* Tstart,Tstop = desired temperature at start/end of run (temperature units)
* Tstart can be a variable (see below)
* damp = damping parameter (time units)
* seed = random number seed to use for white noise (positive integer)
* mlmix_region = region ID of the mlml region to apply the fix to
* zero or more keyword/value pairs may be appended
* keyword = *angmom* or *gjf* or *omega* or *scale* or *tally* or *zero*

  .. parsed-literal::

       *angmom* value = *no* or factor
         *no* = do not thermostat rotational degrees of freedom via the angular momentum
         factor = do thermostat rotational degrees of freedom via the angular momentum and apply numeric scale factor as discussed below
       *gjf* value = *no* or *vfull* or *vhalf*
         *no* = use standard formulation
         *vfull* = use Gronbech-Jensen/Farago formulation
         *vhalf* = use 2GJ formulation
       *omega* value = *no* or *yes*
         *no* = do not thermostat rotational degrees of freedom via the angular velocity
         *yes* = do thermostat rotational degrees of freedom via the angular velocity
       *scale* values = type ratio
         type = atom type (1-N)
         ratio = factor by which to scale the damping coefficient
       *tally* value = *no* or *yes*
         *no* = do not tally the energy added/subtracted to atoms
         *yes* = do tally the energy added/subtracted to atoms
       *zero* value = *no* or *yes*
         *no* = do not set total random force to zero
         *yes* = set total random force to zero

Examples
""""""""

.. code-block:: LAMMPS

  fix mlml_langevin all langevin/mlml 900 900 2.0 12345 2

Description
"""""""""""

Apply a Langevin thermostat, as described in the main LAMMPS documentation, with the difference that this fix applies the thermostat to a specific region of the simulation domain designated by the *mlmix* package. Which region the thermostat is applied to is determined by the *mlmix_region* argument. The *mlmix_region* argument must either be 1 or 2.

Restrictions
""""""""""""
To use this fix, the i2_potential and d2_eval property/atoms must be defined.
The group used with this command must be all. 

Related commands
""""""""""""""""

:doc:`pair_style hybrid/overlay/mlml <pair_hybrid_overlay_mlml>`
:doc:`fix mlml <fix_mlml>`

Default
"""""""

The option defaults are angmom = no, omega = no, scale = 1.0 for all
types, tally = no, zero = no, gjf = no.