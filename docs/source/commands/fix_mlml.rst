.. _fix-mlml:

.. index:: fix mlml

fix mlml
========

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID mlml nevery r_core r_buff r_blend keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* mlml = style name of this fix command
* nevery = rebuild the mlml regions every this many timesteps
* r_core = size of the core region (distance units)
* r_buff = size of the buffer region (distance units)
* r_blend = size of the blending region (distance units)
* keyword = *group* or *fix_classify* or *init_group* or *hysteresis-time*
* one of *group* or *fix_classify* keyword/args must be appended
* *init_group* keyword/args are optional and can only be appended after *fix_classify* args
* *hysteresis-time* keyword/args are optional and can only be appended after *group* or *fix_classify* args

  .. parsed-literal::
       *group* arg = *group-ID*
         group-ID = ID of a group of seed atoms to build region around
       *fix_classify* args = *fix_id* *nfreq* *lb* *ub*
         *fix_id* = ID of a fix that generates output vector for classification
         *nfreq* = querying classifying fix every this many timesteps
         *lb* = lower bound of fix output vector values to select as seed atoms
         *ub* = lower bound of fix output vector values to select as seed atoms
       *init_group* args = *group-ID*
         group-ID = ID of a group of seed atoms to initially build mlml regions around before the first time the fix output vector defined by *fix_classify* is queried
       *hysteresis-time* args = *in_time* *out_time*
         *in_time* = time decay constant (time units) for atoms entering the expensive potential region (region 1).
         *out_time* = time decay constant (time units) for atoms leaving the expensive potential region (region 1).
       *blend* args = *blend_type*
         *blend_type* = type of blending to use for the core region
           *linear* = linear blending (default)
           *cubic* = cubic blending
         
Examples
""""""""

.. code-block:: LAMMPS

    group He_group type 2
    fix mlml_fix all mlml 1 4.0 6.0 4.0 group He_group

    compute ca all coord/atom cutoff 2.0
    fix av_ca all ave/atom 10 10 100 c_ca
    fix mlml_fix all mlml 1 4.0 6.0 4.0 fix_classify av_ca 100 0.4 inf

    compute ca all coord/atom cutoff 2.0
    fix av_ca all ave/atom 10 10 100 c_ca
    group init_defect id 10 11 12
    fix mlml_fix all mlml 1 4.0 6.0 4.0 fix_classify av_ca 100 0.4 inf init_group init_defect

    group He_group type 2
    fix mlml_fix all mlml 1 4.0 6.0 4.0 group He_group hysteresis-time 0.001 0.01



Description
"""""""""""
.. warning::

   **Bug Warning**
    Update: This issue has been resolved in the `17/03/2025 LAMMPS Stable release <https://github.com/lammps/lammps/releases/tag/stable_29Aug2024_update2>`_. Please ensure you are using this version or later to avoid the issue described below.

    In the current version of LAMMPS, if a `pair_style` attempts to build a half neighbor list by pruning a full neighbor list created for a `fix`, it can result in all forces computed by that `pair_style` being zero (atoms effectively have no neighbors).

    This issue affects `ML-MIX` in the following cases:

    - `fix mlml` is defined (which requires a full neighbor list).
    - *All* other `pair_style` definitions require half neighbor lists.  
      (If at least one `pair_style`—such as `ACE` or `UF3`—requires a full neighbor list, then half neighbor lists are constructed correctly, and this issue does not occur.)

This fix command is used to create and maintain a set of regions to be evaluated with 2 different pair_styles using the *hybrid/overlay/mlml* pair style in the *mlmix* package. In *hybrid/overlay/mlml*, different pair_styles are labelled either 1 or 2 to indicate evaluation region.

This fix is used to build an evaluation region around user-specified seed atoms. All atoms which lie inside this region are marked for evaluation with pair_style 1, and all atoms outside it are marked for evaluation with potential 2.

The *group* keyword can be used to specify a fixed LAMMPS group of seed atoms to build regions around. The *fix_classify* keyword can be used to specify a different fix with ID *fix_id* that generates an output per-atom vector, which is used to classify atoms as seed atoms. Atoms which correspond to values between *lb* and *ub* of the fix output vector are used as seed atoms. *lb* can take the value of *-inf* and *ub* can take the value of *inf* to specify no lower or upper bound. *nfreq* specifies how often the fix output vector is queried. The evaluation frequency of the target fix must evenly divide *nfreq*.

.. warning::
    If using the *fix_classify* keyword, this fix must be defined
    after the fix that is queried for the output vector in the input script.

The *init_group* keyword can be optionally used after *fix_classify* to specify a group of seed atoms to initially build the regions around. If this is not specified, then all atoms are flagged to be evaluated with pair_style 1 until the first time the output vector of the fix specified with *fix_classify* is queried.

The regions are rebuilt every *nevery* timesteps. If *nevery* is set to 0, the regions are only built once at the start of the simulation.

Regions are built around all seed atoms with a core region, a blending region
and a buffer region.

The core region is the union of atoms contained within spheres of radius *r_core* around each seed atom.

The blending region is the union of atoms contained within spheres of radius 
*r_blend* around each core atom, which is not already in the core region.
Atoms within the blending region are evaluated with both pair_styles, and the forces are linearly blended between the two. The proportion of pair_style 1 force on any individual atom is determined by

.. math::

   p_{1} = 1.0 - \left(\frac{|\mathbf{r}|}{r_{\text{blend}}}\right)

for *blend linear* and

.. math::

   p_{1} = 1.0 - \left(3\left(\frac{|\mathbf{r}|}{r_{\text{blend}}}\right)^{2} - 2\left(\frac{|\mathbf{r}|}{r_{\text{blend}}}\right)^{3}\right)

**New in 0.2.0** for *blend cubic*, where :math:`\mathbf{r}` is the shortest vector from the atom to any seed atom.

The force on this blended atom is then determined by

.. math::

   \mathbf{F}^{i} = p_{1} \mathbf{F}^{i}_{1} + (1 - p_{1}) \mathbf{F}^{i}_{2}


There are two buffer regions, which are each constructed by taking the union of atoms contained within spheres of radius *r_buff* around blending atoms. The pair_style 1 buffer are atoms external to the blending and core regions, whilst the pair_style 2 buffer is only atoms contained within the core region. Note that if `*r_buff* > *r_core*`, pair_style 2 buffer will contain all core atoms. 

**New in 0.2.0** If *time-decay-hysteresis* is enabled, the proportion of pair_style 1 force on each atom is updated at each region rebuild step using a discretised exponential decay. Rather than immediately transitioning to the :math:`p_{1}` values computed above, an exponentially decaying ramp is used from the previous :math:`p_{1}^{\text{prev}}` values to the newly computed set of :math:`p_{1}` values (the target :math:`p_{1}^{\text{target}}`). This introduces a smooth ramping behaviour.

Two characteristic times control this:

- *in_time* — controls the rate at which atoms ramp *into* the region.
- *out_time* — controls the rate at which atoms ramp *out of* the region.

The new proportion of potential 1 force on each atom is computed each update as:

.. math::

   p_{1} \leftarrow \Delta \cdot \min\left(\frac{dt \cdot nevery}{\tau}, 1\right) + p_{1}^{\text{prev}}

where:

- :math:`\Delta = p_{1}^{\text{target}} - p_{1}^{\text{prev}}`
- :math:`\tau` is either *in_time* or *out_time* depending on the sign of :math:`\Delta`
- *dt* is the simulation timestep
- *nevery* is the number of steps between region rebuilds

Thresholding is then applied such that atoms with :math:`p_{1} < 0.01` are evaluated with only potential 2 and atoms with :math:`p_{1} > 0.99` are evaluated with only potential 1. 

This approach ensures a smooth transition in force evaluation when atoms move across region boundaries, reducing the 'flickering' behaviour that can be seen when atoms move quickly between regions.

Acceleration
""""""""""""
.. warning::

  Kokkos acceleration for this fix is experimental and still in development. Some features may not work as expected, please report any issues you encounter.

This fix has a kokkos enabled version, which exists in order to allow ML/ML simulations to be performed using GPU accelerated kokkos pair_styles. Please see the restrictions section below for information on how the use of this differs from the CPU version.

Restrictions
""""""""""""
**CPU Version:**
To use this fix, the i2_potential and d2_eval property/atoms must be defined.
This fix is designed to be used in conjunction with the *hybrid/overlay/mlml* pair style.

**GPU Version:**
To use this fix, the d_potential_1, d_potential_2, d_eval_1 and d_eval_2 property/atoms must be defined (for device compatibility).


This fix is designed to be used in conjunction with the *hybrid/overlay/mlml* pair style.

Please see :doc:`pair_style hybrid/overlay/mlml <pair>` for more information.

The group specified with this command must be all.

Related commands
""""""""""""""""

:doc:`pair_style hybrid/overlay/mlml <pair_hybrid_overlay_mlml>`

Default
"""""""

none
