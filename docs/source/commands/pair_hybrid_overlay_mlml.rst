.. index:: pair_style hybrid/overlay/mlml

pair_style hybrid/overlay/mlml
==============================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay/mlml keyword value style1 args style2 args ...

* keyword = zero
* style1,style2 = list of one or more pair styles and their arguments
* factor1,factor2 = scale factors for pair styles, may be a variable

* zero or more keyword/value pairs may be appended
* keyword = *angmom* or *gjf* or *omega* or *scale* or *tally* or *zero*

  .. parsed-literal::

       *zero* value = *no* or *yes*
         *no* = do not set total pair_style force to zero
         *yes* = set total pair_style force to zero


Examples
""""""""

.. code-block:: LAMMPS

    pair_style hybrid/overlay/mlml lj/cut 2.5 lj/cut 3.5
    pair_coeff * * lj/cut 1 1 1.0 1.0
    pair_coeff * * lj/cut 2 2 1.5 0.8

    pair_style hybrid/overlay/mlml lj/cut 2.5 eam/alloy
    pair_coeff * * lj/cut 1 1.0 1.0
    pair_coeff * * eam/alloy 2 AlCu.eam.alloy Al Cu

Description
"""""""""""

This documentation should be considered as an extension to the documentation provided for the *hybrid/overlay* pair style. Everything that is valid for the *hybrid/overlay* pair style is also valid for the *hybrid/overlay/mlml* unless stated otherwise here.

The *hybrid/overlay/mlml* style is a variant of the *hybrid/overlay* style which allows for the use of multiple pair styles in one simulation. The *hybrid/overlay/mlml* style is part of the *mlmix* package. With the *hybrid/overlay/mlml* style, different pair styles can be evaluated in different spatial regions of the simulation domain. The full force vector is constructed through force-mixing, a method commonly used in QM/MM simulations.

The assignment of pair styles to type pairs is made via the :doc:`pair_coeff <pair_coeff>` command.  For *hybrid/overlay/mlml*, an additional numeric argument must be specified after the sub-style name, which can either be 1 or 2. This value indicates the *mlmix* region to which the sub-style should be  applied.

For *hybrid/overlay/mlml*, if a pair style is used multiple times in the pair_style command, the additional numeric argument to differentiate between the different instances of the sub-style must be specified after the number which indicates the *mlmix* region.

----------

Forces from force mixing are not the derivatives of a potential energy function, subsequently the 'total energy' of the system is not well defined. Therefore, the *hybrid/overlay/mlml* will set global potential energy and virials to zero.

.. However, per-atom-energies and per-atom-virials are still calculated.

A consequence of force mixing is that the total force on the system may not be zero. To address this, the *zero* keyword can be used. If set to *yes*, the total force on the system will be modified such that the total is zero. This is achieved by averaging the net overall force on the system over all atoms and subtracting this from the force on each atom.

Restrictions
""""""""""""

As well as the restrictions that apply to the *hybrid/overlay* pair style, :doc:`run_style respa <run_style>` is completely disabled for *hybrid/overlay/mlml*. 

To use this pair_style, the i2_potential and d2_eval property/atoms must be defined.

This pair_style is designed to be used in conjunction with *fix mlml*. Please see :doc:`fix mlml <fix>` for more information.


Related commands
""""""""""""""""

:doc:`pair_coeff <pair_coeff>`
:doc:`fix mlml <fix_mlml>`

Default
"""""""

none
