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

    pair_style hybrid/overlay/mlml lj/cut 5.0 pace
    pair_coeff * * lj/cut 2 1.0 1.0
    pair_coeff * * pace 1 Si.yace Si

    pair_style hybrid/overlay/mlml lj/cut 5.0 uf3 3
    pair_coeff * * lj/cut 2 1.0 1.0
    pair_coeff * * uf3 1 W.uf3 W

Description
"""""""""""

.. warning::

   **Bug Warning**  
   In the current version of LAMMPS, if a `pair_style` attempts to build a half neighbor list by pruning a full neighbor list created for a `fix`, it can result in all forces computed by that `pair_style` being zero (atoms effectively have no neighbors).

   This issue affects `ML-MIX` in the following cases:

   - `fix mlml` is defined (which requires a full neighbor list).
   - *All* other `pair_style` definitions require half neighbor lists.  
     (If at least one `pair_style`—such as `ACE` or `UF3`—requires a full neighbor list, then half neighbor lists are constructed correctly, and this issue does not occur.)

   We are actively working to resolve this issue.

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
