.. index:: pair_style eam/mlml
.. index:: pair_style eam/fs/mlml

pair_style eam/mlml command
=============================

For constant precision optimised variant: *eam*

pair_style eam/fs/mlml command
================================

For constant precision optimised variant: *eam/fs*

Syntax
""""""

.. code-block:: LAMMPS

   pair_style eam/mlml
   pair_style eam/fs/mlml

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay/mlml pace eam/fs/mlml
   pair_coeff * * pace 1 W.yace W
   pair_coeff * * eam/fs/full 2 W.eam.fs W


Description
"""""""""""
Style *eam* computes pairwise interactions for metals and metal alloys
using embedded-atom method (EAM) potentials :ref:`(Daw) <Daw3>`.  The total
energy :math:`E_i` of an atom :math:`i` is given by

.. math::

   E_i^\text{EAM} = F_\alpha \left(\sum_{j \neq i}\ \rho_\beta (r_{ij})\right) +
         \frac{1}{2} \sum_{j \neq i} \phi_{\alpha\beta} (r_{ij})

where :math:`F` is the embedding energy which is a function of the atomic
electron density :math:`\rho`, :math:`\phi` is a pair potential interaction,
and :math:`\alpha` and :math:`\beta` are the element types of atoms
:math:`i` and :math:`j`.  The multi-body nature of the EAM potential is a
result of the embedding energy term. Both summations in the formula are over
all neighbors :math:`j` of atom :math:`i` within the cutoff distance.
EAM is documented in detail in :doc:`pair_style eam <pair_eam>`.

This pair style is optimized for the usage with ML-MIX: Communication within
the force computation routine is prevented at the cost of some double
computations. Thereby, one can efficiently balance the computational load
of MLML-EAM simulations via :doc:`fix balance <fix_balance>`.

.. warning::

   The double computations computed by this pair style are disadvantageous
   when force-mixing via :doc:`pair hybrid/overlay/mlml <pair_hybrid_overlay_mlml>`
   is not used.

Force-mixing provides the force :math:`\pmb{F}_i` on an atom :math:`i` that
is interpolated between two models.

.. math::

   \pmb{F}_i^\text{MIX} = \lambda_i \pmb{F}_i^\text{(fast)} + (1-\lambda_i)
   \pmb{F}_i^\text{(precise)}\,,

where the switching parameter :math:`\lambda_i` is computed
dynamically during a simulation by :doc:`fix mlml <fix_mlml>`.

The pair style *eam/fs/mlml* computes the force
:math:`\pmb{F}_i^\text{EAM}` and should be combined
with a precise method like
:doc:`pair_style pace <pair_pace>` via
:doc:`pair_style hybrid/overlay/mlml <pair_hybrid_overlay_mlml>`.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

For atom type pairs I,J and I != J, where types I and J correspond to
two different element types, mixing is performed by LAMMPS as
described above with the individual styles.  You never need to specify
a pair_coeff command with I != J arguments for the eam/mlml styles.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

The eam/mlml pair styles do not write their information to :doc:`binary
restart files <restart>`, since it is stored in tabulated potential
files.  Thus, you need to re-specify the pair_style and pair_coeff
commands in an input script that reads a restart file.

The eam/mlml pair styles can only be used via the *pair* keyword of the
:doc:`run_style respa <run_style>` command.  They do not support the
*inner*, *middle*, *outer* keywords.

----------

Restrictions
""""""""""""

This pair styles are part of the MLML package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""

:doc:`pair_style eam  <pair_eam>`,
:doc:`pair_style hybrid/overlay/mlml <pair_hybrid_overlay_mlml>`,
:doc:`fix mlml <fix_mlml>`

Default
"""""""

none

----------

.. _Daw3:

**(Daw)** Daw, Baskes, Phys Rev Lett, 50, 1285 (1983).
Daw, Baskes, Phys Rev B, 29, 6443 (1984).
