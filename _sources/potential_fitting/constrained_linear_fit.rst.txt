Constrained Potential Fitting
=============================

The constrained potential fitting process aims to fit a cheap potential to a highly localized region of the potential energy surface of an expensive reference model. This approach enables fitting a cheap potential which approximates an expensive potential well in a simulation while maintaining low computational complexity and cost. 

Selecting Constraints
---------------------

The fitting process involves two types of constraints: **soft constraints** and **hard constraints**.

- **Soft Constraints**: These are loose matching conditions that allow flexibility in the fit, suitable for general properties of the potential.
  
- **Hard Constraints**: These are tight matching conditions, often used for critical properties like elastic constants. For example, ensuring the seamless matching of long-range elastic stress fields between the cheap and expensive potentials in solid-state simulations.

The constraints are encoded as design matrices, denoted as :math:`A_H` for hard constraints and :math:`A_S` for soft constraints, with corresponding fitting data vectors :math:`y_H` and :math:`y_S`.

Constrained Fitting
-------------------

The constrained optimization problem can be expressed as:

.. math::

    \min_{\mathbf{c}, ||\mathbf{A}_\mathrm{H} \mathbf{c} - \mathbf{y}_{\mathrm{H}}||^{2} < \alpha} \left(||\mathbf{y}_{\mathrm{S}} - \mathbf{A}_\mathrm{S}\mathbf{c}||^{2} \right)

Where:

- :math:`c` represents the potential parameters.
- :math:`A_H` and :math:`A_S` are the design matrices for the hard and soft constraints, respectively.
- :math:`y_H` and :math:`y_S` are the target data vectors for the hard and soft constraints.
- :math:`\alpha` is the maximum allowable error for the hard constraint data.

The problem can be solved using an augmented Lagrangian approach, where the Lagrange multiplier :math:`\lambda` enforces the hard constraint. The augmented Lagrangian function is:

.. math::

    \mathcal{L}(\mathbf{c}, \lambda) = ||\mathbf{y}_{\mathrm{S}} - \mathbf{A}_\mathrm{S}\mathbf{c}||^{2} + \lambda (||\mathbf{y}_{\mathrm{H}} - \mathbf{A}_\mathrm{H} \mathbf{c}||^{2} - \alpha)

The optimal potential parameters `c` are found by minimizing this objective function with respect to `c`, with the constraint controlled by the Lagrange multiplier :math:`\lambda`.

The solution to this problem requires finding the appropriate value of :math:`\lambda` such that the fit lies on the boundary of the hard constraint subspace, i.e.,

.. math::

    ||\mathbf{A}_\mathrm{H} \mathbf{c} - \mathbf{y}_{\mathrm{H}}||^{2} = \alpha

The constrained fitting process can be performed iteratively by adjusting the value of `lambda`, starting with a logarithmic search to satisfy the hard constraint, followed by interval bisection to refine the solution.

