.. _pipeline-mooring:

=================
Mooring Analysis
=================

The MATLAB function ``compute_mooring_stiffness.m`` calculates the
linearised 6 × 6 mooring stiffness matrix used in the Cummins equation
and parametric sweeps.

Script: ``compute_mooring_stiffness.m``
---------------------------------------

.. rubric:: Location

``SimplePontoon/compute_mooring_stiffness.m``


Function Signature
^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

   [K_moor, F0] = compute_mooring_stiffness(lines, q0)

**Inputs:**

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Argument
     - Description
   * - ``lines``
     - Struct array with one element per mooring line (see below)
   * - ``q0``
     - 6 × 1 equilibrium displacement vector :math:`[x,y,z,\phi,\theta,\psi]^T`

**Outputs:**

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Argument
     - Description
   * - ``K_moor``
     - 6 × 6 linearised mooring stiffness matrix
   * - ``F0``
     - 6 × 1 mooring force/moment at the equilibrium position


Line Structure Fields
^^^^^^^^^^^^^^^^^^^^^

Each element of the ``lines`` struct array must contain:

.. list-table::
   :header-rows: 1
   :widths: 18 12 70

   * - Field
     - Type
     - Description
   * - ``type``
     - string
     - ``'catenary'`` or ``'taut'``
   * - ``anchor``
     - 3 × 1
     - Anchor position :math:`[x, y, z]^T` in global coordinates
   * - ``fairlead``
     - 3 × 1
     - Fairlead attachment in body-fixed coordinates
   * - ``L0``
     - scalar
     - Unstretched line length [m]
   * - ``EA``
     - scalar
     - Axial stiffness [N]
   * - ``w``
     - scalar
     - Weight per unit length in water [N/m] (catenary only)


Example: 3-Line Taut-Wire Configuration
----------------------------------------

.. code-block:: matlab

   depth = 1.0;  R_anch = 1.5;  R_fair = 0.15;  z_fair = -0.01;
   for i = 1:3
       ang = (i-1)*120 * pi/180;
       lines(i).type     = 'taut';
       lines(i).anchor   = [R_anch*cos(ang); R_anch*sin(ang); -depth];
       lines(i).fairlead = [R_fair*cos(ang); R_fair*sin(ang); z_fair];
       lines(i).L0 = 1.6;     % < chord length → pretension
       lines(i).EA = 1000;
   end
   [K_moor, F0] = compute_mooring_stiffness(lines, zeros(6,1));


Example: 3-Line Catenary Configuration
---------------------------------------

.. code-block:: matlab

   depth = 1.0;  R_anch = 4.0;  R_fair = 0.15;
   for i = 1:3
       ang = (i-1)*120 * pi/180;
       lines(i).type     = 'catenary';
       lines(i).anchor   = [R_anch*cos(ang); R_anch*sin(ang); -depth];
       lines(i).fairlead = [R_fair*cos(ang); R_fair*sin(ang); 0];
       lines(i).L0 = 2.5;
       lines(i).w  = 0.2;
       lines(i).EA = 5000;
   end
   [K_moor, F0] = compute_mooring_stiffness(lines, zeros(6,1));


Internal Workflow
-----------------

1. For each perturbation direction :math:`j\in\{1,\dots,6\}`:

   a. Perturb the equilibrium by :math:`+\delta` and :math:`-\delta`.
   b. Transform each fairlead from body-fixed to global coordinates
      using the small-angle rotation matrix.
   c. Solve the line profile (catenary: Newton–Raphson; taut: direct
      elastic formula) to obtain the fairlead force vector.
   d. Map the fairlead force to the body origin as a generalised
      force/moment.
   e. Sum contributions from all lines.

2. Compute central differences: :math:`K_{ij} = -(F_i^+ - F_i^-)/(2\delta)`.


Plotting
--------

An optional ``plot_flag`` activates 3-D visualisation of the mooring
layout, showing anchor positions, fairlead positions, and line profiles
in the global coordinate frame.

.. seealso::

   :ref:`theory-mooring` for the catenary and taut-wire equations, and
   :ref:`pipeline-timedomain` for how the stiffness matrix enters the
   time-domain solver.
