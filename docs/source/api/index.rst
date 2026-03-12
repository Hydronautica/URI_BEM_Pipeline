.. _api:

=============
API Reference
=============

This page provides function signatures, parameter descriptions, and
return values for all scripts in the pipeline.


Python Scripts
--------------

``generate_case.py``
^^^^^^^^^^^^^^^^^^^^

.. rubric:: Module-level script (no function interface)

Generates a complete Nemoh BEM case directory from parametric geometry
definitions.  All configuration is via module-level constants at the top
of the file; see :ref:`pipeline-mesh` for the full parameter table.

**Requires:** Python ≥ 3.6 (standard library only)


``remove_irreg_freq.py``
^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: Module-level script

Reads Nemoh ``.tec`` output files, detects irregular-frequency spikes
using Savitzky–Golay filtering and MAD-based outlier detection, and
writes cleaned coefficient files.

**Requires:** NumPy, SciPy, Matplotlib

.. rubric:: Key internal functions

.. function:: parse_matrix_dat(filepath)

   Parse a whitespace-delimited matrix file (e.g. ``Kh.dat``).

   :param str filepath: Path to the matrix file
   :returns: 2-D NumPy array
   :rtype: numpy.ndarray

.. function:: parse_tec_radiation(filepath, n_dof)

   Parse ``RadiationCoefficients.tec`` into frequency vector, added-mass
   tensor, and damping tensor.

   :param str filepath: Path to the ``.tec`` file
   :param int n_dof: Number of degrees of freedom (6)
   :returns: ``(omega, A, B)`` — frequency array :math:`(N_\omega,)`,
             added mass :math:`(N_\omega, 6, 6)`, damping :math:`(N_\omega, 6, 6)`

.. function:: parse_tec_excitation(filepath, n_dof)

   Parse ``ExcitationForce.tec`` into frequency, magnitude, and phase
   arrays.

   :param str filepath: Path to the ``.tec`` file
   :param int n_dof: Number of degrees of freedom (6)
   :returns: ``(omega, mag, phase)`` — all arrays of shape
             :math:`(N_\omega, 6)`


``nemoh_to_wamit.py``
^^^^^^^^^^^^^^^^^^^^^

.. rubric:: Module-level script

Converts Nemoh-format output to WAMIT ``.1`` / ``.3`` files with
appropriate non-dimensionalisation.

**Requires:** NumPy


MATLAB Scripts
--------------

``run_nemoh_timeseries.m``
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: Main script (no function interface)

Integrates the Cummins equation in the time domain.  Configuration is
via workspace variables defined in the ``USER CONFIGURATION`` section.
See :ref:`pipeline-timedomain` for the full parameter table.

.. rubric:: Key local functions

.. function:: cummins_rhs(t, y, ...)

   Right-hand side of the state-space ODE for ``ode45``.

   :param double t: Current time [s]
   :param double[] y: State vector :math:`[\mathbf{q};\,\mathbf{v};\,\mathbf{I}]`
   :returns: Time derivative :math:`\dot{\mathbf{y}}`
   :rtype: column vector

.. function:: solve_dispersion(omega, d, g)

   Solve the dispersion relation :math:`\omega^2 = gk\tanh(kd)` for
   wavenumber :math:`k` using Newton–Raphson iteration.

   :param double omega: Angular frequency [rad/s]
   :param double d: Water depth [m]
   :param double g: Gravitational acceleration [m/s²]
   :returns: Wavenumber :math:`k` [rad/m]
   :rtype: double


``sweep_sea_states.m``
^^^^^^^^^^^^^^^^^^^^^^

.. rubric:: Main script (no function interface)

Parametric :math:`(H_s, T_p)` sweep; see :ref:`pipeline-sweep`.
Contains the same local functions as ``run_nemoh_timeseries.m``
(``cummins_rhs``, ``solve_dispersion``).


``compute_mooring_stiffness.m``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. function:: compute_mooring_stiffness(lines, q0)

   Compute the linearised 6 × 6 mooring stiffness matrix by central
   differences.

   :param struct[] lines: Array of line structures (fields: ``type``,
          ``anchor``, ``fairlead``, ``L0``, ``EA``, and optionally ``w``)
   :param double[6,1] q0: Equilibrium displacement vector
   :returns: ``[K_moor, F0]`` — stiffness matrix :math:`(6\times6)` and
             equilibrium force vector :math:`(6\times1)`

.. rubric:: Key local functions

.. function:: catenary_solve(FH, FV, w, L0)

   Solve catenary equations for horizontal span :math:`D` and vertical
   span :math:`H` given tensions.

   :param double FH: Horizontal tension [N]
   :param double FV: Vertical tension [N]
   :param double w: Weight per unit length [N/m]
   :param double L0: Unstretched length [m]
   :returns: ``[D, H]`` — spans [m]

.. function:: catenary_inverse(D_target, H_target, w, L0, EA)

   Newton–Raphson solver for fairlead tensions given target span.

   :param double D_target: Horizontal span [m]
   :param double H_target: Vertical span [m]
   :param double w: Weight per unit length [N/m]
   :param double L0: Unstretched length [m]
   :param double EA: Axial stiffness [N]
   :returns: ``[FH, FV]`` — horizontal and vertical tensions [N]

.. function:: taut_force(anchor, fairlead_global, L0, EA)

   Compute taut-wire tension force vector at the fairlead.

   :param double[3,1] anchor: Anchor position (global)
   :param double[3,1] fairlead_global: Fairlead position (global)
   :param double L0: Unstretched length [m]
   :param double EA: Axial stiffness [N]
   :returns: 3 × 1 force vector [N]


``plot_results.m``
^^^^^^^^^^^^^^^^^^

.. rubric:: Utility script

Reads Nemoh output files and generates frequency-domain coefficient
plots (added mass, damping, excitation magnitude/phase).  Useful for
quick inspection before running the time-domain solver.
