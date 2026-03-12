.. _pipeline-timedomain:

======================
Time-Domain Simulation
======================

The MATLAB script ``run_nemoh_timeseries.m`` integrates the Cummins
equation in the time domain using Nemoh hydrodynamic data, Prony
state-space radiation approximation, and user-configurable wave, mooring,
and external-force inputs.

Script: ``run_nemoh_timeseries.m``
----------------------------------

.. rubric:: Location

``SimplePontoon/run_nemoh_timeseries.m``


Configuration Reference
-----------------------

Simulation
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 15 63

   * - Parameter
     - Default
     - Description
   * - ``t_end``
     - 360 s
     - Simulation duration
   * - ``dt_out``
     - 0.01 s
     - Output time step (ODE solver adapts internally)


Wave Input
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Parameter
     - Default
     - Description
   * - ``wave_type``
     - ``'irregular'``
     - ``'regular'`` or ``'irregular'``
   * - ``reg_H``
     - 0.02 m
     - Regular wave height
   * - ``reg_T``
     - 1.0 s
     - Regular wave period
   * - ``irr_Hs``
     - 0.014 m
     - Significant wave height
   * - ``irr_Tp``
     - 3.2 s
     - Peak period
   * - ``irr_gamma``
     - 3.3
     - JONSWAP peak-enhancement factor (1.0 = PM)
   * - ``irr_ncomp``
     - 128
     - Number of spectral components
   * - ``irr_wmin / irr_wmax``
     - 1.0 / 15.0 rad/s
     - Frequency bounds
   * - ``irr_seed``
     - 42
     - Random seed for wave phases

Wave components are discretised uniformly in wavenumber (:math:`k`)
space; see :ref:`theory-spectra` for the rationale.


Environment
^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 15 63

   * - Parameter
     - Default
     - Description
   * - ``water_depth``
     - 1.0 m
     - Water depth (dispersion relation)
   * - ``g_acc``
     - 9.81 m/s²
     - Gravitational acceleration


Prony Fit
^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 22 15 63

   * - Parameter
     - Default
     - Description
   * - ``n_prony``
     - 12
     - Number of exponential terms per DOF pair
   * - ``t_irf_max``
     - 3 s
     - Maximum retardation-function time
   * - ``dt_irf``
     - 0.01 s
     - IRF time step
   * - ``prony_tol``
     - 1e-3
     - Threshold for skipping small DOF pairs


Mooring
^^^^^^^

The mooring stiffness is specified either as a direct 6 × 6 matrix
``K_mooring`` or computed from a ``lines`` struct array via
``compute_mooring_stiffness``.  Three example configurations are
provided in the script:

* **Option A** — 3-line catenary spread at 120°
* **Option B** — 3-line taut wire spread at 120°
* **Option C** — 4-line taut wire spread at 90°


External Force
^^^^^^^^^^^^^^

A user-defined force function handle ``F_ext_func(t, q, v)`` can be
supplied, returning a 6 × 1 vector.  The function receives the current
time, displacement, and velocity, enabling state-dependent forces
(e.g. PTO damping, wind drag).


Active DOFs
^^^^^^^^^^^

The ``active_dofs`` vector (e.g. ``[1, 3, 5]`` for surge, heave, pitch)
selects which degrees of freedom are solved.  All matrices are
sub-indexed accordingly.


Solver Workflow
---------------

1. **Load Nemoh data** — parse ``RadiationCoefficients.tec``,
   ``ExcitationForce.tec``, ``IRF.tec``, ``Kh.dat``, ``Inertia.dat``.
2. **Compute retardation functions** — cosine-transform :math:`B(\omega)`
   to obtain :math:`K_{ij}(t)`.
3. **Prony fitting** — approximate each :math:`K_{ij}(t)` with a sum of
   decaying exponentials (see :ref:`theory-prony`).
4. **Build wave excitation** — synthesise :math:`F_{\text{exc}}(t)` from
   spectral components (see :ref:`theory-spectra`).
5. **Assemble ODE** — form the state vector :math:`[\mathbf{q},\,
   \mathbf{v},\,\mathbf{I}]` and integrate with ``ode45``.
6. **Post-process** — compute forces, statistics, PDFs, exceedance
   curves, and export PDF figures.


Output Figures
--------------

The script generates 11 figures:

1. Added mass :math:`A_{ij}(\omega)`
2. Radiation damping :math:`B_{ij}(\omega)`
3. Excitation force magnitude and phase
4. Retardation functions :math:`K_{ij}(t)` with Prony overlay
5. Wave spectrum and elevation
6. Displacement time histories
7. Velocity time histories
8. Force time histories
9. Probability density functions (PDFs) of steady-state forces
10. Displacement exceedance probability curves
11. Force exceedance probability curves

All figures are exported to ``pdf_outputs/`` as publication-quality PDF
files.


Output Data
-----------

The workspace variable ``timeseries_results`` is saved to a ``.mat``
file containing all time histories, spectral data, Prony coefficients,
and summary statistics including 1 % and 5 % exceedance values.


Transient Handling
------------------

The first ``transient_frac`` (default 20 %) of the time series is
discarded when computing PDFs and exceedance statistics.  This allows
initial transients from zero initial conditions to decay before
statistical analysis.

.. seealso::

   :ref:`theory-cummins` and :ref:`theory-prony` for the mathematical
   formulation, and :ref:`pipeline-sweep` for batch evaluation across
   sea states.
