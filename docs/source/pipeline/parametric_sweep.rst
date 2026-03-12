.. _pipeline-sweep:

================
Parametric Sweep
================

The sweep script ``sweep_sea_states.m`` evaluates the structural loading
across a grid of sea states, producing contour maps of the 1 %
exceedance force for preliminary design screening.

Script: ``sweep_sea_states.m``
------------------------------

.. rubric:: Location

``SimplePontoon/sweep_sea_states.m``


Purpose
-------

For each combination of significant wave height :math:`H_s` and peak
period :math:`T_p` the script:

1. Synthesises an irregular sea state (JONSWAP or PM spectrum).
2. Integrates the Cummins equation in the time domain.
3. Discards the initial transient.
4. Computes the 1 % exceedance level of two force quantities.


Force Definitions
-----------------

**Hydrodynamic net force** (what the mooring must resist):

.. math::

   F_{\text{hydro,net}} = F_{\text{exc}} + F_{\text{ext}}
   - F_{\text{rad}} - \mathbf{K}_h\,\mathbf{q}

**Total net force** (true unbalanced force = :math:`[M+A_\infty]\,\ddot{\mathbf{q}}`):

.. math::

   F_{\text{net}} = F_{\text{exc}} + F_{\text{ext}}
   - F_{\text{rad}} - (\mathbf{K}_h + \mathbf{K}_{\text{moor}})\,\mathbf{q}


Configuration
-------------

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Parameter
     - Default
     - Description
   * - ``Hs_vec``
     - linspace(0.005, 0.050, 5)
     - Significant wave height grid [m]
   * - ``Tp_vec``
     - linspace(0.6, 3.0, 5)
     - Peak period grid [s]
   * - ``irr_gamma``
     - 1.0
     - JONSWAP peak-enhancement factor
   * - ``irr_ncomp``
     - 128
     - Number of wave components
   * - ``t_end``
     - 1800 s
     - Simulation duration per sea state
   * - ``transient_frac``
     - 0.20
     - Fraction discarded as transient

All other parameters (Prony, ODE, mooring) are configured identically to
``run_nemoh_timeseries.m``.


Performance Optimisation
------------------------

The Nemoh data loading, retardation-function computation, and Prony
fitting are performed **once** before the sweep loop begins.  Only the
wave excitation force synthesis and ODE integration are repeated for each
:math:`(H_s,T_p)` pair, substantially reducing total wall-clock time.


Output
------

**Contour maps** — For each active DOF, two contour plots are generated:

* **Figure A** — :math:`F_{\text{hydro,net}}` 1 % exceedance (jet
  colourmap), showing the load the mooring system must resist.
* **Figure B** — :math:`F_{\text{net}}` 1 % exceedance (parula
  colourmap), showing the true net (unbalanced) force on the body.

Both sets of figures are exported to ``pdf_outputs/`` as PDFs.

**Summary tables** — Printed to the MATLAB console with the 1 %
exceedance value for every :math:`(H_s,T_p)` point and active DOF.


Interpretation
--------------

The contour maps allow rapid identification of the critical sea-state
region for structural design.  Peaks typically appear at combinations of
large :math:`H_s` and :math:`T_p` values near the body's heave natural
period.  Comparing the two force metrics reveals how effectively the
mooring system distributes the hydrodynamic load.

.. seealso::

   :ref:`pipeline-timedomain` for the single sea-state solver and
   :ref:`pipeline-mooring` for mooring configuration.
